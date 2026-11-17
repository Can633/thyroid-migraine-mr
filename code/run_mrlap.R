## run_mrlap.R
## ---------------------------------------------------------------------------
## 用已装好的 MRlap 0.0.3.3 + GenomicSEM 0.0.5 对 5 个甲状腺暴露 x 3 个
## FinnGen R12 偏头痛结局做「观测 IVW-MR + 样本重叠/弱工具/Winner's curse 校正」。
##
## 用法:
##   Rscript run_mrlap.R all            # 串行跑全部 15 组合 (default 模式)
##   Rscript run_mrlap.R GD_to_Migraine
##   Rscript run_mrlap.R GD_to_Migraine userIV
##
## mode = default : MRlap 自带的 IV 选择 (p<5e-8, 距离 500kb pruning, 反向剔除)
## mode = userIV  : do_pruning=FALSE, 使用主 MR (TwoSampleMR) 完全相同的 IV 集,
##                  便于直接比较「同一批工具变量下校正前后的变化」
##
## 输出 (每组合独立目录, 避免 GenomicSEM::munge 临时文件互相覆盖):
##   mrlap_run/<mode>/<combo>/result.rds       完整 MRlap 结果
##   mrlap_run/<mode>/<combo>/row.csv          单行汇总
##   mrlap_run/<mode>/<combo>/ivs.csv         使用的工具变量
##   mrlap_run/<mode>/<combo>/run.log         MRlap 详细日志
##   mrlap_run/<mode>/<combo>/*_munge.log     munge 日志 (save_logfiles=TRUE)
##   mrlap_run/<mode>/<combo>/*_ldsc.log      LDSC 日志
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(MRlap)
  library(data.table)
})

source("D:/GWAS/MRlap_ref/mrlap_config.R")

ROOT   <- "D:/GWAS/Migraine"
RUNROOT <- file.path(ROOT, "mrlap_run")
EXPS <- c("GD", "HYPO", "LT", "BNG", "TC")
OUTS <- c("Migraine", "MigraineAura", "MigraineNoAura")

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args) >= 1) args[1] else "all"
mode   <- if (length(args) >= 2) args[2] else "default"
stopifnot(mode %in% c("default", "userIV"))

## userIV 需要主 MR 的完整 IV 集（多数不在 hm3 面板内），用补齐版输入；
## 补齐只影响 run_MR，不影响 munge/LDSC。
INDIR <- if (mode == "userIV") file.path(ROOT, "mrlap_input_full") else file.path(ROOT, "mrlap_input")

combos <- unlist(lapply(EXPS, function(e) paste0(e, "_to_", OUTS)), use.names = FALSE)
todo <- if (target == "all") combos else target
stopifnot(all(todo %in% combos))

## 近似尺度换算因子: alpha_std ~ (h_X/h_Y)*alpha_logOR,  h = median(se*sqrt(N))
sf <- fread(file.path(ROOT, "mrlap_input", "scale_factors.csv"))
h  <- setNames(sf$se_sqrtN_median, sf$trait)

dir.create(file.path(RUNROOT, mode), recursive = TRUE, showWarnings = FALSE)

for (combo in todo) {
  parts <- strsplit(combo, "_to_")[[1]]
  e <- parts[1]; o <- parts[2]
  wd <- file.path(RUNROOT, mode, combo)
  dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  cat("\n============================================================\n")
  cat(sprintf("[%s | %s] %s -> %s\n", format(Sys.time(), "%H:%M:%S"), mode, e, o))
  cat("============================================================\n")

  old <- getwd()
  setwd(wd)
  row <- NULL
  t0 <- proc.time()[3]
  err <- ""
  logcon <- file("run.log", open = "wt")
  sink(logcon, type = "output")
  res <- tryCatch({
    exp_path <- normalizePath(file.path(INDIR, paste0(e, ".tsv")))
    out_path <- normalizePath(file.path(INDIR, paste0(o, ".tsv")))

    extra <- list()
    if (mode == "userIV") {
      hf <- file.path(ROOT, "results", paste0("harmonised_", e, "_to_", o, ".csv"))
      hd <- fread(hf)
      ivs <- if ("mr_keep" %in% names(hd)) hd[mr_keep == TRUE, SNP] else hd$SNP
      ivs <- unique(as.character(ivs[!is.na(ivs)]))
      cat("userIV: ", length(ivs), " SNPs from ", basename(hf), "\n", sep = "")
      extra <- list(do_pruning = FALSE, user_SNPsToKeep = ivs)
    }

    do.call(MRlap, c(list(
      exposure = exp_path,
      exposure_name = paste0(e, "_exp"),
      outcome = out_path,
      outcome_name = paste0(o, "_out"),
      ld = MRlap_ld_dir,
      hm3 = MRlap_hm3,
      MR_threshold = 5e-8,
      MR_pruning_dist = 500,
      MR_pruning_LD = 0,
      MR_reverse = 1e-3,
      save_logfiles = TRUE,
      verbose = TRUE
    ), extra))
  }, error = function(x) {
    err <<- conditionMessage(x)
    cat("ERROR:", conditionMessage(x), "\n")
    NULL
  })
  sink(type = "output")
  close(logcon)
  runtime <- as.numeric(proc.time()[3] - t0)

  if (!is.null(res)) {
    saveRDS(res, "result.rds")
    mc <- res$MRcorrection
    ld <- res$LDSC
    ga <- res$GeneticArchitecture
    fwrite(data.table(SNP = mc$IVs), "ivs.csv")

    ## 近似换算 (仅量级对照)
    ratio <- unname(h[o] / h[e])
    approx_b  <- mc$observed_effect  * ratio
    approx_cb <- mc$corrected_effect * ratio

    row <- data.table(
      exposure = e, outcome = o, mode = mode, combo = combo,
      m_IVs = mc$m_IVs,
      observed_effect = mc$observed_effect,
      observed_effect_se = mc$observed_effect_se,
      observed_effect_p = mc$observed_effect_p,
      corrected_effect = mc$corrected_effect,
      corrected_effect_se = mc$corrected_effect_se,
      corrected_effect_p = mc$corrected_effect_p,
      test_difference = mc$test_difference,
      p_difference = mc$p_difference,
      egger_b = mc$egger_b, egger_se = mc$egger_se,
      egger_intercept_p = mc$egger_intercept_p,
      h2_exp = ld$h2_exp, h2_exp_se = ld$h2_exp_se, int_exp = ld$int_exp,
      h2_out = ld$h2_out, h2_out_se = ld$h2_out_se, int_out = ld$int_out,
      gcov = ld$gcov, gcov_se = ld$gcov_se, rg = ld$rg,
      int_crosstrait = ld$int_crosstrait,
      int_crosstrait_se = ld$int_crosstrait_se,
      polygenicity = ga$polygenicity,
      perSNP_heritability = ga$perSNP_heritability,
      scale_ratio = ratio,
      approx_b_logOR = approx_b, approx_OR = exp(approx_b),
      approx_b_corrected_logOR = approx_cb, approx_OR_corrected = exp(approx_cb),
      runtime_s = round(runtime, 1), error = ""
    )
    cat(sprintf("  IVs=%d | obs=%.5f (p=%.4g) -> corrected=%.5f (p=%.4g) | int_cross=%.5f | %.0fs\n",
                mc$m_IVs, mc$observed_effect, mc$observed_effect_p,
                mc$corrected_effect, mc$corrected_effect_p,
                ld$int_crosstrait, runtime))
  } else {
    row <- data.table(exposure = e, outcome = o, mode = mode, combo = combo,
                      m_IVs = NA_integer_, observed_effect = NA_real_,
                      corrected_effect = NA_real_, runtime_s = round(runtime, 1),
                      error = err)
    cat("  FAILED:", err, "\n")
  }
  fwrite(row, "row.csv")
  setwd(old)
}

cat("\n===== run_mrlap.R DONE =====\n")
