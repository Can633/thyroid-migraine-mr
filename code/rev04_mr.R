# =====================================================================
# rev04_mr.R —— 反向 MR：偏头痛 → 甲状腺疾病
# =====================================================================
# 暴露：FinnGen R12 G6_MIGRAINE / _WITH_AURA / _NO_AURA（IV 见 rev/iv_mig_*.tsv）
# 结局：GD / HYPO / LT / BNG / TC（GCST90627740/744/749/754/763）
# 方法：IVW 主分析 + MR-Egger / 加权中位数 / simple mode / weighted mode
# 敏感性：Cochran Q、Egger 截距、留一法；并做「剔除 MHC IV」的对照
# 输出：results/rev_mr_results.csv, rev_mr_het.csv, rev_mr_ple.csv,
#       rev_mr_summary_or.csv, rev_mr_loo.csv
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(TwoSampleMR)
})
root <- "D:/GWAS/Migraine"
rev  <- file.path(root, "rev")
res  <- file.path(root, "results"); dir.create(res, showWarnings = FALSE)
log  <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

TRAITS <- c("Migraine", "MigraineAura", "MigraineNoAura")
THYR   <- c(GD = "GD", HYPO = "HYPO", LT = "LT", BNG = "BNG", TC = "TC")
MHC_LO <- 25e6; MHC_HI <- 34e6
N_MIG  <- c(Migraine = 395060, MigraineAura = 379299, MigraineNoAura = 376762)  # FinnGen R13

METHODS <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
             "mr_simple_mode", "mr_weighted_mode")

all_res <- list(); all_het <- list(); all_ple <- list(); all_loo <- list()

run_one <- function(tr, th, drop_mhc = FALSE) {
  iv <- fread(file.path(rev, sprintf("iv_mig_%s.tsv", tr)))
  if (drop_mhc) iv <- iv[!(chr == 6 & pos >= MHC_LO & pos <= MHC_HI)]
  if (nrow(iv) < 3) return(NULL)

  iv[, `:=`(pval = pval, eaf = eaf)]
  ex <- data.frame(SNP = iv$rsid, beta = iv$beta, se = iv$se,
                   effect_allele = iv$alt, other_allele = iv$ref,
                   eaf = iv$eaf, pval = iv$pval,
                   samplesize = N_MIG[[tr]], stringsAsFactors = FALSE)
  ex <- format_data(ex, type = "exposure", snp_col = "SNP", beta_col = "beta",
                    se_col = "se", effect_allele_col = "effect_allele",
                    other_allele_col = "other_allele", eaf_col = "eaf",
                    pval_col = "pval", samplesize_col = "samplesize")

  ou <- fread(file.path(rev, sprintf("out_%s.tsv", th)))
  ou <- ou[!is.na(beta) & !is.na(se) & se > 0]
  ou <- ou[!duplicated(rsid)]
  ns <- median(ou$n, na.rm = TRUE)
  ou <- data.frame(SNP = ou$rsid, beta = ou$beta, se = ou$se,
                   effect_allele = ou$ea, other_allele = ou$oa,
                   eaf = ou$eaf, pval = ou$p,
                   samplesize = ns, stringsAsFactors = FALSE)
  ou <- format_data(ou, type = "outcome", snp_col = "SNP", beta_col = "beta",
                    se_col = "se", effect_allele_col = "effect_allele",
                    other_allele_col = "other_allele", eaf_col = "eaf",
                    pval_col = "pval", samplesize_col = "samplesize")

  d <- tryCatch(harmonise_data(ex, ou, action = 2), error = function(e) NULL)
  if (is.null(d) || sum(d$mr_keep) < 3) {
    log("  ", tr, "->", th, ": harmonise failed or <3 SNP"); return(NULL)
  }
  d$id.exposure <- tr; d$id.outcome <- th

  r <- mr(d, method_list = METHODS)
  r <- generate_odds_ratios(r)
  r$direction <- "reverse"
  r$drop_mhc <- drop_mhc
  r$nsnp_harmonised <- sum(d$mr_keep)

  het <- tryCatch(mr_heterogeneity(d), error = function(e) NULL)
  if (!is.null(het)) { het$direction <- "reverse"; het$drop_mhc <- drop_mhc }
  ple <- tryCatch(mr_pleiotropy_test(d), error = function(e) NULL)
  if (!is.null(ple)) { ple$direction <- "reverse"; ple$drop_mhc <- drop_mhc }
  loo <- tryCatch(mr_leaveoneout(d), error = function(e) NULL)
  if (!is.null(loo)) { loo$direction <- "reverse"; loo$drop_mhc <- drop_mhc }

  ivw <- r[r$method == "Inverse variance weighted", ]
  log(sprintf("  %-14s -> %-5s  nSNP=%2d  OR=%.4f (%.4f-%.4f)  p=%.4g%s",
              tr, th, ivw$nsnp, ivw$or, ivw$or_lci95, ivw$or_uci95, ivw$pval,
              ifelse(drop_mhc, "  [no MHC]", "")))
  list(r = r, het = het, ple = ple, loo = loo)
}

for (tr in TRAITS) {
  for (th in THYR) {
    for (dm in c(FALSE, TRUE)) {
      o <- run_one(tr, th, dm)
      if (is.null(o)) next
      all_res[[length(all_res) + 1]] <- o$r
      if (!is.null(o$het)) all_het[[length(all_het) + 1]] <- o$het
      if (!is.null(o$ple)) all_ple[[length(all_ple) + 1]] <- o$ple
      if (!is.null(o$loo)) all_loo[[length(all_loo) + 1]] <- o$loo
    }
  }
}

R <- rbindlist(all_res, fill = TRUE)
H <- rbindlist(all_het, fill = TRUE)
P <- rbindlist(all_ple, fill = TRUE)
L <- rbindlist(all_loo, fill = TRUE)

# 主表（仅主分析，不含 no-MHC 对照）
main <- R[drop_mhc == FALSE]
main[, `:=`(bonferroni = 0.05 / 15,
            sig_bonf = pval < 0.05 / 15)]
fwrite(R,    file.path(res, "rev_mr_results.csv"),  sep = ",")
fwrite(main, file.path(res, "rev_mr_main.csv"),     sep = ",")
fwrite(H,    file.path(res, "rev_mr_het.csv"),      sep = ",")
fwrite(P,    file.path(res, "rev_mr_ple.csv"),      sep = ",")
fwrite(L,    file.path(res, "rev_mr_loo.csv"),      sep = ",")

cat("\n=== 反向 MR 汇总（IVW，含与不含 MHC 两版）===\n")
sm <- R[method == "Inverse variance weighted",
        .(exposure = id.exposure, outcome = id.outcome, nsnp,
          OR = round(or, 4), lo = round(or_lci95, 4), hi = round(or_uci95, 4),
          p = signif(pval, 3), drop_mhc)]
print(sm[order(exposure, outcome, drop_mhc)])

cat("\n=== 5 种方法一致性检查（主分析，p<0.05 的方法数）===\n")
cnt <- main[, .(n_sig_method = sum(pval < 0.05), n_method = .N),
            by = .(exposure = id.exposure, outcome = id.outcome)]
print(cnt[order(-n_sig_method)])
