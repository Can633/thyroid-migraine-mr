# =====================================================================
# rev02_clump.R —— 反向 MR：偏头痛 IV 的 LD clumping
# =====================================================================
# 主路径：ieugwasr::ld_clump（OpenGWAS 1000G EUR，r2<0.001，10,000 kb），
#         与前向 MR 完全一致，保证两个方向的 IV 选择口径可比。
# 退化路径：若 API 不可用，用 10 Mb 距离贪心（按 p 升序，10 Mb 内只留最显著）。
#         理由：clump_kb=10000 且 r2<0.001 时，窗口内几乎任意两两 SNP 的 r2 都
#         会超过 0.001，实际效果近似"每 10 Mb 取一个 top SNP"，故该退化是合理
#         近似；会在输出中标明所用路径，报告中需如实披露。
# 输出：rev/iv_mig_<trait>.tsv（含 F 与 clump_method 列）
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(ieugwasr)
})
root <- "D:/GWAS/Migraine"
rev  <- file.path(root, "rev"); dir.create(rev, showWarnings = FALSE)
jwt  <- tryCatch(trimws(readLines("D:/GWAS/TH+SVD/.opengwas_jwt", warn = FALSE, n = 1)),
                 error = function(e) "")
log  <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

CLUMP_R2 <- 0.001; CLUMP_KB <- 10000

distance_prune <- function(d, kb = CLUMP_KB) {
  d <- d[order(d$pval), ]
  keep <- rep(TRUE, nrow(d))
  for (i in seq_len(nrow(d))) {
    if (!keep[i]) next
    idx <- (i + 1):nrow(d)
    if (length(idx) == 0 || is.na(idx[1])) break
    near <- abs(d$pos[idx] - d$pos[i]) < kb * 1000 & d$chr[idx] == d$chr[i]
    if (any(near)) keep[idx[near]] <- FALSE
  }
  d[keep, ]
}

traits <- c("Migraine", "MigraineAura", "MigraineNoAura")
summary_rows <- list()

for (t in traits) {
  f <- file.path(rev, sprintf("iv_mig_%s.raw.tsv", t))
  d <- fread(f)
  d[, F := (beta / se)^2]
  log(t, ": candidates = ", nrow(d))

  method <- "opengwas_ld_clump"
  keep <- NULL
  if (nzchar(jwt)) {
    x <- data.frame(rsid = d$rsid, pval = d$pval, stringsAsFactors = FALSE)
    x <- x[!duplicated(x$rsid) & !is.na(x$pval), ]
    keep <- tryCatch(
      ieugwasr::ld_clump(x, clump_r2 = CLUMP_R2, clump_kb = CLUMP_KB,
                         pop = "EUR", opengwas_jwt = jwt),
      error = function(e) { log("  ld_clump ERROR: ", conditionMessage(e)); NULL })
  }
  if (is.null(keep) || nrow(keep) == 0) {
    log("  -> fallback: 10 Mb distance pruning")
    keep <- distance_prune(d)
    method <- "distance_10Mb"
  } else {
    keep <- d[match(keep$rsid, d$rsid), ]
    keep <- keep[!is.na(rsid)]
    method <- "opengwas_ld_clump"
  }

  keep[, clump_method := method]
  keep <- keep[order(chr, pos)]
  fwrite(keep, file.path(rev, sprintf("iv_mig_%s.tsv", t)), sep = "\t")
  log("  IVs kept: ", nrow(keep), " (", method, ")")
  summary_rows[[t]] <- data.frame(trait = t, n_candidate = nrow(d),
                                  n_iv = nrow(keep), method = method)
}

sm <- rbindlist(summary_rows)
fwrite(sm, file.path(rev, "iv_clump_summary.tsv"), sep = "\t")
print(sm)
