# -*- coding: utf-8 -*-
# 15 组合的 MR-PRESSO 全局检验 + 离群值 + 校正后估计；并做 radial MR（IVW radial）
suppressMessages({
  library(MRPRESSO); library(RadialMR); library(data.table)
})
root <- "D:/GWAS/Migraine"
res  <- file.path(root, "results")
files <- list.files(res, pattern = "^harmonised_.*_to_.*\\.csv$", full.names = TRUE)
files <- files[!grepl("_LOO_", basename(files))]
stopifnot(length(files) == 15)

out <- list()
for (f in files) {
  nm <- sub("^harmonised_(.*)_to_(.*)\\.csv$", "\\1|\\2", basename(f))
  d  <- fread(f)
  d  <- d[toupper(as.character(mr_keep)) %in% c("TRUE","1","T")]
  d  <- d[is.finite(beta.exposure) & is.finite(beta.outcome) &
          is.finite(se.exposure)   & is.finite(se.outcome)]
  if (nrow(d) < 4) {
    out[[length(out)+1]] <- data.table(combo = nm, nsnp = nrow(d), note = "too few SNPs")
    next
  }
  # ---- MR-PRESSO ----
  pr <- tryCatch(
    mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
              SdOutcome = "se.outcome", SdExposure = "se.exposure",
              data = as.data.frame(d), OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
              NbDistribution = 1000, SignifThreshold = 0.05, seed = 20261116),
    error = function(e) NULL)

  g_p <- NA_real_; n_out <- NA_integer_; b_adj <- NA_real_; se_adj <- NA_real_
  if (!is.null(pr)) {
    g_p <- pr$`MR-PRESSO results`$`Global Test`$Pvalue
    or  <- pr$`MR-PRESSO results`$`Outlier Test`
    if (!is.null(or)) {
      n_out <- sum(or$Pvalue < 0.05, na.rm = TRUE)
      b_adj <- unname(pr$`Main MR results`$`Causal Estimate`[2])
      se_adj<- unname(pr$`Main MR results`$`Sd`[2])
    } else {
      n_out <- 0L
    }
  }
  # ---- radial MR（IVW radial，Q 检验与离群值）----
  rad_q <- NA_real_; rad_nout <- NA_integer_
  rr <- tryCatch({
    r <- RadialMR::format_radial(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, d$SNP)
    ivw_radial(r, alpha = 0.05, weights = 1)
  }, error = function(e) NULL)
  if (!is.null(rr)) {
    rad_q <- rr$Qstatistic; rad_nout <- sum(rr$outliers$Pvalue < 0.05, na.rm = TRUE)
  }
  out[[length(out)+1]] <- data.table(
    combo = nm, nsnp = nrow(d),
    presso_global_p = g_p, presso_outliers = n_out,
    presso_b = b_adj, presso_se = se_adj,
    presso_OR = ifelse(is.na(b_adj), NA_real_, exp(b_adj)),
    radial_Q = rad_q, radial_outliers = rad_nout)
  cat(sprintf("%-26s k=%3d  PRESSO p=%s  out=%s  radial out=%s\n", nm, nrow(d),
              ifelse(is.na(g_p), "NA", formatC(g_p, format="g", digits=3)),
              ifelse(is.na(n_out), "NA", n_out), ifelse(is.na(rad_nout), "NA", rad_nout)))
  flush.console()
}
ans <- rbindlist(out, fill = TRUE)
fwrite(ans, file.path(res, "mr_presso_radial.csv"))
cat("\n写出 results/mr_presso_radial.csv\n")
print(ans)
