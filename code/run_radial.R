suppressMessages({library(RadialMR); library(data.table)})
res <- "D:/GWAS/Migraine/results"
files <- list.files(res, pattern = "^harmonised_.*_to_.*\\.csv$", full.names = TRUE)
files <- files[!grepl("_LOO_", basename(files))]
out <- list()
for (f in files) {
  nm <- sub("^harmonised_(.*)_to_(.*)\\.csv$", "\\1|\\2", basename(f))
  d <- fread(f)
  d <- d[toupper(as.character(mr_keep)) %in% c("TRUE","1","T")]
  r <- tryCatch({
    rr <- RadialMR::format_radial(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, d$SNP)
    RadialMR::ivw_radial(rr, alpha = 0.05, weights = 1)
  }, error = function(e) NULL)
  if (is.null(r)) { out[[length(out)+1]] <- data.table(combo=nm, k=nrow(d), Q=NA_real_, Q_p=NA_real_, n_out=NA_integer_); next }
  out[[length(out)+1]] <- data.table(combo=nm, k=nrow(d), Q=r$Qstatistic, Q_p=r$Qpval,
                                     n_out=sum(r$outliers$Pvalue < 0.05, na.rm=TRUE))
}
a <- rbindlist(out)
fwrite(a, file.path(res, "radial_mr.csv"))
print(a)
cat("\n合计离群变异:", sum(a$n_out, na.rm=TRUE), "\n")
