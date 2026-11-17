# -*- coding: utf-8 -*-
# Colocalization (coloc.abf) of thyroid exposures vs migraine outcomes,
# locus by locus, using the long input table produced by coloc_prepare.py.
suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

ROOT <- "D:/GWAS/Migraine"
INDIR <- file.path(ROOT, "coloc_input")
OUTDIR <- file.path(ROOT, "coloc_results")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

f <- file.path(INDIR, "coloc_long.tsv")
if (!file.exists(f)) stop("missing ", f)
d <- fread(f)
cat("rows:", nrow(d), " blocks:", nrow(unique(d[, .(locus, exposure, outcome)])), "\n")

# fixed per-dataset sample sizes (same as LDSC / MR pipeline)
N_EXP <- c(GD = 1866840, HYPO = 1771600, LT = 1405560, BNG = 1792740, TC = 2091280)
N_OUT <- c(Migraine = 395060, MigraineAura = 379299, MigraineNoAura = 376762)

args <- commandArgs(trailingOnly = TRUE)
p12 <- if (length(args) >= 1) as.numeric(args[[1]]) else 1e-5
tag <- if (length(args) >= 2) args[[2]] else "main"
cat("prior p12 =", p12, " tag =", tag, "\n")

blocks <- unique(d[, .(locus, exposure, outcome)])

# Keep only (locus, exposure) pairs where that exposure actually contributes an
# IV to the locus. Otherwise we would be running coloc on regions that were
# never selected by this exposure's signal, which makes PP.H4 uninterpretable.
loci <- fread(file.path(INDIR, "loci.tsv"))
# NOTE: must group by locus here. Without `by = locus`, unlist() of the
# comma-split exposures recycles across rows and silently duplicates pairs.
valid <- loci[, .(exposure = unlist(strsplit(exposures, ","))), by = locus][
  , .(locus, exposure)]
valid <- unique(valid)
blocks <- merge(blocks, valid, by = c("locus", "exposure"))
blocks <- unique(blocks)
cat("blocks after keeping only originating exposures:", nrow(blocks), "\n")
res <- vector("list", nrow(blocks))
for (i in seq_len(nrow(blocks))) {
  b <- blocks[i]
  sub <- d[locus == b$locus & exposure == b$exposure & outcome == b$outcome]
  sub <- sub[!duplicated(snp)]
  if (nrow(sub) < 30) {
    res[[i]] <- data.table(locus = b$locus, exposure = b$exposure,
                           outcome = b$outcome, nsnps = nrow(sub),
                           PP.H0 = NA_real_, PP.H1 = NA_real_, PP.H2 = NA_real_,
                           PP.H3 = NA_real_, PP.H4 = NA_real_, status = "skipped_small")
    next
  }
  out <- tryCatch({
    d1 <- list(snp = as.character(sub$snp), beta = as.numeric(sub$beta_x),
               varbeta = as.numeric(sub$se_x)^2, MAF = as.numeric(sub$maf),
               N = unname(N_EXP[[b$exposure]]), type = "cc")
    d2 <- list(snp = as.character(sub$snp), beta = as.numeric(sub$beta_y),
               varbeta = as.numeric(sub$se_y)^2, MAF = as.numeric(sub$maf),
               N = unname(N_OUT[[b$outcome]]), type = "cc")
    # signal strength in the region: coloc PP.H0 is uninformative when neither
    # trait has a signal here, so these columns decide which blocks to trust
    zx <- abs(as.numeric(sub$beta_x) / as.numeric(sub$se_x))
    zy <- abs(as.numeric(sub$beta_y) / as.numeric(sub$se_y))
    r <- suppressWarnings(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = p12))
    s <- r$summary
    data.table(locus = b$locus, exposure = b$exposure, outcome = b$outcome,
               nsnps = as.integer(s[["nsnps"]]),
               min_z_exp = max(zx), min_z_out = max(zy),
               PP.H0 = as.numeric(s[["PP.H0.abf"]]), PP.H1 = as.numeric(s[["PP.H1.abf"]]),
               PP.H2 = as.numeric(s[["PP.H2.abf"]]), PP.H3 = as.numeric(s[["PP.H3.abf"]]),
               PP.H4 = as.numeric(s[["PP.H4.abf"]]),
               top_snp = as.character(r$results$snp[which.max(r$results$SNP.PP.H4)]),
               top_snp_pp = as.numeric(max(r$results$SNP.PP.H4)),
               status = "ok")
  }, error = function(e) {
    data.table(locus = b$locus, exposure = b$exposure, outcome = b$outcome,
               nsnps = nrow(sub), min_z_exp = NA_real_, min_z_out = NA_real_,
               PP.H0 = NA_real_, PP.H1 = NA_real_, PP.H2 = NA_real_,
               PP.H3 = NA_real_, PP.H4 = NA_real_,
               top_snp = NA_character_, top_snp_pp = NA_real_,
               status = paste("error:", conditionMessage(e)))
  })
  res[[i]] <- out
}
out <- rbindlist(res, fill = TRUE)
out[, `:=`(prior_p12 = p12)]
fwrite(out, file.path(OUTDIR, sprintf("coloc_%s.tsv", tag)), sep = "\t")
cat("done blocks:", nrow(out), " ok:", sum(out$status == "ok"),
    " PP.H4>0.8:", sum(out$status == "ok" & out$PP.H4 > 0.8),
    " PP.H4>0.5:", sum(out$status == "ok" & out$PP.H4 > 0.5), "\n")
