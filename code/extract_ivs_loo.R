suppressPackageStartupMessages({library(data.table); library(dplyr); library(ieugwasr)})
jwt <- trimws(readLines("D:/GWAS/TH+SVD/.opengwas_jwt", warn=FALSE, n=1))
thyr <- "D:/GWAS/TH+SVD/thyroid/harmonised"
# Consortium-provided leave-one-biobank-out EUR meta-analyses (FinnGen REMAINS in
# GD/HYPO/LT/TC versions; only BNG-LOO also drops FinnGen among others).
EXPOS_LOO <- c(BNG_LOO="GCST90627765", GD_LOO="GCST90627767",
               HYPO_LOO="GCST90627769", LT_LOO="GCST90627771", TC_LOO="GCST90627775")
P_THRESH <- 5e-8
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(),"%H:%M:%S")), ..., "\n", sep="")

read_exposure <- function(acc) {
  f <- file.path(thyr, acc, paste0(acc,".h.tsv.gz"))
  dt <- fread(f, select=c("rsid","chromosome","base_pair_location","effect_allele",
                          "other_allele","beta","standard_error","effect_allele_frequency",
                          "p_value","n"), showProgress=FALSE)
  dt <- dt[!is.na(p_value) & p_value < P_THRESH & !is.na(beta) & !is.na(standard_error) &
             !is.na(rsid) & rsid != ""]
  ea <- toupper(dt$effect_allele); oa <- toupper(dt$other_allele)
  pal <- ((ea=="A"&oa=="T")|(ea=="T"&oa=="A")|(ea=="C"&oa=="G")|(ea=="G"&oa=="C"))
  eaf <- dt$effect_allele_frequency
  dt <- dt[!(pal & !is.na(eaf) & (eaf>0.42 | eaf<0.08))]
  dt[, `:=`(beta=as.numeric(beta), se=as.numeric(standard_error),
            eaf=as.numeric(effect_allele_frequency), pval=as.numeric(p_value),
            n=as.numeric(n), chr=as.character(chromosome))]
  dt[, F := (beta/se)^2]
  dt <- dt[F>10]
  as.data.frame(dt[, c("rsid","chr","base_pair_location","effect_allele",
                       "other_allele","beta","se","eaf","pval","n","F")])
}
for (e in names(EXPOS_LOO)) {
  d <- read_exposure(EXPOS_LOO[e])
  x <- data.frame(rsid=d$rsid, pval=d$pval, stringsAsFactors=FALSE)
  cld <- ieugwasr::ld_clump(x, clump_r2=0.001, clump_kb=10000, pop="EUR", opengwas_jwt=jwt)
  keep <- d[match(cld$rsid, d$rsid), ]
  keep <- keep[!is.na(keep$rsid), ]
  log(e, ": ", nrow(keep), " IVs (from ", nrow(d), " hits)")
  fwrite(keep, file.path("D:/GWAS/Migraine", paste0("iv_", e, ".tsv")), sep="\t")
}
log("done LOO IV extraction")
