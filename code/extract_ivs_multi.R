suppressPackageStartupMessages({library(data.table); library(dplyr); library(ieugwasr)})
jwt <- trimws(readLines("D:/GWAS/TH+SVD/.opengwas_jwt", warn=FALSE, n=1))
thyr <- "D:/GWAS/TH+SVD/thyroid/harmonised"
root <- "D:/GWAS/Migraine"
P_THRESH <- 5e-8
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(),"%H:%M:%S")), ..., "\n", sep="")

# label -> GWAS Catalog accession (single-ancestry EUR where available)
EXPOS <- c(GD="GCST90627744", LT="GCST90627754",
           BNG="GCST90627740", HYPO="GCST90627749", THC="GCST90627763")

read_exposure <- function(acc) {
  f <- file.path(thyr, acc, paste0(acc,".h.tsv.gz"))
  if (!file.exists(f)) { log("MISSING ",f); return(NULL) }
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

allrs <- c()
for (e in names(EXPOS)) {
  d <- read_exposure(EXPOS[e])
  if (is.null(d)) next
  x <- data.frame(rsid=d$rsid, pval=d$pval, stringsAsFactors=FALSE)
  cld <- tryCatch(ieugwasr::ld_clump(x, clump_r2=0.001, clump_kb=10000, pop="EUR", opengwas_jwt=jwt),
                  error=function(err){log("CLUMP ERR ",e,": ",conditionMessage(err));NULL})
  if (is.null(cld)) { log(e,": no clump result, skipping"); next }
  keep <- d[match(cld$rsid, d$rsid), ]
  keep <- keep[!is.na(keep$rsid), ]
  log(e, " (", EXPOS[e], "): ", nrow(keep), " IVs")
  fwrite(keep, file.path(root, paste0("iv_", e, ".tsv")), sep="\t")
  allrs <- union(allrs, keep$rsid)
}
writeLines(allrs, file.path(root, "iv_rsids.txt"))
log("saved union IV rsids: ", length(allrs))
