suppressPackageStartupMessages({library(data.table); library(TwoSampleMR)})
root <- "D:/GWAS/Migraine"; resdir <- file.path(root,"results_ukb")
dir.create(resdir, showWarnings=FALSE, recursive=TRUE)
log <- function(...) cat(sprintf("[%s] ",format(Sys.time(),"%H:%M:%S")),...,"\n",sep="")

EXPOS <- c(BNG="iv_BNG.tsv", GD="iv_GD.tsv", HYPO="iv_HYPO.tsv",
          LT="iv_LT.tsv", TC="iv_TC.tsv")
OUTS <- c(MigraineUKB="outcome_ukb.tsv")
# FinnGen R12 sample sizes (cases + controls), verified via official API 2026-11-12
OUT_N <- c(MigraineUKB=458440)  # UK Biobank WGS, GCST90473326

run_battery <- function(expf, outf, elabel, olabel) {
  exp <- fread(file.path(root,expf))
  oc  <- fread(file.path(root,outf))
  setnames(exp, "base_pair_location", "pos", skip_absent=TRUE)
  m <- merge(exp, oc, by="rsid")
  log(elabel,"->",olabel,": merged ",nrow(m)," SNPs")
  if (nrow(m)<3) return(NULL)
  expdat <- data.frame(
    SNP=m$rsid, beta.exposure=as.numeric(m$beta), se.exposure=as.numeric(m$se),
    pval.exposure=as.numeric(m$pval),
    effect_allele.exposure=toupper(m$effect_allele),
    other_allele.exposure=toupper(m$other_allele),
    eaf.exposure=as.numeric(m$eaf), samplesize.exposure=as.numeric(m$n),
    id.exposure=elabel, exposure=elabel, stringsAsFactors=FALSE)
  outdat <- data.frame(
    SNP=m$rsid, beta.outcome=as.numeric(m$beta_out), se.outcome=as.numeric(m$se_out),
    pval.outcome=as.numeric(m$pval_out),
    effect_allele.outcome=toupper(m$alt), other_allele.outcome=toupper(m$ref),
    eaf.outcome=as.numeric(m$eaf_out), samplesize.outcome=as.numeric(OUT_N[olabel]),
    id.outcome=olabel, outcome=olabel, stringsAsFactors=FALSE)
  hd <- harmonise_data(expdat, outdat, action=2)
  log("  harmonised -> ",nrow(hd))
  if (nrow(hd)<3) return(NULL)
  ste <- tryCatch(steiger_filtering(hd), error=function(e) NULL)
  if (!is.null(ste) && "steiger_dir" %in% names(ste)) {
    ok <- sum(ste$steiger_dir, na.rm=TRUE)
    log("  Steiger ",ok,"/",nrow(hd))
    if (ok>=3) hd <- ste[ste$steiger_dir,]
  } else {
    log("  Steiger filtering unavailable (no directionality column)")
  }
  list(nsnp=nrow(hd), result=mr(hd),
       het=tryCatch(mr_heterogeneity(hd),error=function(e) NULL),
       ple=tryCatch(mr_pleiotropy_test(hd),error=function(e) NULL),
       singles=tryCatch(mr_singlesnp(hd),error=function(e) NULL), hd=hd)
}

allres <- list()
for (on in names(OUTS)) {
  of <- file.path(root, OUTS[on])
  if (!file.exists(of) || file.size(of)<100) { log("skip ",on,": file missing/empty"); next }
  for (en in names(EXPOS)) {
    k <- paste0(en,"_to_",on)
    r <- tryCatch(run_battery(EXPOS[en], OUTS[on], en, on),
                  error=function(err){log("ERROR ",k,": ",conditionMessage(err));NULL})
    if (!is.null(r)) allres[[k]] <- r
  }
}
if (length(allres)) {
  tab <- do.call(rbind, lapply(names(allres), function(k){
    r<-allres[[k]]; z<-r$result
    z$key<-k  # avoid duplicate exposure/outcome columns from mr() output
    z}))
  fwrite(tab, file.path(resdir,"mr_results.csv"))
  for (nm in c("het","ple")) {
    z <- do.call(rbind, lapply(names(allres), function(k){
      v<-allres[[k]][[nm]]; if(is.null(v)) return(NULL); cbind(data.frame(key=k),v)}))
    if (!is.null(z)) fwrite(z, file.path(resdir,paste0("mr_",nm,".csv")))
  }
  for (k in names(allres)) fwrite(allres[[k]]$hd, file.path(resdir,paste0("harmonised_",k,".csv")))
  log("SAVED ",length(allres)," combinations")
  print(as.data.frame(tab[,c("key","nsnp","method","b","se","pval")], stringsAsFactors=FALSE))
} else log("NO RESULTS")
