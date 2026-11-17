# =====================================================================
# GD / LT (GBMI 2026) -> Migraine & subtypes (FinnGen R12) : two-sample MR
# =====================================================================
# Exposure : GCST90627744 (GD, EUR), GCST90627754 (LT, EUR); GBMI/VTB
#            Nat Genet 2026 PMID 41644669, GRCh38, harmonised GWAS-SSF
# Outcome  : finngen_R13_G6_MIGRAINE / _WITH_AURA / _NO_AURA (GRCh38)
# IVs      : OpenGWAS 1000G EUR LD clump r2<0.001/10Mb, F>10, pal-filter
#           Steiger directionality; IVW primary + Egger/WM/mode/singlesnp
# Output   : D:/GWAS/Migraine/results/
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(TwoSampleMR); library(ieugwasr)
})
jwt <- trimws(readLines("D:/GWAS/TH+SVD/.opengwas_jwt", warn = FALSE, n = 1))
root <- "D:/GWAS/Migraine"
resdir <- file.path(root, "results"); dir.create(resdir, showWarnings = FALSE, recursive = TRUE)
thyr <- "D:/GWAS/TH+SVD/thyroid/harmonised"
EXPOS <- c(GD = "GCST90627744", LT = "GCST90627754")
OUTCOMES <- c(
  Migraine       = "G6_MIGRAINE",
  MigraineAura   = "G6_MIGRAINE_WITH_AURA",
  MigraineNoAura = "G6_MIGRAINE_NO_AURA")
P_THRESH <- 5e-8; CLUMP_R2 <- 0.001; CLUMP_KB <- 10000
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
log("TwoSampleMR ", as.character(packageVersion("TwoSampleMR")), " | JWT len ", nchar(jwt))

# ---- 1. exposure instruments (stream gz, keep significant + F + pal filters) ----
read_exposure <- function(acc) {
  f <- file.path(thyr, acc, paste0(acc, ".h.tsv.gz"))
  log("Exposure ", acc)
  dt <- fread(f, select = c("rsid","chromosome","base_pair_location",
                            "effect_allele","other_allele","beta",
                            "standard_error","effect_allele_frequency","p_value","n"),
              showProgress = FALSE)
  dt <- dt[!is.na(p_value) & p_value < P_THRESH & !is.na(beta) & !is.na(standard_error) &
             !is.na(rsid) & rsid != ""]
  ea <- toupper(dt$effect_allele); oa <- toupper(dt$other_allele)
  pal <- ((ea=="A"&oa=="T")|(ea=="T"&oa=="A")|(ea=="C"&oa=="G")|(ea=="G"&oa=="C"))
  eaf <- dt$effect_allele_frequency
  dt <- dt[!(pal & !is.na(eaf) & (eaf > 0.42 | eaf < 0.08))]
  dt[, `:=`(beta = as.numeric(beta), se = as.numeric(standard_error),
            eaf = as.numeric(effect_allele_frequency), pval = as.numeric(p_value),
            n = as.numeric(n), chr = as.character(chromosome))]
  dt[, F := (beta/se)^2]
  dt <- dt[F > 10]
  log("  raw IVs after p/F/pal filters: ", nrow(dt))
  as.data.frame(dt[, c("rsid","chr","base_pair_location","effect_allele",
                       "other_allele","beta","se","eaf","pval","n","F")])
}

# ---- 2. LD clump via OpenGWAS (1000G EUR) ----
# ieugwasr::ld_clump takes data.frame(rsid, pval) and returns data.frame(rsid,...)
clump_ivs <- function(dat, label) {
  x <- data.frame(rsid = dat$rsid, pval = dat$pval, stringsAsFactors = FALSE)
  x <- x[!is.na(x$rsid) & x$rsid != "" & !is.na(x$pval), ]
  cld <- tryCatch(ieugwasr::ld_clump(x, clump_r2 = CLUMP_R2, clump_kb = CLUMP_KB,
                                     pop = "EUR", opengwas_jwt = jwt),
                  error = function(e) { log("  clump ERROR: ", conditionMessage(e)); NULL })
  if (is.null(cld) || nrow(cld) == 0) return(NULL)
  keep <- dat[match(cld$rsid, dat$rsid), ]
  log("  ", label, ": clumped ", nrow(x), " -> ", nrow(keep), " IVs")
  keep
}

# ---- 3. read outcome for a union of rsIDs (single pass per outcome) ----
read_outcome <- function(code, want) {
  f <- file.path(root, paste0("finngen_R13_", code, ".gz"))
  if (!file.exists(f)) stop("Missing: ", f)
  log("Outcome ", code, ": scanning for ", length(want), " SNPs")
  res <- fread(f, select = c("rsids","ref","alt","pval","beta","sebeta"),
               header = TRUE, showProgress = FALSE)
  res <- res[rsids %in% want]
  setnames(res, c("rsids","pval","beta","sebeta"),
           c("rsid","pval_out","beta_out","se_out"))
  res <- unique(res, by = "rsid")
  log("  matched ", nrow(res))
  as.data.frame(res)
}

# ---- 4. build exposure/outcome formatted dat then harmonise ----
run_mr_battery <- function(exp_df, oc_df, exp_label, out_label) {
  m <- merge(exp_df, oc_df, by = "rsid")
  if (nrow(m) < 3) return(NULL)
  # TwoSampleMR expects .exposure / .outcome suffixed columns
  exp_dat <- data.frame(
    SNP = m$rsid,
    beta.exposure = as.numeric(m$beta), se.exposure = as.numeric(m$se),
    pval.exposure = as.numeric(m$pval),
    effect_allele.exposure = toupper(m$effect_allele),
    other_allele.exposure  = toupper(m$other_allele),
    eaf.exposure = as.numeric(m$eaf),
    id.exposure = exp_label, exposure = exp_label,
    stringsAsFactors = FALSE)
  out_dat <- data.frame(
    SNP = m$rsid,
    beta.outcome = as.numeric(m$beta_out), se.outcome = as.numeric(m$se_out),
    pval.outcome = as.numeric(m$pval_out),
    effect_allele.outcome = toupper(m$alt),
    other_allele.outcome  = toupper(m$ref),
    id.outcome = out_label, outcome = out_label,
    stringsAsFactors = FALSE)
  keep <- !is.na(exp_dat$beta.exposure) & !is.na(out_dat$beta.outcome) &
          !is.na(exp_dat$se.exposure) & !is.na(out_dat$se.outcome)
  exp_dat <- exp_dat[keep, ]; out_dat <- out_dat[keep, ]
  hd <- harmonise_data(exp_dat, out_dat, action = 2)
  hd <- hd[!is.na(hd$beta.exposure) & !is.na(hd$beta.outcome) &
             !is.na(hd$se.exposure) & !is.na(hd$se.outcome), ]
  if (nrow(hd) < 3) return(NULL)
  ste <- tryCatch(steiger_filtering(hd), error = function(e) NULL)
  if (!is.null(ste) && "steiger_dir" %in% names(ste)) {
    ok <- sum(ste$steiger_dir, na.rm = TRUE)
    log("  Steiger kept ", ok, "/", nrow(hd))
    if (ok >= 3) hd <- ste[ste$steiger_dir, ] else log("  Steiger disabled (too few)")
  }
  log("  ", exp_label, "->", out_label, ": ", nrow(hd), " SNPs in MR")
  list(nsnp = nrow(hd),
       result = mr(hd),
       het = tryCatch(mr_heterogeneity(hd), error = function(e) NULL),
       ple = tryCatch(mr_pleiotropy_test(hd), error = function(e) NULL),
       singles = tryCatch(mr_singlesnp(hd), error = function(e) NULL),
       hd = hd)
}

# ---- main: clump all exposures first, then one outcome pass each ----
exps <- list(); ivsets <- list()
for (e in names(EXPOS)) {
  exps[[e]] <- read_exposure(EXPOS[e])
  ivsets[[e]] <- clump_ivs(exps[[e]], e)
}

allres <- list()
for (o in names(OUTCOMES)) {
  code <- OUTCOMES[o]
  want <- unique(unlist(lapply(ivsets, function(x) if (!is.null(x)) x$rsid else character())))
  if (length(want) == 0) { log("No IVs for ", o); next }
  oc <- read_outcome(code, want)
  for (e in names(EXPOS)) {
    ivs <- ivsets[[e]]
    if (is.null(ivs)) next
    key <- paste0(e, "_to_", o)
    r <- tryCatch(run_mr_battery(ivs, oc, e, o),
                  error = function(err) { log("ERROR ", key, ": ", conditionMessage(err)); NULL })
    if (!is.null(r)) allres[[key]] <- r
  }
}

# ---- save ----
if (length(allres)) {
  tab <- do.call(rbind, lapply(names(allres), function(k) {
    r <- allres[[k]]; cbind(data.frame(exposure = r$exposure, outcome = r$outcome,
                                       nsnp = r$nsnp), r$result)
  }))
  # fix labels (run_mr_battery has no labels in result)
  tab$exposure <- gsub("_to_.*", "", rep(names(allres), each = 1))
  # simpler: rebuild with key
  tab <- do.call(rbind, lapply(names(allres), function(k) {
    r <- allres[[k]]; parts <- strsplit(k, "_to_")[[1]]
    cbind(data.frame(exposure = parts[1], outcome = parts[2], nsnp = r$nsnp), r$result)
  }))
  fwrite(tab, file.path(resdir, "mr_results.csv"))
  for (nm in c("het","ple")) {
    tl <- do.call(rbind, lapply(names(allres), function(k) {
      z <- allres[[k]][[nm]]; if (is.null(z)) return(NULL)
      cbind(data.frame(key = k), z)
    }))
    if (!is.null(tl)) fwrite(tl, file.path(resdir, paste0("mr_", nm, ".csv")))
  }
  for (k in names(allres)) fwrite(allres[[k]]$hd, file.path(resdir, paste0("harmonised_", k, ".csv")))
  log("Saved to ", resdir)
  print(as.data.frame(tab[, c("exposure","outcome","nsnp","method",
                              "b","se","pval","lo_ci","hi_ci")], stringsAsFactors = FALSE))
} else log("NO RESULTS")
