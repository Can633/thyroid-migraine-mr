# extract_ivs_hypo.R
# ---------------------------------------------------------------------------
# 数据修复：GCST90627749 (HYPO EUR) 原始下载被截断（352,321,536 / 477,246,985 字节），
# 导致旧 iv_HYPO.tsv 只覆盖 chr1-12（248 个 IV），而同病 LOO 版覆盖 chr1-22（337 个）。
# 补齐文件后，用与 extract_ivs.R 完全相同的规则重新提取 HYPO 的工具变量，
# 并重建 iv_rsids.txt（所有 iv_*.tsv 的并集，含 LOO 版），供结局 SNP 提取使用。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table); library(dplyr); library(ieugwasr)})
jwt  <- trimws(readLines("D:/GWAS/TH+SVD/.opengwas_jwt", warn = FALSE, n = 1))
thyr <- "D:/GWAS/TH+SVD/thyroid/harmonised"
root <- "D:/GWAS/Migraine"
ACC <- "GCST90627749"
P_THRESH <- 5e-8
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

read_exposure <- function(acc) {
  f <- file.path(thyr, acc, paste0(acc, ".h.tsv.gz"))
  dt <- fread(f, select = c("rsid","chromosome","base_pair_location","effect_allele",
                            "other_allele","beta","standard_error","effect_allele_frequency",
                            "p_value","n"), showProgress = FALSE)
  log("read ", acc, ": ", nrow(dt), " rows")
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
  dt <- dt[F > 10]
  log("genome-wide significant, F>10: ", nrow(dt),
      " ; chromosomes: ", paste(sort(unique(dt$chr)), collapse=","))
  as.data.frame(dt[, c("rsid","chr","base_pair_location","effect_allele",
                       "other_allele","beta","se","eaf","pval","n","F")])
}

d <- read_exposure(ACC)
x <- data.frame(rsid = d$rsid, pval = d$pval, stringsAsFactors = FALSE)
cld <- ieugwasr::ld_clump(x, clump_r2 = 0.001, clump_kb = 10000,
                          pop = "EUR", opengwas_jwt = jwt)
keep <- d[match(cld$rsid, d$rsid), ]
keep <- keep[!is.na(keep$rsid), ]
log("HYPO clumped IVs: ", nrow(keep),
    " ; chromosomes: ", paste(sort(unique(keep$chr)), collapse=","))
fwrite(keep, file.path(root, "iv_HYPO.tsv"), sep = "\t")

# 重建并集（包含主版 + LOO 版，避免遗漏）
fs <- list.files(root, pattern = "^iv_.*\\.tsv$", full.names = TRUE)
rs <- unique(unlist(lapply(fs, function(f) fread(f, select="rsid")[[1]])))
writeLines(rs, file.path(root, "iv_rsids.txt"))
log("saved IV rsids: ", length(rs))
