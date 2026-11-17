## mrlap_ldcheck.R
## 检验 MRlap 默认 pruning (500kb 距离, 无 LD clumping) 选出的 IV 是否因
## LD 相关而低估标准误 -> 把同样的 IV 用主 MR 的参数 (r2=0.001, 10Mb, EUR)
## 重新 clumping, 再比较 IVW 效应与 p 值。
suppressPackageStartupMessages({library(data.table); library(ieugwasr)})
ROOT <- "D:/GWAS/Migraine"
e <- "GD"; o <- "Migraine"
jwt <- trimws(readLines("D:/GWAS/TH+SVD/.opengwas_jwt", warn=FALSE, n=1))

iv <- fread(file.path(ROOT, "mrlap_run/default", paste0(e,"_to_",o), "ivs.csv"))$SNP
g  <- fread(file.path(ROOT, "mrlap_input", paste0(e, ".tsv")))
g  <- g[rsid %in% iv]
g[, pval := 2*pnorm(-abs(beta/se))]
cat(sprintf("MRlap 默认 IV: %d ; 在输入表中找到 %d\n", length(iv), nrow(g)))

cl <- tryCatch(
  ieugwasr::ld_clump(data.frame(rsid=g$rsid, pval=g$pval, stringsAsFactors=FALSE),
                     clump_r2=0.001, clump_kb=10000, pop="EUR", opengwas_jwt=jwt),
  error=function(x){cat("ld_clump ERROR:", conditionMessage(x), "\n"); NULL})
if (is.null(cl)) quit(status=1)
keep <- cl$rsid
cat(sprintf("LD clumping (r2=0.001, 10Mb, EUR) 后剩余: %d / %d\n", length(keep), nrow(g)))

res <- readRDS(file.path(ROOT, "mrlap_run/default", paste0(e,"_to_",o), "result.rds"))
h <- as.data.table(res$harmonised_mr_data)

ivw_std <- function(d) {
  w <- 1/d$std_SE.out^2
  b <- sum(w*d$std_beta.out*d$std_beta.exp)/sum(w*d$std_beta.exp^2)
  se <- sqrt(1/sum(w*d$std_beta.exp^2))
  c(b=b, se=se, p=2*pnorm(-abs(b/se)))
}
cat("\n--- MRlap 标准化尺度 IVW ---\n")
a0 <- ivw_std(h); cat(sprintf("  全部 %d IV : b=%.5f se=%.5f p=%.4g\n", nrow(h), a0[1], a0[2], a0[3]))
a1 <- ivw_std(h[rsid %in% keep]); cat(sprintf("  clump 后 %d IV: b=%.5f se=%.5f p=%.4g\n",
                                              sum(h$rsid %in% keep), a1[1], a1[2], a1[3]))

## 反推 log OR 尺度
se_e <- g[, setNames(se, rsid)]; N_e <- g[, setNames(N, rsid)]
ou <- fread(file.path(ROOT, "mrlap_input", paste0(o, ".tsv")))[rsid %in% iv]
se_o <- ou[, setNames(se, rsid)]; N_o <- ou[, setNames(N, rsid)]
d <- h[, .(rsid, std_beta.exp, std_beta.out)]
d[, `:=`(be = std_beta.exp * se_e[rsid] * sqrt(N_e[rsid]),
         bo = std_beta.out * se_o[rsid] * sqrt(N_o[rsid]),
         seo = se_o[rsid])]
ivw_logor <- function(dd) {
  w <- 1/dd$seo^2
  b <- sum(w*dd$bo*dd$be)/sum(w*dd$be^2)
  se <- sqrt(1/sum(w*dd$be^2))
  c(b=b, se=se, p=2*pnorm(-abs(b/se)), OR=exp(b))
}
cat("\n--- 反推 log OR 尺度, 标准 IVW ---\n")
b0 <- ivw_logor(d); cat(sprintf("  全部 %d IV : b=%.5f se=%.5f p=%.4g OR=%.4f\n", nrow(d), b0[1], b0[2], b0[3], b0[4]))
b1 <- ivw_logor(d[rsid %in% keep]); cat(sprintf("  clump 后 %d IV: b=%.5f se=%.5f p=%.4g OR=%.4f\n",
                                                sum(d$rsid %in% keep), b1[1], b1[2], b1[3], b1[4]))

## 主 MR 的 IV 有多少落在 MRlap 输入(hm3 过滤)里
main <- fread(file.path(ROOT, "results", paste0("harmonised_",e,"_to_",o,".csv")))[mr_keep==TRUE]
cat(sprintf("\n主 MR IV: %d ; 在 mrlap 输入(hm3 集)中: %d ; 与 MRlap 默认 IV 交集: %d\n",
            nrow(main), sum(main$SNP %in% g$rsid), length(intersect(main$SNP, iv))))
fwrite(data.table(SNP=keep), file.path(ROOT, "results", "mrlap_GD_Migraine_clumped_IVs.csv"))
