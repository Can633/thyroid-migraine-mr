## mrlap_diagnose.R
## 诊断 MRlap 观测效应与主 MR IVW 差异的来源:
##   (A) IV 集不同 (500kb 距离 pruning vs LD clumping r2=0.001/10Mb)
##   (B) 权重不同 (MRlap: 1/std_SE.out^2 = N_out 常数; 主MR: 1/se_out^2)
##   (C) 尺度不同 (MRlap: std_beta = beta/(se*sqrt(N)), 二分类 log OR 下非常数)
suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/GWAS/Migraine"
e <- "GD"; o <- "Migraine"

res <- readRDS(file.path(ROOT, "mrlap_run/default", paste0(e, "_to_", o), "result.rds"))
h <- as.data.table(res$harmonised_mr_data)
cat("harmonised_mr_data 列:\n"); print(names(h))
cat("行数:", nrow(h), "\n\n")

## (1) 复现 MRlap 的 IVW (标准化尺度, 权重 = 1/std_SE.out^2)
w <- 1 / h$std_SE.out^2
a_std <- sum(w * h$std_beta.out * h$std_beta.exp) / sum(w * h$std_beta.exp^2)
se_std <- sqrt(1 / sum(w * h$std_beta.exp^2))
cat(sprintf("复现 MRlap IVW: %.6f (se %.6f, p %.4g)  [MRlap 报 %.6f]\n",
            a_std, se_std, 2*pnorm(-abs(a_std/se_std)), res$MRcorrection$observed_effect))

## (2) 补回 se / N, 并由 std_beta 反推【已按等位基因对齐】的 log OR
##     beta_aligned = std_beta * se * sqrt(N)  (MRlap 的 std_beta 已含方向翻转)
ex <- fread(file.path(ROOT, "mrlap_input", paste0(e, ".tsv")))
ou <- fread(file.path(ROOT, "mrlap_input", paste0(o, ".tsv")))
ex <- ex[, .(rsid, se_e=se, N_e=N)]
ou <- ou[, .(rsid, se_o=se, N_o=N)]
d <- merge(h[, .(rsid, std_beta.exp, std_beta.out, std_SE.out)], ex, by="rsid")
d <- merge(d, ou, by="rsid")
d[, `:=`(be = std_beta.exp * se_e * sqrt(N_e),
         bo = std_beta.out * se_o * sqrt(N_o))]

## (3) 在 MRlap 的 IV 上做两种 IVW
w_std <- 1 / d$se_o^2
b_ivw <- sum(w_std * d$bo * d$be) / sum(w_std * d$be^2)
se_ivw <- sqrt(1/sum(w_std * d$be^2))
b_eq  <- sum(d$bo * d$be) / sum(d$be^2)
cat(sprintf("同一批 MRlap IV, 反推 log OR 尺度:\n"))
cat(sprintf("  标准 IVW (w=1/se_out^2): b=%.6f se=%.6f p=%.4g  OR=%.4f\n",
            b_ivw, se_ivw, 2*pnorm(-abs(b_ivw/se_ivw)), exp(b_ivw)))
cat(sprintf("  等权重回归            : b=%.6f  OR=%.4f\n", b_eq, exp(b_eq)))

## (4) 尺度因子 c = 1/(se*sqrt(N)) 的变异 (MRlap 假定近似常数)
c_e <- 1/(d$se_e*sqrt(d$N_e)); c_o <- 1/(d$se_o*sqrt(d$N_o))
cat(sprintf("\n尺度因子 c_exp: median=%.5f CV=%.3f range=[%.5f, %.5f]\n",
            median(c_e), sd(c_e)/mean(c_e), min(c_e), max(c_e)))
cat(sprintf("尺度因子 c_out: median=%.5f CV=%.3f range=[%.5f, %.5f]\n",
            median(c_o), sd(c_o)/mean(c_o), min(c_o), max(c_o)))
cat(sprintf("理论 alpha_std / alpha_logOR = c_out/c_exp (median) = %.4f\n",
            median(c_o)/median(c_e)))
cat(sprintf("实测 alpha_std / alpha_IVW_logOR = %.4f\n\n", a_std/b_ivw))

## (5) 主 MR 自己的 IVW 与 IV 集重叠
main <- fread(file.path(ROOT, "results", paste0("harmonised_", e, "_to_", o, ".csv")))
main <- main[mr_keep == TRUE]
w_m <- 1/main$se.outcome^2
b_m <- sum(w_m*main$beta.outcome*main$beta.exposure)/sum(w_m*main$beta.exposure^2)
se_m <- sqrt(1/sum(w_m*main$beta.exposure^2))
cat(sprintf("主 MR IV 数: %d | MRlap IV 数: %d | 交集: %d\n",
            nrow(main), nrow(d), length(intersect(main$SNP, d$rsid))))
cat(sprintf("主 MR 自己的 IVW: b=%.6f se=%.6f p=%.4g OR=%.4f\n",
            b_m, se_m, 2*pnorm(-abs(b_m/se_m)), exp(b_m)))
cat(sprintf("IV 集效应 (标准 IVW 尺度): %.6f -> %.6f  (差 %.6f)\n", b_m, b_ivw, b_ivw-b_m))
