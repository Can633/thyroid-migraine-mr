#!/usr/bin/env Rscript
# 跨种族敏感性分析: 甲状腺病(GBMI, 欧洲) -> 偏头痛(Genes&Health, 南亚 GCST90727115)
# 对 5 个暴露读取已按位置对齐的 harmonized 文件, 手动计算 MR 方法
suppressPackageStartupMessages(library(data.table))

expos <- c("GD","LT","BNG","HYPO","TC")
labels <- c(GD="Graves disease (GD)", LT="Lymphocytic thyroiditis (LT)",
            BNG="Benign thyroid goitre (BNG)", HYPO="Hypothyroidism (HYPO)",
            TC="Thyroid cancer (TC)")
out_label <- "Migraine (Genes&Health, South Asian, n=44,026)"

res <- list()
for (E in expos) {
  f <- paste0("D:/GWAS/Migraine/harm_small_", E, ".csv")
  if (!file.exists(f)) { cat("MISSING", f, "\n"); next }
  d <- fread(f)
  if (nrow(d) == 0) { cat(E, "no instruments matched\n"); next }
  b_exp <- d$b_exp; se_exp <- d$se_exp
  b_out <- d$b_out; se_out <- d$se_out
  # ratio estimate + delta-method SE
  r   <- b_out / b_exp
  se_r <- sqrt((se_out / b_exp)^2 + (b_out * se_exp / b_exp^2)^2)
  w   <- 1 / se_r^2

  # IVW
  b_ivw <- sum(r * w) / sum(w)
  se_ivw <- sqrt(1 / sum(w))
  p_ivw <- 2 * pnorm(-abs(b_ivw / se_ivw))

  # MR-Egger (weights = 1/se_out^2)
  we <- 1 / se_out^2
  fit <- lm(b_out ~ b_exp, weights = we)
  b_eg <- coef(fit)[2]; se_eg <- sqrt(vcov(fit)[2,2]); p_eg <- 2*pnorm(-abs(b_eg/se_eg))
  int_eg <- coef(fit)[1]; se_int <- sqrt(vcov(fit)[1,1]); p_int <- 2*pnorm(-abs(int_eg/se_int))

  # Weighted median (simple)
  ord <- order(w, decreasing = TRUE)
  cum <- cumsum(w[ord]) / sum(w)
  med_idx <- ord[which(cum >= 0.5)[1]]
  b_wm <- r[med_idx]; se_wm <- se_r[med_idx]
  p_wm <- 2 * pnorm(-abs(b_wm / se_wm))

  # Cochran Q (IVW)
  Q <- sum(w * (r - b_ivw)^2); df <- nrow(d) - 1
  p_Q <- pchisq(Q, df, lower.tail = FALSE)

  res[[E]] <- data.table(
    exposure = labels[[E]], code = E, outcome = out_label,
    nsnp = nrow(d),
    ivw_b = b_ivw, ivw_se = se_ivw, ivw_p = p_ivw,
    ivw_OR = exp(b_ivw), ivw_lo = exp(b_ivw - 1.96*se_ivw), ivw_hi = exp(b_ivw + 1.96*se_ivw),
    egger_b = b_eg, egger_se = se_eg, egger_p = p_eg, egger_int = int_eg, egger_int_p = p_int,
    wmed_OR = exp(b_wm), wmed_p = p_wm,
    Q = Q, Q_p = p_Q
  )
  cat(sprintf("%-32s n=%2d  IVW OR=%.3f (%.3f-%.3f) p=%.3g | Egger p_int=%.3g | WMed OR=%.3f p=%.3g | Q p=%.3g\n",
              labels[[E]], nrow(d), exp(b_ivw), exp(b_ivw-1.96*se_ivw), exp(b_ivw+1.96*se_ivw), p_ivw,
              p_int, exp(b_wm), p_wm, p_Q))
}

tab <- rbindlist(res)
fwrite(tab, "D:/GWAS/Migraine/results/mr_small_summary.csv")
cat("\nSAVED -> D:/GWAS/Migraine/results/mr_small_summary.csv\n")
