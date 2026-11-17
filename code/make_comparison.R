suppressPackageStartupMessages({library(data.table); library(ggplot2)})
resdir <- "D:/GWAS/Migraine/results"
log <- function(...) cat(sprintf("[%s] ",format(Sys.time(),"%H:%M:%S")),...,"\n",sep="")

ivw <- function(f) {
  d <- fread(f)
  d <- d[method=="Inverse variance weighted"]
  d[, trait := sub("_to_.*","",key)]
  d[, trait := sub("_LOO$","",trait)]
  d[, outcome := sub(".*_to_","",key)]
  d[, `:=`(OR=exp(b), lo=exp(b-1.96*se), hi=exp(b+1.96*se))]
  d[, .(trait, outcome, nsnp, b, se, pval, OR, lo, hi)]
}
A <- ivw(file.path(resdir,"mr_results.csv"))
L <- ivw(file.path(resdir,"mr_results_loo.csv"))
setnames(L, c("nsnp","b","se","pval","OR","lo","hi"),
            c("nsnp.loo","b.loo","se.loo","pval.loo","OR.loo","lo.loo","hi.loo"))
cmp <- merge(A[, .(trait, outcome, nsnp, b, se, pval, OR, lo, hi)],
             L[, .(trait, outcome, nsnp.loo, b.loo, se.loo, pval.loo, OR.loo, lo.loo, hi.loo)],
             by=c("trait","outcome"), all=TRUE)
cmp <- cmp[order(trait, outcome)]
fwrite(cmp, file.path(resdir,"mr_compare.csv"))
log("comparison rows: ", nrow(cmp))

# Long format for forest plot (paired Main vs LOO)
long <- rbind(
  cmp[, .(trait, outcome, version="Main (GBMI full)", OR, lo, hi, nsnp)],
  cmp[, .(trait, outcome, version="LOO (exclude biobank)", OR=OR.loo, lo=lo.loo, hi=hi.loo, nsnp=nsnp.loo)]
)[!is.na(OR)]
long[, outcome := factor(outcome, levels=c("Migraine","MigraineAura","MigraineNoAura"),
                         labels=c("Migraine (all)","Migraine w/ aura","Migraine w/o aura"))]
long[, trait := factor(trait, levels=c("GD","HYPO","LT","BNG","TC"),
                       labels=c("GD","Hypothyroid","LT","Benign goiter","Thyroid cancer"))]
p <- ggplot(long, aes(x=OR, y=trait, xmin=lo, xmax=hi, color=version)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey50") +
  geom_pointrange(position=position_dodge(width=0.45), size=0.7, linewidth=0.7) +
  scale_x_continuous(limits=c(0.94,1.10), breaks=seq(0.94,1.10,0.02)) +
  scale_color_manual(values=c("Main (GBMI full)"="#1f77b4","LOO (exclude biobank)"="#d62728")) +
  facet_wrap(~outcome, ncol=3) +
  labs(x="Odds ratio (IVW, 95% CI)", y="Thyroid exposure",
       title="Thyroid disease -> migraine: Main vs leave-one-biobank-out exposure",
       subtitle="Note: LOO versions exclude CCPM/KCPS (FinnGen retained); FinnGen overlap NOT removed",
       color="Exposure") +
  theme_bw(base_size=11) + theme(legend.position="bottom", strip.background=element_rect(fill="#eee"))
ggsave(file.path(resdir,"forest_main_vs_loo.png"), p, width=11, height=4.6, dpi=180)
log("forest saved: forest_main_vs_loo.png")

# Console summary
cmp[, `:=`(A_OR=round(OR,3), A_p=signif(pval,3), L_OR=round(OR.loo,3), L_p=signif(pval.loo,3))]
print(cmp[, .(trait, outcome, nsnp, A_OR, A_p, nsnp.loo, L_OR, L_p)])
