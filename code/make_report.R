suppressPackageStartupMessages({library(data.table); library(ggplot2); library(dplyr)})
root<-"D:/GWAS/Migraine"; res<-file.path(root,"results")

tab<-fread(file.path(res,"mr_results.csv"))
het<-fread(file.path(res,"mr_het.csv"))
ple<-fread(file.path(res,"mr_ple.csv"))

# keep primary + sensitivity methods, compute OR and CI
tab<-tab[method %in% c("Inverse variance weighted","MR Egger","Weighted median","Simple mode","Weighted mode")]
if (!"hi_ci" %in% names(tab)) { tab[, hi_ci:=b+1.96*se]; tab[, lo_ci:=b-1.96*se] }
tab[, `:=`(OR=exp(b), OR_lo=exp(lo_ci), OR_hi=exp(hi_ci))]
# pretty labels
tab[exposure=="GD", exp_label:="Graves disease (GD)"]
tab[exposure=="LT", exp_label:="Lymphocytic thyroiditis (LT)"]
tab[exposure=="BNG", exp_label:="Benign nodular goiter (BNG)"]
tab[exposure=="HYPO", exp_label:="Hypothyroidism (HYPO)"]
tab[exposure=="TC", exp_label:="Thyroid cancer (TC)"]
tab[outcome=="Migraine", out_label:="Migraine (overall)"]
tab[outcome=="MigraineAura", out_label:="Migraine with aura"]
tab[outcome=="MigraineNoAura", out_label:="Migraine without aura"]
tab[, combo:=paste(exp_label,"->",out_label)]
tab[, method_lab:=factor(method, levels=c("Inverse variance weighted","Weighted median","Weighted mode","Simple mode","MR Egger"),
                          labels=c("IVW (primary)","Weighted median","Weighted mode","Simple mode","MR Egger"))]
tab[, sig:=ifelse(pval<0.05,"p<0.05","p>=0.05")]

fwrite(tab, file.path(res,"mr_summary_or.csv"))

# forest plot: IVW + weighted median per combination
plot_dat <- tab[method %in% c("Inverse variance weighted","Weighted median")]
ivw<-tab[method=="Inverse variance weighted"]
p<-ggplot(plot_dat, aes(x=OR, y=combo, color=sig)) +
  geom_point(size=3.2) +
  geom_errorbarh(aes(xmin=OR_lo, xmax=OR_hi), height=0.22, linewidth=0.9) +
  geom_vline(xintercept=1, linetype="dashed", color="grey40") +
  scale_x_log10(breaks=c(0.8,1,1.25,1.6,2,2.5)) +
  facet_wrap(~method_lab, nrow=1) +
  labs(title="Thyroid disease vs Migraine: IVW (primary) and sensitivity estimates",
       subtitle="GBMI/VTB 2026 GD & LT (EUR) to FinnGen R12 migraine endpoints (EUR); OR per SD of genetically predicted disease liability",
       x="Odds ratio (95% CI)", y=NULL, color="Significance") +
  theme_bw(base_size=11) +
  theme(strip.text=element_text(size=10), legend.position="bottom",
        axis.text.y=element_text(size=9))
ggsave(file.path(res,"forest_mr.png"), p, width=13, height=9.0, dpi=180)

# full 5-method plot
p2<-ggplot(tab, aes(x=OR, y=method_lab, color=sig)) +
  geom_point(size=2.6) + geom_errorbarh(aes(xmin=OR_lo,xmax=OR_hi),height=0.18,linewidth=0.8) +
  geom_vline(xintercept=1,linetype="dashed",color="grey40") +
  scale_x_log10() + facet_wrap(~combo, nrow=2, scales="free_x") +
  labs(title="All MR methods by exposure-outcome pair", x="Odds ratio (95% CI)", y=NULL) +
  theme_bw(base_size=10) + theme(legend.position="bottom", strip.text=element_text(size=8))
ggsave(file.path(res,"forest_all_methods.png"), p2, width=13, height=8, dpi=160)

cat("Saved forest plots and summary.\n")
print(as.data.frame(ivw[,.(combo,nsnp,OR=round(OR,3),CI95=paste0(round(OR_lo,3),"-",round(OR_hi,3)),pval=signif(pval,3))]))
