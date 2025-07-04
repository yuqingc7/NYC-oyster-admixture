setwd("/workdir/yc2644/CV_NYC_lcWGS/results_ngsadmix/")

target="HHSV85FI"
pop <- read.table(paste0("ALL_ds_Re0.2_",target,".info"), as.is=T) %>% 
  dplyr::rename(p=V1, ind=V2)

K=2

all <- read.table(paste0("ngsadmix_K",K,"_run1_ALL_ds_Re0.2_",
                         target,"_ALLsnplist_LDpruned_maf0.05_pval1e-6_pind0.86_cv30_noinver.qopt"))
lcWGS_q <- cbind(pop, all) %>% 
  dplyr::select(p,ind,V1) %>% 
  dplyr::rename(all_lcWGS=V1) %>% 
  mutate(ind = str_replace_all(ind, "003a", "003")) %>% 
  mutate(ind = str_replace_all(ind, "008a", "008"))
quantile(subset(lcWGS_q, p == "HH13a")$all_lcWGS, 0.90) # 0.006789585
dim(subset(lcWGS_q,p=="SV0512a" & all_lcWGS > 0.0068)) #23
dim(subset(lcWGS_q,p=="CT5785" & all_lcWGS > 0.0068)) #29
table(lcWGS_q$p)
# CT5785  FI1012   HH13a SV0512a 
# 30      17      25      24 

mean(subset(lcWGS_q,p=="SV0512a" & all_lcWGS > 0.0068)$all_lcWGS)
mean(subset(lcWGS_q,p=="CT5785" & all_lcWGS > 0.0068)$all_lcWGS)

all_ddRAD <- read.table(paste0("ngsadmix_K",K,"_run1_ALL_ds_Re0.2_",
                         target,"_ALLsnplist_LDpruned_maf0.05_pval1e-6_pind0.86_cv30_noinver.ddRAD.qopt"))
lcWGS_ddRAD_q <- cbind(pop, all_ddRAD) %>% 
  dplyr::select(p,ind,V2) %>% 
  dplyr::rename(all_lcWGS_in_ddRAD=V2) %>% 
  mutate(ind = str_replace_all(ind, "003a", "003")) %>% 
  mutate(ind = str_replace_all(ind, "008a", "008"))

all_ddRAD_f <- read.table(paste0("ngsadmix_K",K,"_run1_ALL_ds_Re0.2_",
                               target,"_ALLsnplist_LDpruned_maf0.05_pval1e-6_pind0.86_cv30_noinver.ddRAD_sites.qopt"))
lcWGS_ddRAD_f_q <- cbind(pop, all_ddRAD_f) %>% 
  dplyr::select(p,ind,V1) %>% 
  dplyr::rename(all_lcWGS_in_ddRAD_f=V1) %>% 
  mutate(ind = str_replace_all(ind, "003a", "003")) %>% 
  mutate(ind = str_replace_all(ind, "008a", "008"))

all_coding <- read.table(paste0("ngsadmix_K",K,"_run1_ALL_ds_Re0.2_",
                                target,"_ALLsnplist_LDpruned_maf0.05_pval1e-6_pind0.86_cv30_noinver.codingSNP36868.qopt"))
lcWGS_coding_q <- cbind(pop, all_coding) %>% 
  dplyr::select(p,ind,V1) %>% 
  dplyr::rename(all_lcWGS_in_CDS=V1) %>% 
  mutate(ind = str_replace_all(ind, "003a", "003")) %>% 
  mutate(ind = str_replace_all(ind, "008a", "008"))

all_noncoding <- read.table(paste0("ngsadmix_K",K,"_run1_ALL_ds_Re0.2_",
                                target,"_ALLsnplist_LDpruned_maf0.05_pval1e-6_pind0.86_cv30_noinver.noncodingSNP36868.qopt"))
lcWGS_noncoding_q <- cbind(pop, all_noncoding) %>% 
  dplyr::select(p,ind,V2) %>% 
  dplyr::rename(all_lcWGS_not_in_CDS=V2) %>% 
  mutate(ind = str_replace_all(ind, "003a", "003")) %>% 
  mutate(ind = str_replace_all(ind, "008a", "008"))


ddRAD_q <- read_csv("/local/workdir/yc2644/CV_NYC_ddRAD/structure/n380_structure.csv") %>% 
  dplyr::select(Group,Ind,P2) %>% 
  mutate(Ind = str_replace_all(Ind, "-", "_")) %>% 
  mutate(Group = str_replace_all(Group, "HH1013a", "HH13a")) %>% 
  mutate(Ind = str_replace_all(Ind, "HH1013a_0", "HH13a_")) %>% 
  dplyr::rename(p=Group,ind=Ind,all_ddRAD=P2)

table(ddRAD_q$p)
# AB812s   FI1012 FIS1013a  HH0512a  HH1013s    HH13a   HH812s   HH912s  IRV1012 NTC1013a   PM1012 
# 22       13       14       17       23       22       22       23       23       20       16 
# PRA1013a   RI812s  SV0512a  SV1013a   SV812s   SV912s  TPZ0713  WFM812s 
# 23       13       20       23       23       23       20       20 

ddRAD_snp_760_q <- read_csv("/local/workdir/yc2644/CV_NYC_ddRAD/structure/n380_SNP760_structure.csv") %>% 
  dplyr::select(Group,Ind,P2) %>% 
  mutate(Ind = str_replace_all(Ind, "-", "_")) %>% 
  mutate(Group = str_replace_all(Group, "HH1013a", "HH13a")) %>% 
  mutate(Ind = str_replace_all(Ind, "HH1013a_0", "HH13a_")) %>% 
  dplyr::rename(p=Group,ind=Ind,snp760_ddRAD=P2)

real_Q <- full_join(full_join(full_join(full_join(full_join(full_join(ddRAD_q,lcWGS_q),lcWGS_ddRAD_q),lcWGS_ddRAD_f_q),
                              lcWGS_coding_q),lcWGS_noncoding_q),ddRAD_snp_760_q)
real_Q

## only HH, SV empirical
# real_Q_subset <- real_Q %>%
#   filter(p %in% c("FI1012","HH13a","SV0512a")) %>% 
#   mutate(p =  factor(p, levels = c("FI1012","HH13a","SV0512a"))) %>%
#   arrange(p)

# real_Q_subset <- real_Q %>%
#   filter(p %in% c("HH13a","SV0512a")) %>% 
#   mutate(p =  factor(p, levels = c("HH13a","SV0512a"))) %>%
#   arrange(p)

real_Q_subset <- real_Q %>%
  filter(p %in% c("SV0512a")) %>% 
  mutate(p =  factor(p, levels = c("SV0512a"))) %>%
  arrange(p)

real_Q_subset <- real_Q_subset[complete.cases(real_Q_subset), ]
real_Q_subset
# ddRAD, Q=0.0018 
# lcWGS when only FI, Q=0.001676593
subset(real_Q_subset,all_ddRAD>0.0018)
subset(real_Q_subset,all_lcWGS>0.0017)

# correlation test --------------------------------------------------------
library("ggpubr")

# shapiro.test(real_Q_subset$all_ddRAD) #data is not normal
# shapiro.test(real_Q_subset$all_lcWGS) #data is not normal
# ggqqplot(real_Q_subset$all_ddRAD)
# # Note that, if the data are not normally distributed, 
# # it’s recommended to use the non-parametric correlation, 
# # including Spearman and Kendall rank-based correlation tests.
# 
# cor.test(real_Q_subset$all_ddRAD, real_Q_subset$all_lcWGS, method=c("kendall"),
#          exact=FALSE)
# 
# cor.test(real_Q_subset$all_ddRAD, real_Q_subset$all_lcWGS, method=c("spearman"),
#          exact=FALSE)
# 
ggscatter(real_Q_subset, x = "all_ddRAD", y = "all_lcWGS", 
          add = "reg.line", conf.int = TRUE, xlim=c(0,0.225), ylim=c(0,0.225),
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "ddRAD estimated Q", ylab = "lcWGS estimated Q") +
          #title = "Spearman Test for\nInferred Aquaculture Admixture Levels")+
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed")

dim(real_Q_subset)
is.na(real_Q_subset)

ggscatter(real_Q_subset, x = "all_ddRAD", y = "all_lcWGS_in_ddRAD", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman", xlim=c(0,0.225),ylim=c(0,0.225),
          xlab = "ddRAD estimated Q", ylab = "lcWGS in ddRAD contigs estimated Q") +
          #title = "Spearman Test for\nInferred Aquaculture Admixture Levels")+
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed")

ggscatter(real_Q_subset, x = "all_lcWGS_in_ddRAD", y = "all_lcWGS", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",xlim=c(0,0.225),ylim=c(0,0.225),
          xlab = "lcWGS in ddRAD contigs estimated Q", ylab = "lcWGS estimated Q") +
          #title = "Spearman Test for\nInferred Aquaculture Admixture Levels")+
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed")

ggscatter(real_Q_subset, x = "all_lcWGS_in_ddRAD_f", y = "all_lcWGS", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",xlim=c(0,0.225),ylim=c(0,0.225),
          xlab = "lcWGS in filtered ddRAD contigs estimated Q", ylab = "lcWGS estimated Q") +
  #title = "Spearman Test for\nInferred Aquaculture Admixture Levels")+
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed")

ggscatter(real_Q_subset, x = "all_ddRAD", y = "all_lcWGS_in_ddRAD", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",xlim=c(0,0.225),ylim=c(0,0.225),
          xlab = "ddRAD estimated Q", ylab = "lcWGS in ddRAD contigs estimated Q") +
  #title = "Spearman Test for\nInferred Aquaculture Admixture Levels")+
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed")

ggscatter(real_Q_subset, x = "all_ddRAD", y = "all_lcWGS_in_ddRAD_f", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",xlim=c(0,0.225),ylim=c(0,0.225),
          xlab = "ddRAD estimated Q", ylab = "lcWGS in filtered ddRAD contigs estimated Q") +
  #title = "Spearman Test for\nInferred Aquaculture Admixture Levels")+
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed")

ggscatter(real_Q_subset, x = "snp760_ddRAD", y = "all_lcWGS_in_ddRAD_f", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",xlim=c(0,0.225),ylim=c(0,0.225),
          xlab = "Randomly downsampled ddRAD estimated Q", ylab = "lcWGS in filtered ddRAD contigs estimated Q") +
  #title = "Spearman Test for\nInferred Aquaculture Admixture Levels")+
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed")

ggscatter(real_Q_subset, x = "snp760_ddRAD", y = "all_ddRAD", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",xlim=c(0,0.225),ylim=c(0,0.225),
          xlab = "Randomly downsampled ddRAD estimated Q", ylab = "ddRAD estimated Q") +
  #title = "Spearman Test for\nInferred Aquaculture Admixture Levels")+
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed")


ggscatter(real_Q_subset, x = "all_lcWGS_in_CDS", y = "all_lcWGS", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "lcWGS in coding regions", ylab = "lcWGS",
          title = "Spearman Test for\nInferred Aquaculture Admixture Levels") +
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed")

ggscatter(real_Q_subset, x = "all_lcWGS_in_CDS", y = "all_lcWGS_not_in_CDS", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "lcWGS in coding regions", ylab = "lcWGS not in coding regions",
          title = "Spearman Test for\nInferred Aquaculture Admixture Levels") +
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed")

wilcox.test(real_Q_subset$all_lcWGS_in_CDS, real_Q_subset$all_lcWGS_not_in_CDS,
            paired = TRUE)

# RMSE
sqrt(mean((real_Q_subset$all_ddRAD-real_Q_subset$all_lcWGS_in_ddRAD)[1:17]^2))
# 0.03952208
sqrt(mean(na.omit((real_Q_subset$all_lcWGS-real_Q_subset$all_lcWGS_in_ddRAD))^2))
# 0.02252953
sqrt(mean(na.omit((real_Q_subset$all_lcWGS-real_Q_subset$all_ddRAD))^2))
# 0.02484003


# ALL exDEBY --------------------------------------------------------------

setwd("/workdir/yc2644/CV_NYC_lcWGS/results_ngsadmix/")

target="exDEBY"
pop <- read.table(paste0("ALL_ds_Re0.2_",target,".info"), as.is=T) %>% 
  dplyr::rename(p=V1, ind=V2)

## Best K --------------------------------------------------------
K=5
all_coding <- read.table(paste0("ngsadmix_K",K,"_run1_ALL_ds_Re0.2_",
                                target,"_ALLsnplist_LDpruned_maf0.05_pval1e-6_pind0.86_cv30_noinver.codingSNP36868.qopt"))
lcWGS_coding_q <- cbind(pop, all_coding) %>% 
  dplyr::mutate(all_lcWGS_in_CDS=1-V1) %>% 
  dplyr::select(p,ind,all_lcWGS_in_CDS) %>% 
  mutate(ind = str_replace_all(ind, "003a", "003")) %>% 
  mutate(ind = str_replace_all(ind, "008a", "008"))

K=5
all_noncoding <- read.table(paste0("ngsadmix_K",K,"_run1_ALL_ds_Re0.2_",
                                   target,"_ALLsnplist_LDpruned_maf0.05_pval1e-6_pind0.86_cv30_noinver.noncodingSNP36868.qopt"))
lcWGS_noncoding_q <- cbind(pop, all_noncoding) %>% 
  dplyr::mutate(all_lcWGS_not_in_CDS=1-V4) %>% 
  dplyr::select(p,ind,all_lcWGS_not_in_CDS) %>% 
  mutate(ind = str_replace_all(ind, "003a", "003")) %>% 
  mutate(ind = str_replace_all(ind, "008a", "008"))

real_Q <- full_join(lcWGS_coding_q,lcWGS_noncoding_q)

real_Q_subset <- real_Q %>%
  filter(p %in% c("SV0512a","CT5785")) %>% 
  mutate(p =  factor(p, levels = c("SV0512a","CT5785"))) %>%
  arrange(p)
dim(real_Q_subset)

## correlation test --------------------------------------------------------
library("ggpubr")

ggscatter(real_Q_subset, x = "all_lcWGS_in_CDS", y = "all_lcWGS_not_in_CDS", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "lcWGS in coding regions (K=5)", ylab = "lcWGS not in coding regions (K=5)") +
  #title = "Spearman Test for\nInferred Aquaculture Admixture Levels")+
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed")

## K=2 --------------------------------------------------------
K=2
all_coding <- read.table(paste0("ngsadmix_K",K,"_run1_ALL_ds_Re0.2_",
                                target,"_ALLsnplist_LDpruned_maf0.05_pval1e-6_pind0.86_cv30_noinver.codingSNP36868.qopt"))
lcWGS_coding_q <- cbind(pop, all_coding) %>% 
  dplyr::mutate(all_lcWGS_in_CDS=1-V2) %>% 
  dplyr::select(p,ind,all_lcWGS_in_CDS) %>% 
  mutate(ind = str_replace_all(ind, "003a", "003")) %>% 
  mutate(ind = str_replace_all(ind, "008a", "008"))

all_noncoding <- read.table(paste0("ngsadmix_K",K,"_run1_ALL_ds_Re0.2_",
                                   target,"_ALLsnplist_LDpruned_maf0.05_pval1e-6_pind0.86_cv30_noinver.noncodingSNP36868.qopt"))
lcWGS_noncoding_q <- cbind(pop, all_noncoding) %>% 
  dplyr::mutate(all_lcWGS_not_in_CDS=1-V1) %>% 
  dplyr::select(p,ind,all_lcWGS_not_in_CDS) %>% 
  mutate(ind = str_replace_all(ind, "003a", "003")) %>% 
  mutate(ind = str_replace_all(ind, "008a", "008"))

real_Q <- full_join(lcWGS_coding_q,lcWGS_noncoding_q)

real_Q_subset <- real_Q %>%
  filter(p %in% c("SV0512a","CT5785")) %>% 
  mutate(p =  factor(p, levels = c("SV0512a","CT5785"))) %>%
  arrange(p)
dim(real_Q_subset)

## correlation test --------------------------------------------------------
library("ggpubr")

ggscatter(real_Q_subset, x = "all_lcWGS_in_CDS", y = "all_lcWGS_not_in_CDS", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "lcWGS in coding regions", ylab = "lcWGS not in coding regions") +
  #title = "Spearman Test for\nInferred Aquaculture Admixture Levels")+
  xlim(0,0.275)+
  ylim(-0.01,0.275)+
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed")
