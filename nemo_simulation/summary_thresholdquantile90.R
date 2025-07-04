setwd("/local/workdir/yc2644/CV_NYC_nemo_sim/nemo_rep_sum")

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
#devtools::install_github('smin95/smplot2')
library(smplot2) #https://smin95.github.io/dataviz/raincloud-and-forest-plots.html

# empirical admixed ---------------------------------------------------------------
str_ddRAD_all <- read_csv("/local/workdir/yc2644/CV_NYC_ddRAD/structure/n380_structure.csv") %>% 
  dplyr::select(Group,P2) %>% 
  mutate(location = c(rep("Aquaculture",27), 
                      rep("Hudson River",189), 
                      rep("ddRAD\nEast River",164))) %>% dplyr::rename(p=Group,AQ=P2)
#check <- subset(str_ddRAD_all, location == "Hudson River"&p != "PRA1013a")
zero_ddRAD <- quantile(subset(str_ddRAD_all, location == "Hudson River"&p != "PRA1013a")$AQ,
                       probs=0.90)
HH_check <- subset(str_ddRAD_all, location == "Hudson River"&p != "PRA1013a")
dim(HH_check)
mean(HH_check$AQ)
dim(HH_check[HH_check$AQ>0.004,])
16/166

ggplot(HH_check, aes(x = AQ)) +
  geom_histogram(bins = 15, fill = "skyblue", color = "black") +
  geom_vline(xintercept = zero_ddRAD, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Admixture Level (Q)", y = "Count") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    panel.grid = element_blank()
  )

str_ddRAD <- str_ddRAD_all %>% 
  filter(location=="ddRAD\nEast River")
str_ddRAD %>% dplyr::group_by(location) %>% 
  dplyr::summarise(mean=mean(AQ))
# location              mean
# <chr>                <dbl>
#   1 "ddRAD\nEast River" 0.0235

str_ddRAD %>% dplyr::group_by(p) %>% 
  dplyr::summarise(max=max(AQ),min=min(AQ))
# p          max   min
# <chr>    <dbl> <dbl>
# 1 AB812s   0.144 0.001
# 2 NTC1013a 0.213 0.001
# 3 RI812s   0.243 0.001
# 4 SV0512a  0.219 0.001
# 5 SV1013a  0.154 0.001
# 6 SV812s   0.078 0    
# 7 SV912s   0.086 0    
# 8 WFM812s  0.223 0 


table(str_ddRAD$p)
# AB812s NTC1013a   RI812s  SV0512a  SV1013a   SV812s   SV912s  WFM812s 
# 22       20       13       20       23       23       23       20 
dim(subset(str_ddRAD, AQ>zero_ddRAD))[1]/dim(str_ddRAD)[1] #0.4512195
dim(subset(str_ddRAD, p=="SV0512a"&AQ>zero_ddRAD))[1]/20 #0.55
dim(subset(str_ddRAD, p=="SV1013a"&AQ>zero_ddRAD))[1]/23 #0.4782609
dim(subset(str_ddRAD, p=="SV812s"&AQ>zero_ddRAD))[1]/23 #0.3913043
dim(subset(str_ddRAD, p=="SV912s"&AQ>zero_ddRAD))[1]/23 #0.3913043
(11+11+9+9)/(20+23+23+23) #0.4494382

dim(subset(str_ddRAD, p=="AB812s"&AQ>zero_ddRAD))[1]/22
dim(subset(str_ddRAD, p=="NTC1013a"&AQ>zero_ddRAD))[1]/20
dim(subset(str_ddRAD, p=="RI812s"&AQ>zero_ddRAD))[1]/13
dim(subset(str_ddRAD, p=="WFM812s"&AQ>zero_ddRAD))[1]/20

pop <- read.table( "/workdir/yc2644/CV_NYC_lcWGS/results_ngsadmix/ALL_ds_Re0.2_exDEBY.info", as.is=T) %>% 
  dplyr::rename(p=V1, ind=V2)
q <- read.table("/workdir/yc2644/CV_NYC_lcWGS/results_ngsadmix/ngsadmix_K5_run1_ALL_ds_Re0.2_exDEBY_ALLsnplist_LDpruned_maf0.05_pval1e-6_pind0.86_cv30_noinver.qopt")
str_lcWGS_all <- cbind(pop, q) %>% 
  mutate(AQ = 1-V3) %>% filter(p=="SV0512a"|p=="CT5785"|p=="HH13a") %>%
  mutate(location = case_when(
    p == "SV0512a" ~ "lcWGS\nEast River",
    p == "CT5785" ~ "lcWGS\nLIS",
    p == "HH13a" ~ "Hudson River",
    TRUE ~ p
  )) %>% dplyr::select(p,AQ,location) %>% mutate(across(all_of("location"), as_factor))
  
zero_lcWGS <- quantile(subset(str_lcWGS_all, location == "Hudson River")$AQ,
                       probs=0.90)
HH_CHECK <- subset(str_lcWGS_all, location == "Hudson River")
range(HH_CHECK$AQ)
median(HH_CHECK$AQ)
dim(HH_CHECK)
HH_CHECK[HH_CHECK$AQ>zero_lcWGS,]
3/25
ggplot(HH_CHECK, aes(x = AQ)) +
  geom_histogram(bins = 15, fill = "skyblue", color = "black") +
  geom_vline(xintercept = zero_lcWGS, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Admixture Level (Q)", y = "Count") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    panel.grid = element_blank()
  )

str_lcWGS <- str_lcWGS_all %>% 
  filter(location=="lcWGS\nEast River"|location=="lcWGS\nLIS")
str_lcWGS %>% dplyr::group_by(location) %>% 
  dplyr::summarise(mean=mean(AQ))
# # A tibble: 2 × 2
# location              mean
# <fct>                <dbl>
# 1 "lcWGS\nEast River" 0.0785
# 2 "lcWGS\nLIS"        0.0751

str_lcWGS_ER <- str_lcWGS %>% filter(location=="lcWGS\nEast River")
str_lcWGS_CT <- str_lcWGS %>% filter(location=="lcWGS\nLIS")

range(str_lcWGS_ER$AQ)
range(str_lcWGS_CT$AQ)
sort(str_lcWGS_ER$AQ, decreasing = TRUE)
sort(str_lcWGS_CT$AQ, decreasing = TRUE)

dim(subset(str_lcWGS, AQ>zero_lcWGS))[1]/dim(str_lcWGS)[1] #0.9074074
dim(subset(str_lcWGS_ER, AQ>zero_lcWGS))[1]/dim(str_lcWGS_ER)[1] #0.9166667
dim(subset(str_lcWGS_CT, AQ>zero_lcWGS))[1]/dim(str_lcWGS_CT)[1] #0.9

str_real_nozero <- rbind(subset(str_ddRAD, AQ>zero_ddRAD),subset(str_lcWGS, AQ>zero_lcWGS))

str_real_nozero %>% dplyr::group_by(location) %>% 
  dplyr::summarise(mean=mean(AQ))
# location              mean
# <chr>                <dbl>
#   1 "ddRAD\nEast River" 0.0498
# 2 "lcWGS\nEast River" 0.0854
# 3 "lcWGS\nLIS"        0.0830

str_real_nozero %>% dplyr::group_by(location,p) %>% 
  dplyr::summarise(mean=mean(AQ))
# location            p          mean
# <chr>               <chr>     <dbl>
#   1 "ddRAD\nEast River" AB812s   0.035 
# 2 "ddRAD\nEast River" NTC1013a 0.0478
# 3 "ddRAD\nEast River" RI812s   0.0686
# 4 "ddRAD\nEast River" SV0512a  0.058 
# 5 "ddRAD\nEast River" SV1013a  0.0523
# 6 "ddRAD\nEast River" SV812s   0.0343
# 7 "ddRAD\nEast River" SV912s   0.0469
# 8 "ddRAD\nEast River" WFM812s  0.0586
# 9 "lcWGS\nEast River" SV0512a  0.0854
# 10 "lcWGS\nLIS"        CT5785   0.0830

dim(str_real_nozero %>% filter(location=="ddRAD\nEast River")) #48
74/164
p1 <- ggplot(str_real_nozero,
       mapping = aes(x = location, y = AQ,fill=location)) +
  sm_raincloud(sep_level = 2,
               boxplot.params = list(outlier.shape = NA,color=NA,fill=NA),
               #err.params = list(size = 1.2),
               violin.params = list(scale="width",color="black",alpha=0.2), 
               #errorbar_type = "sd",
               point_jitter_width = 0.3,
               point.params = list(shape = 1, color = "lightslateblue",
                                   size = 1.5, alpha=0.45))+
  scale_fill_manual(values=rep("lightslateblue",8))+
  geom_point(stat = "summary", fun = "median", size = 2, color = "red",na.rm = F) +
  ylab(label = "Q values for admixed individuals") +
  xlab(label = NULL) +
  #labs(title = "Empirical Data") +
  theme_minimal() +
  theme(axis.text.x = element_text(size=12,hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size=14,hjust = 0.5, vjust = 0.5),
        axis.title = element_text(size = 14),
        strip.text = element_text(size=14),
        plot.title = element_text(size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  coord_cartesian(ylim = c(0, 1))
p1

# cons high admixed box plots -------------------------------------------------------

# Create an empty list to store data frames
all_tot_Q <- list()

directory_path <- "/local/workdir/yc2644/CV_NYC_nemo_sim/nemo_rep_sum/"
file_pattern <- "^2pop_ntrl_eq10400_s3270_cons_(g1|g6|g11|g21|g31)_total_rep10_tidy_Q\\.tsv$"
matching_files <- list.files(path = directory_path, pattern = file_pattern, full.names = TRUE)
print(matching_files) #1scen*4gen

# Function to read and combine TSV files
read_and_combine_tsv <- function(file_path) {
  read_tsv(file_path)
}

# Use purrr::map_dfr to read and combine the TSV files into one data frame
combined_cons_tot_Q <- map_dfr(matching_files, read_and_combine_tsv) %>% 
  mutate(gen_to_admix=gen_to_admix-1) %>% 
  mutate(gen_to_admix = as.factor(gen_to_admix))
# Print or work with the combined data
print(combined_cons_tot_Q) #21600=1scen*4gen*10rep*6rate*90ind

immigration_rate_list <- c(`cons0.001`="m=0.0001",`cons0.002`="m=0.0002",
                           `cons0.005`="m=0.0005",`cons0.01`="m=0.001",
                           `cons0.02`="m=0.002",`cons0.05`="m=0.005")

zero <- quantile(subset(combined_cons_tot_Q,pop == "W2")$AQ_assignment_rate,
                 probs=0.90)

fst=0.360
ggplot(subset(combined_cons_tot_Q,pop=="W1"&gen_to_admix!=0&AQ_assignment_rate>zero&admix_scenario%in%c("cons0.002","cons0.02","cons0.05")),
       mapping = aes(x = gen_to_admix, 
                     y = AQ_assignment_rate,
                     fill=as.factor(gen_to_admix))) +
  sm_raincloud(sep_level = 2,
               boxplot.params = list(outlier.shape = NA,color=NA,fill=NA),
               #err.params = list(size = 1.2),
               violin.params = list(scale="width",color="black",alpha=0.2), 
               #errorbar_type = "sd",
               point_jitter_width = 0.3,
               point.params = list(shape = 1, color = "lightslateblue",
                                   size = 1.5, alpha=0.45))+
  scale_fill_manual(values=rep("lightslateblue",8))+
  geom_point(stat = "summary", fun = "median", size = 2, color = "red",na.rm = F) +
  geom_line(stat = "summary", fun = "median", aes(group = 1), 
            color = "red", size = 0.8, na.rm = F,linetype = "dashed")+
  ylim(c(0,1))+
  #geom_hline(yintercept=0.08, color="red") +
  # annotate("segment", x = 0, xend = 0.25, y = 0.08, yend = 0.08, 
  #          arrow = arrow(length = unit(0.2, "cm")), color = "red")+
  # geom_errorbar(stat = "summary",fun.data = "mean_cl_boot",width = 0.2,
  #   size = 1,color = "red") +
  xlab(label = "Generation(s) since the first migration event") +
  ylab(label = "Q values for admixed individuals") +
  labs(title = paste0("Initial Mean Fst=",fst)) +
  theme_minimal() +
  theme(axis.text = element_text(size=14,hjust = 0.5, vjust = 0.5),
        axis.title = element_text(size = 14),
        strip.text = element_text(size=14),
        plot.title = element_text(size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  facet_wrap(~admix_scenario, nrow=1,
             #scales = "free_y",
             labeller = as_labeller(immigration_rate_list))


# cons high non-admixed ---------------------------------------------------

combined_cons_non <- subset(combined_cons_tot_Q,pop=="W1"&admix_scenario%in%c("cons0.002","cons0.02","cons0.05")) %>% 
  group_by(gen_to_admix,admix_scenario,simulation_replicate) %>% 
  dplyr::summarise(num_zeros = sum(AQ_assignment_rate<=zero)) %>% 
  ungroup() %>% 
  mutate(freq_unadmixed = num_zeros/30) %>% 
  dplyr::select(gen_to_admix,admix_scenario,simulation_replicate,freq_unadmixed)
  
p11 <- ggplot(subset(combined_cons_non, gen_to_admix!=0),
       mapping=aes(x=gen_to_admix, y=1-freq_unadmixed),color="lightslateblue")+
  geom_point(stat = "summary", fun = "mean", size = 2, na.rm = F) +
  geom_errorbar(stat = "summary",fun.data = "mean_cl_boot",width = 0.1,size = 1)+
  xlab(label = "Generation(s) since the first migration event") +
  labs(y = "Proportions of admixed individuals\nin native population") +
  labs(title = paste0("Initial Mean Fst=",fst)) +
  theme_minimal() +
  theme(axis.text = element_text(size=14, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(size = 16),
        strip.text = element_text(size=16),
        plot.title = element_text(size = 16))+
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.background = element_rect(fill = "white", color = "black", size = 0.5))+
  coord_cartesian(ylim = c(0, 1))+
  facet_wrap(~admix_scenario, nrow=1,
             #scales = "free_y",
             labeller = as_labeller(immigration_rate_list))


# cons low admixed box plots -------------------------------------------------------

# Create an empty list to store data frames
all_tot_Q <- list()

directory_path <- "/local/workdir/yc2644/CV_NYC_nemo_sim/nemo_rep_sum/"
file_pattern <- "^2pop_ntrl_eq10400_s280_cons_(g1|g6|g11|g21|g31)_total_rep10_tidy_Q\\.tsv$"
matching_files <- list.files(path = directory_path, pattern = file_pattern, full.names = TRUE)
print(matching_files) #1scen*4gen

# Function to read and combine TSV files
read_and_combine_tsv <- function(file_path) {
  read_tsv(file_path)
}

# Use purrr::map_dfr to read and combine the TSV files into one data frame
combined_cons_tot_Q <- map_dfr(matching_files, read_and_combine_tsv) %>% 
  mutate(gen_to_admix=gen_to_admix-1) %>% 
  mutate(gen_to_admix = as.factor(gen_to_admix))
# Print or work with the combined data
print(combined_cons_tot_Q) #21600=1scen*4gen*10rep*6rate*90ind

immigration_rate_list <- c(`cons0.001`="m=0.0001",`cons0.002`="m=0.0002",
                           `cons0.005`="m=0.0005",`cons0.01`="m=0.001",
                           `cons0.02`="m=0.002",`cons0.05`="m=0.005")

zero <- quantile(subset(combined_cons_tot_Q,pop == "W2")$AQ_assignment_rate, 
                 probs=0.90)
subset(combined_cons_tot_Q,pop == "W2") %>% 
  dplyr::group_by(simulation_replicate) %>% 
  dplyr::summarise(max=max(AQ_assignment_rate),
                   mean=mean(AQ_assignment_rate))

fst=0.070
p2 <- ggplot(subset(combined_cons_tot_Q,
                    pop=="W1"&gen_to_admix!=0 & AQ_assignment_rate>zero & admix_scenario%in%c("cons0.002","cons0.02","cons0.05")),
       mapping = aes(x = gen_to_admix, 
                     y = AQ_assignment_rate,
                     fill=as.factor(gen_to_admix))) +
    sm_raincloud(sep_level = 2,
               boxplot.params = list(outlier.shape = NA,color=NA,fill=NA),
               #err.params = list(size = 1.2),
               violin.params = list(scale="width",color="black",alpha=0.2), 
               #errorbar_type = "sd",
               point_jitter_width = 0.3,
               point.params = list(shape = 1, color = "lightslateblue",
                                   size = 1.5, alpha=0.45))+
  scale_fill_manual(values=rep("lightslateblue",8))+
  geom_point(stat = "summary", fun = "median", size = 2, color = "red",na.rm = F) +
  geom_line(stat = "summary", fun = "median", aes(group = 1), 
            color = "red", size = 0.8, na.rm = F,linetype = "dashed")+
  ylim(c(0,1))+
  xlab(label = "Generation(s) since the first migration event") +
  ylab(label = NULL) +
  labs(title = paste0("Initial Mean Fst=",fst)) +
  theme_minimal() +
  theme(axis.text = element_text(size=14,hjust = 0.5, vjust = 0.5),
        axis.title = element_text(size = 14),
        strip.text = element_text(size=14),
        plot.title = element_text(size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  facet_wrap(~admix_scenario, nrow=1,
             #scales = "free_y",
             labeller = as_labeller(immigration_rate_list))
p2

# cons low non-admixed ---------------------------------------------------

combined_cons_non <- subset(combined_cons_tot_Q,pop=="W1"&admix_scenario%in%c("cons0.002","cons0.02","cons0.05")) %>% 
  group_by(gen_to_admix,admix_scenario,simulation_replicate) %>% 
  dplyr::summarise(num_zeros = sum(AQ_assignment_rate<=zero)) %>% 
  ungroup() %>% 
  mutate(freq_unadmixed = num_zeros/30) %>% 
  dplyr::select(gen_to_admix,admix_scenario,simulation_replicate,freq_unadmixed)

p12 <- ggplot(subset(combined_cons_non,gen_to_admix!=0), 
       mapping=aes(x=gen_to_admix, y=1-freq_unadmixed),color="lightslateblue")+
  geom_point(stat = "summary", fun = "mean", size = 2, na.rm = F) +
  geom_errorbar(stat = "summary",fun.data = "mean_cl_boot",width = 0.1,size = 1)+
  xlab(label = "Generation(s) since the first migration event") +
  labs(y = "Proportions of admixed individuals\nin native population") +
  labs(title = paste0("Initial Mean Fst=",fst)) +
  theme_minimal() +
  theme(axis.text = element_text(size=14,hjust = 0.5, vjust = 0.5),
        axis.title = element_text(size = 16),
        strip.text = element_text(size=16),
        plot.title = element_text(size = 16))+
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.background = element_rect(fill = "white", color = "black", size = 0.5))+
  coord_cartesian(ylim = c(0, 1))+
  facet_wrap(~admix_scenario, nrow=1,
             #scales = "free_y",
             labeller = as_labeller(immigration_rate_list))


# single high - admixed box plots -------------------------------------------------------

# Create an empty list to store data frames
all_tot_Q <- list()

directory_path <- "/local/workdir/yc2644/CV_NYC_nemo_sim/nemo_rep_sum/"
#file_pattern <- "^2pop_ntrl_eq10400_s1000_single_(g1|g2|g3|g5|g10|g15|g20|g30)_total_rep10_tidy_Q\\.tsv$"
file_pattern <- "^2pop_ntrl_eq10400_s3270_single_(g1|g2|g3|g4|g6|g11|g16|g21)_total_rep10_tidy_Q\\.tsv$"
#file_pattern <- "^2pop_ntrl_eq10400_s280_single_(g1|g2|g3|g4|g6|g11|g16|g21)_total_rep10_tidy_Q\\.tsv$"
fst=0.360 #0.070 or 0.360

matching_files <- list.files(path = directory_path, pattern = file_pattern, full.names = TRUE)
print(matching_files) #1scen*8gen

# Function to read and combine TSV files
read_and_combine_tsv <- function(file_path) {
  read_tsv(file_path)
}

# Use purrr::map_dfr to read and combine the TSV files into one data frame
combined_cons_tot_Q <- map_dfr(matching_files, read_and_combine_tsv) %>% 
  mutate(gen_to_admix=gen_to_admix-1) %>% 
  mutate(gen_to_admix = as.factor(gen_to_admix))
# Print or work with the combined data
print(combined_cons_tot_Q) #21600=1scen*4gen*10rep*6rate*90ind

immigration_rate_list <- c(`single0.1`="m=0.01",`single0.2`="m=0.02",`single0.5`="m=0.05",
                           `single0.8`="m=0.08",`single1`="m=0.1")

zero <- quantile(subset(combined_cons_tot_Q,pop == "W2")$AQ_assignment_rate,probs=0.90)

p3 <- ggplot(subset(combined_cons_tot_Q,pop=="W1"&AQ_assignment_rate>zero&admix_scenario%in%c("single0.1","single0.5","single1")),
       mapping = aes(x = gen_to_admix, 
                     y = AQ_assignment_rate,
                     fill=as.factor(gen_to_admix))) +
  sm_raincloud(sep_level = 2,
               boxplot.params = list(outlier.shape = NA,color=NA,fill=NA),
               #err.params = list(size = 1.2),
               violin.params = list(scale="width",color="black",alpha=0.2), 
               #errorbar_type = "sd",
               point_jitter_width = 0.3,
               point.params = list(shape = 1, color = "lightslateblue",
                                   size = 1.5, alpha=0.45))+
  scale_fill_manual(values=rep("lightslateblue",8))+
  geom_point(data=subset(combined_cons_tot_Q,
                         pop=="W1"&gen_to_admix!=0&AQ_assignment_rate>zero&admix_scenario%in%c("single0.1","single0.5","single1")),
             stat = "summary", fun = "median", size = 2, color = "red",na.rm = F) +
  geom_line(data=subset(combined_cons_tot_Q,
                        pop=="W1"&gen_to_admix!=0&AQ_assignment_rate>zero&admix_scenario%in%c("single0.1","single0.5","single1")),
            stat = "summary", fun = "median", aes(group = 1), 
            color = "red", size = 0.8, na.rm = F,linetype = "dashed")+
  # geom_errorbar(stat = "summary",fun.data = "mean_cl_boot",width = 0.2,
  #   size = 1,color = "red") +
  xlab(label = "Generation(s) since the first migration event") +
  ylab(label = "Q values for admixed individuals") +
  labs(title = paste0("Initial Mean Fst=",fst)) +
  theme_minimal() +
  theme(axis.text = element_text(size=14,hjust = 0.5, vjust = 0.5),
        axis.title = element_text(size = 14),
        strip.text = element_text(size=14),
        plot.title = element_text(size = 14),
        legend.position = "none",
        #axis.text.x = element_text(size=12,hjust = 0.5, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  facet_wrap(~admix_scenario, nrow=1,
             #scales = "free_y",
             labeller = as_labeller(immigration_rate_list)
  )
p3

# single high non-admixed ---------------------------------------------------

combined_cons_non <- subset(combined_cons_tot_Q,pop=="W1"&admix_scenario%in%c("single0.1","single0.5","single1")) %>% 
  group_by(gen_to_admix,admix_scenario,simulation_replicate) %>% 
  dplyr::summarise(num_zeros = sum(AQ_assignment_rate<=zero)) %>% 
  ungroup() %>% 
  mutate(freq_unadmixed = num_zeros/30) %>% 
  dplyr::select(gen_to_admix,admix_scenario,simulation_replicate,freq_unadmixed)

p13 <- ggplot(combined_cons_non, 
       mapping=aes(x=gen_to_admix, y=1-freq_unadmixed),color="lightslateblue")+
  geom_point(stat = "summary", fun = "mean", size = 2, na.rm = F) +
  geom_errorbar(stat = "summary",fun.data = "mean_cl_boot",width = 0.1,size = 1)+
  xlab(label = "Generation(s) since the first migration event") +
  labs(y = "Proportions of admixed individuals\nin native population") +
  labs(title = paste0("Initial Mean Fst=",fst)) +
  theme_minimal() +
  theme(axis.text = element_text(size=14,hjust = 0.5, vjust = 0.5),
        axis.title = element_text(size = 16),
        strip.text = element_text(size=16),
        plot.title = element_text(size = 16))+
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.background = element_rect(fill = "white", color = "black", size = 0.5))+
  coord_cartesian(ylim = c(0, 1))+
  facet_wrap(~admix_scenario, nrow=1,
             #scales = "free_y",
             labeller = as_labeller(immigration_rate_list))



# single low - admixed box plots -------------------------------------------------------

file_pattern <- "^2pop_ntrl_eq10400_s280_single_(g1|g2|g3|g4|g6|g11|g16|g21)_total_rep10_tidy_Q\\.tsv$"
fst=0.070 #0.070 or 0.360

matching_files <- list.files(path = directory_path, pattern = file_pattern, full.names = TRUE)
print(matching_files) #1scen*8gen

# Function to read and combine TSV files
read_and_combine_tsv <- function(file_path) {
  read_tsv(file_path)
}

# Use purrr::map_dfr to read and combine the TSV files into one data frame
combined_cons_tot_Q <- map_dfr(matching_files, read_and_combine_tsv) %>% 
  mutate(gen_to_admix=gen_to_admix-1) %>% 
  mutate(gen_to_admix = as.factor(gen_to_admix))
# Print or work with the combined data
print(combined_cons_tot_Q) #21600=1scen*4gen*10rep*6rate*90ind

immigration_rate_list <- c(`single0.1`="m=0.01",`single0.2`="m=0.02",`single0.5`="m=0.05",
                           `single0.8`="m=0.08",`single1`="m=0.1")

zero <- quantile(subset(combined_cons_tot_Q,pop == "W2")$AQ_assignment_rate,probs=0.90)

p4 <- ggplot(subset(combined_cons_tot_Q,pop=="W1"&AQ_assignment_rate>zero&admix_scenario%in%c("single0.1","single0.5","single1")),
             mapping = aes(x = gen_to_admix, 
                           y = AQ_assignment_rate,
                           fill=as.factor(gen_to_admix))) +
  sm_raincloud(sep_level = 2,
               boxplot.params = list(outlier.shape = NA,color=NA,fill=NA),
               #err.params = list(size = 1.2),
               violin.params = list(scale="width",color="black",alpha=0.2), 
               #errorbar_type = "sd",
               point_jitter_width = 0.3,
               point.params = list(shape = 1, color = "lightslateblue",
                                   size = 1.5, alpha=0.45))+
  ylim(c(0,1))+
  scale_fill_manual(values=rep("lightslateblue",8))+
  geom_point(data=subset(combined_cons_tot_Q,
                         pop=="W1"&gen_to_admix!=0&AQ_assignment_rate>zero&admix_scenario%in%c("single0.1","single0.5","single1")),
             stat = "summary", fun = "median", size = 2, color = "red",na.rm = F) +
  geom_line(data=subset(combined_cons_tot_Q,
                        pop=="W1"&gen_to_admix!=0&AQ_assignment_rate>zero&admix_scenario%in%c("single0.1","single0.5","single1")),
            stat = "summary", fun = "median", aes(group = 1), 
            color = "red", size = 0.8, na.rm = F,linetype = "dashed")+
  # geom_errorbar(stat = "summary",fun.data = "mean_cl_boot",width = 0.2,
  #   size = 1,color = "red") +
  xlab(label = "Generation(s) since the first migration event") +
  ylab(label = "Q values for admixed individuals") +
  labs(title = paste0("Initial Mean Fst=",fst,")")) +
  theme_minimal() +
  theme(axis.text = element_text(size=14,hjust = 0.5, vjust = 0.5),
        axis.title = element_text(size = 14),
        strip.text = element_text(size=14),
        plot.title = element_text(size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  facet_wrap(~admix_scenario, nrow=1,
             #scales = "free_y",
             labeller = as_labeller(immigration_rate_list)
  )
p4

# single low non-admixed ---------------------------------------------------

combined_cons_non <- subset(combined_cons_tot_Q,pop=="W1"&admix_scenario%in%c("single0.1","single0.5","single1")) %>% 
  group_by(gen_to_admix,admix_scenario,simulation_replicate) %>% 
  dplyr::summarise(num_zeros = sum(AQ_assignment_rate<=zero)) %>% 
  ungroup() %>% 
  mutate(freq_unadmixed = num_zeros/30) %>% 
  dplyr::select(gen_to_admix,admix_scenario,simulation_replicate,freq_unadmixed)

p14 <- ggplot(combined_cons_non, 
       mapping=aes(x=gen_to_admix, y=1-freq_unadmixed),color="lightslateblue")+
  geom_point(stat = "summary", fun = "mean", size = 2, na.rm = F) +
  geom_errorbar(stat = "summary",fun.data = "mean_cl_boot",width = 0.1,size = 1)+
  xlab(label = "Generation(s) since the first migration event") +
  labs(y = "Proportions of admixed individuals\nin native population") +
  labs(title = paste0("Initial Mean Fst=",fst)) +
  theme_minimal() +
  theme(axis.text = element_text(size=14,hjust = 0.5, vjust = 0.5),
        axis.title = element_text(size = 16),
        strip.text = element_text(size=16),
        plot.title = element_text(size = 16))+
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.background = element_rect(fill = "white", color = "black", size = 0.5))+
  coord_cartesian(ylim = c(0, 1))+
  facet_wrap(~admix_scenario, nrow=1,
             #scales = "free_y",
             labeller = as_labeller(immigration_rate_list))

library(patchwork)

layout <-" 
AABBBB
CCCCCC
DDDDDD
"

#p1+p2+p4+p3+ plot_layout(design=layout, axis_titles = "collect")#,axes = "collect")

layout1 <-" 
AABBBB
"
p5 <- p1+p2+ plot_layout(design=layout1,axes = "collect")
layout2 <-" 
AAAAAA
BBBBBB
"
p6 <- p4+p3+ plot_layout(design=layout2,axes = "collect")
layout3 <-" 
AAAAAA
BBBBBB
BBBBBB
"
p5/p6+ plot_layout(design=layout3)

layout4 <-" 
AABB
CCDD
"
p14+p13+p12+p11+ plot_layout(design=layout4,axes = "collect")
