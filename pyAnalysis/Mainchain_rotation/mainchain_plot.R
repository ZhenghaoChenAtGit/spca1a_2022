# install.packages("tidyverse")
# install.packages("rstatix")
# install.packages("ggpubr")
library(tidyverse)
library(readr)
library(ggpubr)
library(rstatix)

# Change the path here if any errors.
df_sum <- read_csv("/Mainchain_rotation/main_chain_rotation_angle.csv")

dodge <- position_dodge(width = 1)

df_sum_filtered <- df_sum %>% filter(domain != "TM11")

domain = c("A_doamin", "N_domain", "P_domain","TM1-2","TM3-4","TM5-6","TM7-10")

# Wilcox test
stat_test_bonferroni <- df_sum %>%
  group_by(domain) %>%
  pairwise_wilcox_test(angle ~ entry,) %>% 
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

# Define a palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# plot p corrected with bonferroni.


violin_bonferroni <- 
  ggviolin(df_sum_filtered, x = "entry", y = "angle", fill = "domain", 
           trim = TRUE, 
           palette = cbPalette,
           add = "boxplot", 
           add.params = list(width = 0.06,fill = "white")) +
  xlab("")+
  ylab("Angle(°)")+
  scale_x_discrete(labels=c("s2b" = "SERCA2b", 
                            "s1a" = "SERCA1a",
                            "spca1Ca" = paste("SPCA1a","Ca",sep = "\n")))+
  stat_pvalue_manual(stat_test_bonferroni,
                     label = "{p.adj.signif}",
                     hide.ns = TRUE,
                     y.position = 130,
                     step.increase = 0.15,
                     step.group.by = "domain")
facet(violin_bonferroni, facet.by = "domain", nrow = 2)

# save the plot
# ggsave("violin_bonferroni_rest.pdf", width = 16, height = 9)

