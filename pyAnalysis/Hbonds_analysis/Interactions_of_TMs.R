library(tidyverse)
library(readr)

# Change the path to Interactions_of_TMs_sum.csv
df_sum_twoState <- read_csv("pyAnalysis/Hbonds_Analysis/output/Interactions_of_TMs_sum.csv")

df_inter_both <- df_sum_twoState %>% filter(interdomain == "TRUE")

dodge <- position_dodge(width = 1)

list_tm = c("TM1",
            "TM2",
            "TM3",
            "TM4",
            "TM5",
            "TM6",
            "TM7",
            "TM8",
            "TM9",
            "TM10")
color_tm = c("#670A1F",
             "#053061",
             "#B2182B",
             "#D6604D",
             "#F4A582",
             "#FDDBC7",
             "#D1E5F0",
             "#92C5DE",
             "#4393C3",
             "#2166AC"
             )


# plot two states

df_inter_both %>%
  filter(domain %in% list_tm) %>% 
  ggplot(aes(x = entry, y = total_count, fill = state)) +
  geom_bar(stat = "identity",
           position = dodge,
           ) +
  geom_text(aes(y = total_count,
                label=round(total_count, 2)), 
            position = dodge,
            vjust = -0.5)+
  xlab("") +
  ylab("Interdomain hydrogen bonds") +
  theme_classic() +
  facet_wrap(. ~ factor(domain, levels = list_tm),nrow = 2) +
  scale_fill_manual(values = c("#5698CF","#FFB190"))+
  theme(axis.text.x = element_text(angle = 45, vjust = -0.01))



