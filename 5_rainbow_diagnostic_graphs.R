#### Workspace ####

rm(list=ls())
setwd("C:/Users/s222165681/OneDrive - Deakin University/Rainbow trout allometry/Experiment 1/Results/Processed") 


#### Packages ####

library(tidyverse)
library(glue)
library(mclust)
library(lme4)
library(car)
library(readxl)
library(lubridate)
library(hms)
library(modelr)
library(viridis)
library(glue)
library(sjPlot)
library(MuMIn)
library(patchwork)
library(ggpubr)
library(gridExtra)


#### Import Processed Data ####

tb_respirometry <- read.csv("tb_respirometry.csv")
tb_smr <- read.csv("tb_smr.csv") 
tb_MR_master <- read.csv("tb_MR_master.csv")


#### Data Preparation ####

tb_respirometry_n25 <- tb_respirometry[which(tb_respirometry$Treatment=="25°C, normoxia"),] %>%
  arrange(desc(mass))

tb_respirometry_h25 <- tb_respirometry[which(tb_respirometry$Treatment=="25°C, hyperoxia"),] %>%
  arrange(desc(mass))

tb_respirometry_h17 <- tb_respirometry[which(tb_respirometry$Treatment=="17°C, hyperoxia"),] %>%
  arrange(desc(mass))

tb_respirometry_n17 <- tb_respirometry[which(tb_respirometry$Treatment=="17°C, normoxia"),] %>%
  arrange(desc(mass))

tb_respirometry_n21 <- tb_respirometry[which(tb_respirometry$Treatment=="21°C, normoxia"),] %>%
  arrange(desc(mass))


n17_SMR <- left_join(tb_respirometry_n17, tb_smr %>% select(ID_fish, SMR, SMR_mlnd), by = "ID_fish", copy = FALSE,
                     keep = FALSE,
                     na_matches = "na")

h17_SMR <- left_join(tb_respirometry_h17, tb_smr %>% select(ID_fish, SMR, SMR_mlnd), by = "ID_fish", copy = FALSE,
                     keep = FALSE,
                     na_matches = "na")

n21_SMR <- left_join(tb_respirometry_n21, tb_smr %>% select(ID_fish, SMR, SMR_mlnd), by = "ID_fish", copy = FALSE,
                     keep = FALSE,
                     na_matches = "na")

n25_SMR <- left_join(tb_respirometry_n25, tb_smr %>% select(ID_fish, SMR, SMR_mlnd), by = "ID_fish", copy = FALSE,
                     keep = FALSE,
                     na_matches = "na")

h25_SMR <- left_join(tb_respirometry_h25, tb_smr %>% select(ID_fish, SMR, SMR_mlnd), by = "ID_fish", copy = FALSE,
                     keep = FALSE,
                     na_matches = "na")


#### Respirometry graphs ####

p_n17_SMR <- ggplot(data = n17_SMR) +
  geom_line(aes(Phase, SMR, colour = "4 SMR-lowest10%"), linetype = "dashed") +
  geom_point(aes(Phase, MR_wBR, colour = "1 Metabolic rate with background")) +
  geom_point(aes(Phase, MR, colour = "2 Metabolic Rate")) +
  geom_line(aes(Phase, SMR_mlnd, colour = "3 SMR-mlnd")) +
  scale_x_continuous(breaks = seq(0, max(tb_respirometry$Phase), 10)) +  # Use Phase variable for breaks
  facet_wrap(~factor(mass), scales = "free_y", ncol = 4) +
  labs(x = "Phase (20 min)", y = expression(Metabolic~rate~(mg~O[2]~min^-1)), title = "17°C, normoxia") +
  scale_color_manual(values = c("black", "red", "blue", "darkgreen")) +
  theme_classic() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )

p_n17_SMR <- p_n17_SMR +
  guides(
    colour = guide_legend(
      title.position = "top",  
      title.hjust = 0.5,  
      title = "Legend"  
    )
  )

print(p_n17_SMR)

p_h17_SMR <- ggplot(data = h17_SMR) +
  geom_line(aes(Phase, SMR, colour = "4 SMR-lowest10%"), linetype = "dashed") +
  geom_point(aes(Phase, MR_wBR, colour = "1 Metabolic rate with background")) +
  geom_point(aes(Phase, MR, colour = "2 Metabolic Rate")) +
  geom_line(aes(Phase, SMR_mlnd, colour = "3 SMR-mlnd")) +
  scale_x_continuous(breaks = seq(0, max(tb_respirometry$Phase), 10)) +  # Use Phase variable for breaks
  facet_wrap(~factor(mass), scales = "free_y", ncol = 4) +
  labs(x = "Phase (20 min)", y = expression(Metabolic~rate~(mg~O[2]~min^-1)), title = "17°C, hyperoxia") +
  scale_color_manual(values = c("black", "red", "blue", "darkgreen")) +
  theme_classic() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )

p_h17_SMR <- p_h17_SMR +
  guides(
    colour = guide_legend(
      title.position = "top", 
      title.hjust = 0.5,  
      title = "Legend"  
    )
  )

print(p_h17_SMR)

p_n21_SMR <- ggplot(data = n21_SMR) +
  geom_line(aes(Phase, SMR, colour = "4 SMR-lowest10%"), linetype = "dashed") +
  geom_point(aes(Phase, MR_wBR, colour = "1 Metabolic rate with background")) +
  geom_point(aes(Phase, MR, colour = "2 Metabolic Rate")) +
  geom_line(aes(Phase, SMR_mlnd, colour = "3 SMR-mlnd")) +
  scale_x_continuous(breaks = seq(0, max(tb_respirometry$Phase), 10)) +  # Use Phase variable for breaks
  facet_wrap(~factor(mass), scales = "free_y", ncol = 4) +
  labs(x = "Phase (20 min)", y = expression(Metabolic~rate~(mg~O[2]~min^-1)), title = "21°C, normoxia") +
  scale_color_manual(values = c("black", "red", "blue", "darkgreen")) +
  theme_classic() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )

p_n21_SMR <- p_n21_SMR +
  guides(
    colour = guide_legend(
      title.position = "top", 
      title.hjust = 0.5,  
      title = "Legend"  
    )
  )

print(p_n21_SMR)

p_n25_SMR <- ggplot(data = n25_SMR) +
  geom_line(aes(Phase, SMR, colour = "4 SMR-lowest10%"), linetype = "dashed") +
  geom_point(aes(Phase, MR_wBR, colour = "1 Metabolic rate with background")) +
  geom_point(aes(Phase, MR, colour = "2 Metabolic Rate")) +
  geom_line(aes(Phase, SMR_mlnd, colour = "3 SMR-mlnd")) +
  scale_x_continuous(breaks = seq(0, max(tb_respirometry$Phase), 10)) +  # Use Phase variable for breaks
  facet_wrap(~factor(mass), scales = "free_y", ncol = 4) +
  labs(x = "Phase (20 min)", y = expression(Metabolic~rate~(mg~O[2]~min^-1)), title = "25°C, normoxia") +
  scale_color_manual(values = c("black", "red", "blue", "darkgreen")) +
  theme_classic() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )

p_n25_SMR <- p_n25_SMR +
  guides(
    colour = guide_legend(
      title.position = "top", 
      title.hjust = 0.5,  
      title = "Legend"  
    )
  )

print(p_n25_SMR)

p_h25_SMR <- ggplot(data = h25_SMR) +
  geom_line(aes(Phase, SMR, colour = "4 SMR-lowest10%"), linetype = "dashed") +
  geom_point(aes(Phase, MR_wBR, colour = "1 Metabolic rate with background")) +
  geom_point(aes(Phase, MR, colour = "2 Metabolic Rate")) +
  geom_line(aes(Phase, SMR_mlnd, colour = "3 SMR-mlnd")) +
  scale_x_continuous(breaks = seq(0, max(tb_respirometry$Phase), 10)) +  # Use Phase variable for breaks
  facet_wrap(~factor(mass), scales = "free_y", ncol = 4) +
  labs(x = "Phase (20 min)", y = expression(Metabolic~rate~(mg~O[2]~min^-1)), title = "25°C, normoxia") +
  scale_color_manual(values = c("black", "red", "blue", "darkgreen")) +
  theme_classic() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )

p_h25_SMR <- p_h25_SMR +
  guides(
    colour = guide_legend(
      title.position = "top", 
      title.hjust = 0.5,  
      title = "Legend"  
    )
  )

print(p_h25_SMR)


#### Comparing SMR trials ####

n25_master <- tb_MR_master[which(tb_MR_master$Treatment=="25°C, normoxia"),] 
h25_master <- tb_MR_master[which(tb_MR_master$Treatment=="25°C, hyperoxia"),] 
h17_master <- tb_MR_master[which(tb_MR_master$Treatment=="17°C, hyperoxia"),]
n17_master <- tb_MR_master[which(tb_MR_master$Treatment=="17°C, normoxia"),] 
n21_master <- tb_MR_master[which(tb_MR_master$Treatment=="21°C, normoxia"),] 


#### SMR trial comparison ####

p_n17_SMR_comp <- ggplot(data=n17_master, aes(mass, SMR, colour = as.factor(Trial))) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  scale_x_log10(limits = c(10,500)) +
  scale_y_log10(limits = c(0.01,1)) +
  labs(x="Mass (g)", y="SMR (mgO2/min)", title = "SMR.comparison - normoxia 17" ) +
  theme_classic()

plot(p_n17_SMR_comp)

p_n21_SMR_comp <- ggplot(data=n21_master, aes(mass, SMR, colour = as.factor(Trial))) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  scale_x_log10(limits = c(10,500)) +
  scale_y_log10(limits = c(0.01,2)) +
  labs(x="Mass (g)", y="SMR (mgO2/min)", title = "SMR.comparison - normoxia 21" ) +
  theme_classic()

plot(p_n21_SMR_comp)

p_n25_SMR_comp <- ggplot(data=n25_master, aes(mass, SMR, colour = as.factor(Trial))) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  scale_x_log10(limits = c(10,500)) +
  scale_y_log10(limits = c(0.01,2)) +
  labs(x="Mass (g)", y="SMR (mgO2/min)", title = "SMR.comparison - normoxia 25" ) +
  theme_classic()

plot(p_n25_SMR_comp)

p_h17_SMR_comp <- ggplot(data=h17_master, aes(mass, SMR, colour = as.factor(Trial))) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  scale_x_log10(limits = c(10,500)) +
  scale_y_log10(limits = c(0.01,2)) +
  labs(x="Mass (g)", y="SMR (mgO2/min)", title = "SMR.comparison - hyperoxia 17" ) +
  theme_classic()

plot(p_h17_SMR_comp)

p_h25_SMR_comp <- ggplot(data=h25_master, aes(mass, SMR, colour = as.factor(Trial))) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  scale_x_log10(limits = c(10,500)) +
  scale_y_log10(limits = c(0.01,2)) +
  labs(x="Mass (g)", y="SMR (mgO2/min)", title = "SMR.comparison - hyperoxia 25" ) +
  theme_classic()

plot(p_h25_SMR_comp)
