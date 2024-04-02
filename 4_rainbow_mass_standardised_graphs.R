#### Workspace ####

rm(list=ls())
setwd("C:/Users/s222165681/OneDrive - Deakin University/Rainbow trout allometry/Experiment 1/Results/Processed") 


#### Packages ####

library(tidyverse) #loads ggplot2, dplyr, tidyr, tibble, stringr
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

tb_smr <- read.csv("tb_smr.csv")
tb_mmr <- read.csv("tb_mmr.csv")
tb_MR_master <- read.csv("tb_MR_master.csv") 
 

#### Mass-standardised SMR calculation ####

lm_log_smr_17_norm <- lm(log10(SMR) ~ log10(mass), data = filter(tb_smr, Treatment == "17°C, normoxia"))
lm_log_smr_25_norm <- lm(log10(SMR) ~ log10(mass), data = filter(tb_smr, Treatment == "25°C, normoxia"))
lm_log_smr_17_hyper <- lm(log10(SMR) ~ log10(mass), data = filter(tb_smr, Treatment == "17°C, hyperoxia"))
lm_log_smr_25_hyper <- lm(log10(SMR) ~ log10(mass), data = filter(tb_smr, Treatment == "25°C, hyperoxia"))
lm_log_smr_21_norm <- lm(log10(SMR) ~ log10(mass), data = filter(tb_smr, Treatment == "21°C, normoxia"))

mass = data.frame(mass = mean(tb_smr$mass)) 

SMR_mean_predicted_17_norm = predict(lm_log_smr_17_norm, newdata = mass) 
SMR_mean_predicted_25_norm = predict(lm_log_smr_25_norm, newdata = mass)
SMR_mean_predicted_17_hyper = predict(lm_log_smr_17_hyper, newdata = mass)
SMR_mean_predicted_25_hyper = predict(lm_log_smr_25_hyper, newdata = mass)
SMR_mean_predicted_21_norm = predict(lm_log_smr_21_norm, newdata = mass)

tb_smr_17_norm <- tb_smr[which(tb_smr$Treatment == "17°C, normoxia"),] %>% 
  add_residuals(lm_log_smr_17_norm) 

tb_smr_17_norm <- 
  tb_smr_17_norm %>%
  mutate(
    SMR_mass_17_norm = 10^(SMR_mean_predicted_17_norm + resid),
    resid = NULL
  )

tb_smr_25_norm <- tb_smr[which(tb_smr$Treatment == "25°C, normoxia"),] %>% 
  add_residuals(lm_log_smr_25_norm) 

tb_smr_25_norm <- 
  tb_smr_25_norm %>%
  mutate(
    SMR_mass_25_norm = 10^(SMR_mean_predicted_25_norm + resid),
    resid = NULL
  )

tb_smr_17_hyper <- tb_smr[which(tb_smr$Treatment == "17°C, hyperoxia"),] %>% 
  add_residuals(lm_log_smr_17_hyper) 

tb_smr_17_hyper <- 
  tb_smr_17_hyper %>%
  mutate(
    SMR_mass_17_hyper = 10^(SMR_mean_predicted_17_hyper + resid),
    resid = NULL
  )


tb_smr_25_hyper <- tb_smr[which(tb_smr$Treatment == "25°C, hyperoxia"),] %>% 
  add_residuals(lm_log_smr_25_hyper) 

tb_smr_25_hyper <- 
  tb_smr_25_hyper %>%
  mutate(
    SMR_mass_25_hyper = 10^(SMR_mean_predicted_25_hyper + resid),
    resid = NULL
  )

tb_smr_21_norm <- tb_smr[which(tb_smr$Treatment == "21°C, normoxia"),] %>% 
  add_residuals(lm_log_smr_21_norm) 

tb_smr_21_norm <- 
  tb_smr_21_norm %>%
  mutate(
    SMR_mass_21_norm = 10^(SMR_mean_predicted_21_norm + resid),
    resid = NULL
  )


#### Mass-standardised SMR graph ####

p_SMR_mass <- ggplot() +
  geom_boxplot(aes(x = Treatment, y = SMR_mass_17_norm, fill = "17°C, normoxia"), data = tb_smr_17_norm, outlier.shape = NA, width = .15, alpha = 0.7, position = position_nudge(x = .25)) +
  gghalves::geom_half_point(aes(x = Treatment, y = SMR_mass_17_norm, color = "17°C, normoxia", fill = "17°C, normoxia"), data = tb_smr_17_norm, side = "l", range_scale = .4, alpha = .3, position = position_nudge(x = .25)) +
  
  geom_boxplot(aes(x = Treatment, y = SMR_mass_25_norm, fill = "25°C, normoxia"), data = tb_smr_25_norm, outlier.shape = NA, width = .15, alpha = 0.5, position = position_nudge(x = .25)) +
  gghalves::geom_half_point(aes(x = Treatment, y = SMR_mass_25_norm, color = "25°C, normoxia", fill = "25°C, normoxia"), data = tb_smr_25_norm, side = "l", range_scale = .4, alpha = .3, position = position_nudge(x = .25)) +
  
  geom_boxplot(aes(x = Treatment, y = SMR_mass_17_hyper, fill = "17°C, hyperoxia"), data = tb_smr_17_hyper, outlier.shape = NA, width = .15, alpha = 0.5, position = position_nudge(x = -.25)) +
  gghalves::geom_half_point(aes(x = Treatment, y = SMR_mass_17_hyper, color = "17°C, hyperoxia", fill = "17°C, hyperoxia"), data = tb_smr_17_hyper, side = "l", range_scale = .4, alpha = .3, position = position_nudge(x = -.25)) +
  
  geom_boxplot(aes(x = Treatment, y = SMR_mass_25_hyper, fill = "25°C, hyperoxia"), data = tb_smr_25_hyper, outlier.shape = NA, width = .15, alpha = 0.5, position = position_nudge(x = -.25)) +
  gghalves::geom_half_point(aes(x = Treatment, y = SMR_mass_25_hyper, color = "25°C, hyperoxia", fill = "25°C, hyperoxia"), data = tb_smr_25_hyper, side = "l", range_scale = .4, alpha = .3, position = position_nudge(x = -.25)) +
  
  geom_boxplot(aes(x = Treatment, y = SMR_mass_21_norm, fill = "21°C, normoxia"), data = tb_smr_21_norm, outlier.shape = NA, width = .15, alpha = 0.5) +
  gghalves::geom_half_point(aes(x = Treatment, y = SMR_mass_21_norm, color = "21°C, normoxia", fill = "21°C, normoxia"), data = tb_smr_21_norm, side = "l", range_scale = .4, alpha = .3) +
  
  labs(x = "Treatment", y = expression(SMR~(mg~O[2]~min^-1~88.5~g^-1))) +
  theme_classic() +
  scale_x_discrete(limits = c("17°C, normoxia", "17°C, hyperoxia", "21°C, normoxia", "25°C, normoxia", "25°C, hyperoxia")) +
  scale_fill_manual(name = "Treatment", values = c(
    "17°C, normoxia" = "slateblue1",
    "25°C, normoxia" = "darkorange2",
    "17°C, hyperoxia" = "blue",
    "25°C, hyperoxia" = "red3",
    "21°C, normoxia" = "orange"),
    breaks = c("17°C, normoxia", "17°C, hyperoxia", "21°C, normoxia", "25°C, normoxia", "25°C, hyperoxia")) +
  scale_color_manual(name = "Treatment", values = c(
    "17°C, normoxia" = "slateblue1",
    "25°C, normoxia" = "darkorange2",
    "17°C, hyperoxia" = "blue",
    "25°C, hyperoxia" = "red3",
    "21°C, normoxia" = "orange"),
    breaks = c("17°C, normoxia", "17°C, hyperoxia", "21°C, normoxia", "25°C, normoxia", "25°C, hyperoxia"))


#### Mass-standardised MMR calculation ####

lm_log_mmr_17_norm <- lm(log10(MMR) ~ log10(mass), data = filter(tb_mmr, Treatment == "17°C, normoxia"))
lm_log_mmr_25_norm <- lm(log10(MMR) ~ log10(mass), data = filter(tb_mmr, Treatment == "25°C, normoxia"))
lm_log_mmr_17_hyper <- lm(log10(MMR) ~ log10(mass), data = filter(tb_mmr, Treatment == "17°C, hyperoxia"))
lm_log_mmr_25_hyper <- lm(log10(MMR) ~ log10(mass), data = filter(tb_mmr, Treatment == "25°C, hyperoxia"))
lm_log_mmr_21_norm <- lm(log10(MMR) ~ log10(mass), data = filter(tb_mmr, Treatment == "21°C, normoxia"))

mass = data.frame(mass = mean(tb_mmr$mass)) 

MMR_mean_predicted_17_norm = predict(lm_log_mmr_17_norm, newdata = mass) 
MMR_mean_predicted_25_norm = predict(lm_log_mmr_25_norm, newdata = mass)
MMR_mean_predicted_17_hyper = predict(lm_log_mmr_17_hyper, newdata = mass)
MMR_mean_predicted_25_hyper = predict(lm_log_mmr_25_hyper, newdata = mass)
MMR_mean_predicted_21_norm = predict(lm_log_mmr_21_norm, newdata = mass)

tb_mmr_17_norm <- tb_mmr[which(tb_mmr$Treatment == "17°C, normoxia"),] %>% 
  add_residuals(lm_log_mmr_17_norm) 

tb_mmr_17_norm <- 
  tb_mmr_17_norm %>%
  mutate(
    MMR_mass_17_norm = 10^(MMR_mean_predicted_17_norm + resid),
    resid = NULL
  )

tb_mmr_25_norm <- tb_mmr[which(tb_mmr$Treatment == "25°C, normoxia"),] %>% 
  add_residuals(lm_log_mmr_25_norm) 

tb_mmr_25_norm <- 
  tb_mmr_25_norm %>%
  mutate(
    MMR_mass_25_norm = 10^(MMR_mean_predicted_25_norm + resid),
    resid = NULL
  )
tb_mmr_17_hyper <- tb_mmr[which(tb_mmr$Treatment == "17°C, hyperoxia"),] %>% 
  add_residuals(lm_log_mmr_17_hyper) 

tb_mmr_17_hyper <- 
  tb_mmr_17_hyper %>%
  mutate(
    MMR_mass_17_hyper = 10^(MMR_mean_predicted_17_hyper + resid),
    resid = NULL
  )


tb_mmr_25_hyper <- tb_mmr[which(tb_mmr$Treatment == "25°C, hyperoxia"),] %>% 
  add_residuals(lm_log_mmr_25_hyper) 

tb_mmr_25_hyper <- 
  tb_mmr_25_hyper %>%
  mutate(
    MMR_mass_25_hyper = 10^(MMR_mean_predicted_25_hyper + resid),
    resid = NULL
  )

tb_mmr_21_norm <- tb_mmr[which(tb_mmr$Treatment == "21°C, normoxia"),] %>% 
  add_residuals(lm_log_mmr_21_norm) 

tb_mmr_21_norm <- 
  tb_mmr_21_norm %>%
  mutate(
    MMR_mass_21_norm = 10^(MMR_mean_predicted_21_norm + resid),
    resid = NULL
  )


#### Mass-standardised MMR graph ####

p_MMR_mass <- ggplot() +
  geom_boxplot(aes(x = Treatment, y = MMR_mass_17_norm, fill = "17°C, normoxia"), data = tb_mmr_17_norm, outlier.shape = NA, width = .15, alpha = 0.7, position = position_nudge(x = .25)) +
  gghalves::geom_half_point(aes(x = Treatment, y = MMR_mass_17_norm, color = "17°C, normoxia", fill = "17°C, normoxia"), data = tb_mmr_17_norm, side = "l", range_scale = .4, alpha = .3, position = position_nudge(x = .25)) +
  
  geom_boxplot(aes(x = Treatment, y = MMR_mass_25_norm, fill = "25°C, normoxia"), data = tb_mmr_25_norm, outlier.shape = NA, width = .15, alpha = 0.5, position = position_nudge(x = .25)) +
  gghalves::geom_half_point(aes(x = Treatment, y = MMR_mass_25_norm, color = "25°C, normoxia", fill = "25°C, normoxia"), data = tb_mmr_25_norm, side = "l", range_scale = .4, alpha = .3, position = position_nudge(x = .25)) +
  
  geom_boxplot(aes(x = Treatment, y = MMR_mass_17_hyper, fill = "17°C, hyperoxia"), data = tb_mmr_17_hyper, outlier.shape = NA, width = .15, alpha = 0.5, position = position_nudge(x = -.25)) +
  gghalves::geom_half_point(aes(x = Treatment, y = MMR_mass_17_hyper, color = "17°C, hyperoxia", fill = "17°C, hyperoxia"), data = tb_mmr_17_hyper, side = "l", range_scale = .4, alpha = .3, position = position_nudge(x = -.25)) +
  
  geom_boxplot(aes(x = Treatment, y = MMR_mass_25_hyper, fill = "25°C, hyperoxia"), data = tb_mmr_25_hyper, outlier.shape = NA, width = .15, alpha = 0.5, position = position_nudge(x = -.25)) +
  gghalves::geom_half_point(aes(x = Treatment, y = MMR_mass_25_hyper, color = "25°C, hyperoxia", fill = "25°C, hyperoxia"), data = tb_mmr_25_hyper, side = "l", range_scale = .4, alpha = .3, position = position_nudge(x = -.25)) +
  
  geom_boxplot(aes(x = Treatment, y = MMR_mass_21_norm, fill = "21°C, normoxia"), data = tb_mmr_21_norm, outlier.shape = NA, width = .15, alpha = 0.5) +
  gghalves::geom_half_point(aes(x = Treatment, y = MMR_mass_21_norm, color = "21°C, normoxia", fill = "21°C, normoxia"), data = tb_mmr_21_norm, side = "l", range_scale = .4, alpha = .3) +
  
  labs(x = "Treatment", y = expression(MMR~(mg~O[2]~min^-1~88.5~g^-1))) +
  theme_classic() +
  scale_x_discrete(limits = c("17°C, normoxia", "17°C, hyperoxia", "21°C, normoxia", "25°C, normoxia", "25°C, hyperoxia")) +
  scale_fill_manual(name = "Treatment", values = c(
    "17°C, normoxia" = "slateblue1",
    "25°C, normoxia" = "darkorange2",
    "17°C, hyperoxia" = "blue",
    "25°C, hyperoxia" = "red3",
    "21°C, normoxia" = "orange"),
    breaks = c("17°C, normoxia", "17°C, hyperoxia", "21°C, normoxia", "25°C, normoxia", "25°C, hyperoxia")) +
  scale_color_manual(name = "Treatment", values = c(
    "17°C, normoxia" = "slateblue1",
    "25°C, normoxia" = "darkorange2",
    "17°C, hyperoxia" = "blue",
    "25°C, hyperoxia" = "red3",
    "21°C, normoxia" = "orange"),
    breaks = c("17°C, normoxia", "17°C, hyperoxia", "21°C, normoxia", "25°C, normoxia", "25°C, hyperoxia"))


#### Mass-standardised AAS calculation ####

lm_log_aas_17_norm <- lm(log10(AAS) ~ log10(mass), data = filter(tb_MR_master, Treatment == "17°C, normoxia"))
lm_log_aas_25_norm <- lm(log10(AAS) ~ log10(mass), data = filter(tb_MR_master, Treatment == "25°C, normoxia"))
lm_log_aas_17_hyper <- lm(log10(AAS) ~ log10(mass), data = filter(tb_MR_master, Treatment == "17°C, hyperoxia"))
lm_log_aas_25_hyper <- lm(log10(AAS) ~ log10(mass), data = filter(tb_MR_master, Treatment == "25°C, hyperoxia"))
lm_log_aas_21_norm <- lm(log10(AAS) ~ log10(mass), data = filter(tb_MR_master, Treatment == "21°C, normoxia"))

mass <- data.frame(mass = mean(tb_MR_master$mass))

AAS_mean_predicted_17_norm <- predict(lm_log_aas_17_norm, newdata = mass)
AAS_mean_predicted_25_norm <- predict(lm_log_aas_25_norm, newdata = mass)
AAS_mean_predicted_17_hyper <- predict(lm_log_aas_17_hyper, newdata = mass)
AAS_mean_predicted_25_hyper <- predict(lm_log_aas_25_hyper, newdata = mass)
AAS_mean_predicted_21_norm <- predict(lm_log_aas_21_norm, newdata = mass)

tb_MR_master_17_norm <- tb_MR_master[which(tb_MR_master$Treatment == "17°C, normoxia"),] %>% 
  add_residuals(lm_log_aas_17_norm)

tb_MR_master_17_norm <- 
  tb_MR_master_17_norm %>%
  mutate(
    AAS_mass_17_norm = 10^(AAS_mean_predicted_17_norm + resid),
    resid = NULL
  )

tb_MR_master_25_norm <- tb_MR_master[which(tb_MR_master$Treatment == "25°C, normoxia"),] %>% 
  add_residuals(lm_log_aas_25_norm)

tb_MR_master_25_norm <- 
  tb_MR_master_25_norm %>%
  mutate(
    AAS_mass_25_norm = 10^(AAS_mean_predicted_25_norm + resid),
    resid = NULL
  )

tb_MR_master_17_hyper <- tb_MR_master[which(tb_MR_master$Treatment == "17°C, hyperoxia"),] %>% 
  add_residuals(lm_log_aas_17_hyper)

tb_MR_master_17_hyper <- 
  tb_MR_master_17_hyper %>%
  mutate(
    AAS_mass_17_hyper = 10^(AAS_mean_predicted_17_hyper + resid),
    resid = NULL
  )

tb_MR_master_25_hyper <- tb_MR_master[which(tb_MR_master$Treatment == "25°C, hyperoxia"),] %>% 
  add_residuals(lm_log_aas_25_hyper)

tb_MR_master_25_hyper <- 
  tb_MR_master_25_hyper %>%
  mutate(
    AAS_mass_25_hyper = 10^(AAS_mean_predicted_25_hyper + resid),
    resid = NULL
  )

tb_MR_master_21_norm <- tb_MR_master[which(tb_MR_master$Treatment == "21°C, normoxia"),] %>% 
  add_residuals(lm_log_aas_21_norm)

tb_MR_master_21_norm <- 
  tb_MR_master_21_norm %>%
  mutate(
    AAS_mass_21_norm = 10^(AAS_mean_predicted_21_norm + resid),
    resid = NULL
  )


#### Mass-standardised AAS graph ####

p_AAS_mass <- ggplot() +
  geom_boxplot(aes(x = Treatment, y = AAS_mass_17_norm, fill = "17°C, normoxia"), data = tb_MR_master_17_norm, outlier.shape = NA, width = .15, alpha = 0.7, position = position_nudge(x = .25)) +
  gghalves::geom_half_point(aes(x = Treatment, y = AAS_mass_17_norm, color = "17°C, normoxia", fill = "17°C, normoxia"), data = tb_MR_master_17_norm, side = "l", range_scale = .4, alpha = .3, position = position_nudge(x = .25)) +
  
  geom_boxplot(aes(x = Treatment, y = AAS_mass_25_norm, fill = "25°C, normoxia"), data = tb_MR_master_25_norm, outlier.shape = NA, width = .15, alpha = 0.5, position = position_nudge(x = .25)) +
  gghalves::geom_half_point(aes(x = Treatment, y = AAS_mass_25_norm, color = "25°C, normoxia", fill = "25°C, normoxia"), data = tb_MR_master_25_norm, side = "l", range_scale = .4, alpha = .3, position = position_nudge(x = .25)) +
  
  geom_boxplot(aes(x = Treatment, y = AAS_mass_17_hyper, fill = "17°C, hyperoxia"), data = tb_MR_master_17_hyper, outlier.shape = NA, width = .15, alpha = 0.5, position = position_nudge(x = -.25)) +
  gghalves::geom_half_point(aes(x = Treatment, y = AAS_mass_17_hyper, color = "17°C, hyperoxia", fill = "17°C, hyperoxia"), data = tb_MR_master_17_hyper, side = "l", range_scale = .4, alpha = .3, position = position_nudge(x = -.25)) +
  
  geom_boxplot(aes(x = Treatment, y = AAS_mass_25_hyper, fill = "25°C, hyperoxia"), data = tb_MR_master_25_hyper, outlier.shape = NA, width = .15, alpha = 0.5, position = position_nudge(x = -.25)) +
  gghalves::geom_half_point(aes(x = Treatment, y = AAS_mass_25_hyper, color = "25°C, hyperoxia", fill = "25°C, hyperoxia"), data = tb_MR_master_25_hyper, side = "l", range_scale = .4, alpha = .3, position = position_nudge(x = -.25)) +
  
  geom_boxplot(aes(x = Treatment, y = AAS_mass_21_norm, fill = "21°C, normoxia"), data = tb_MR_master_21_norm, outlier.shape = NA, width = .15, alpha = 0.5) +
  gghalves::geom_half_point(aes(x = Treatment, y = AAS_mass_21_norm, color = "21°C, normoxia", fill = "21°C, normoxia"), data = tb_MR_master_21_norm, side = "l", range_scale = .4, alpha = .3) +
  
  labs(y = expression(AAS~(mg~O[2]~min^-1~88.5~g^-1))) +
  theme_classic() +
  scale_x_discrete(limits = c("17°C, normoxia", "17°C, hyperoxia", "21°C, normoxia", "25°C, normoxia", "25°C, hyperoxia")) +
  scale_fill_manual(name = "Treatment", values = c(
    "17°C, normoxia" = "slateblue1",
    "25°C, normoxia" = "darkorange2",
    "17°C, hyperoxia" = "blue",
    "25°C, hyperoxia" = "red3",
    "21°C, normoxia" = "orange"),
    breaks = c("17°C, normoxia", "17°C, hyperoxia", "21°C, normoxia", "25°C, normoxia", "25°C, hyperoxia")) +
  scale_color_manual(name = "Treatment", values = c(
    "17°C, normoxia" = "slateblue1",
    "25°C, normoxia" = "darkorange2",
    "17°C, hyperoxia" = "blue",
    "25°C, hyperoxia" = "red3",
    "21°C, normoxia" = "orange"),
    breaks = c("17°C, normoxia", "17°C, hyperoxia", "21°C, normoxia", "25°C, normoxia", "25°C, hyperoxia"))



#### Final mass-standardised graph ####

p_final <- ggarrange(p_SMR_mass + rremove("xlab") , 
                     p_MMR_mass + rremove("xlab"), 
                     p_AAS_mass + rremove("xlab"),
                     nrow = 3,
                     align = "v",
                     common.legend = TRUE,
                     heights = c(1, 1, 1))

plot(p_final)


#tb_mass_std <- tb_smr_17_norm %>% 
  full_join(tb_smr_17_hyper, by = "ID_fish") %>%
  full_join(tb_smr_25_norm, by = "ID_fish") %>%
  full_join(tb_smr_25_hyper, by = "ID_fish") %>%
  full_join(tb_smr_21_norm, by = "ID_fish") %>%
  mutate(SMR_mass = coalesce(SMR_mass_17_norm, SMR_mass_17_hyper),
         SMR_mass_17_norm = NULL,
         SMR_mass_17_hyper = NULL
  )


