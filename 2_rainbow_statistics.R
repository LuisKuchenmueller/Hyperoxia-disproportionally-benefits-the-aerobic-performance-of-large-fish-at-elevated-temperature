#### Workspace ####

rm(list=ls())
setwd("C:/Users/s222165681/OneDrive - Deakin University/Rainbow trout allometry/Experiment 1/Results/Processed") 


#### Packages ####

library(ggplot2)
library(emmeans)
library(lme4)
library(multcomp)
library(multcompView)
library(lsmeans)
library(tidyverse)
library(car)
library(AICcmodavg)


#### Import Processed Data ####

tb_MR_master <- read.csv("tb_MR_master.csv") 


#### Data Preparation ####

grand_mean <- mean(tb_MR_master$mass)
log_grand_mean <- mean(tb_MR_master$log_mass)

tb_MR_master$Temp_class <-as.factor(tb_MR_master$Temp_class) #sets as factor
tb_MR_master$O2 <-as.factor(tb_MR_master$O2)


#### Test Assumptions ####

qqnorm(log10(tb_MR_master$mass))
qqline(log10(tb_MR_master$mass))
hist(log10(tb_MR_master$mass))

qqnorm(log10(tb_MR_master$SMR))
qqline(log10(tb_MR_master$SMR))
hist(log10(tb_MR_master$SMR))

qqnorm(log10(tb_MR_master$MMR))
qqline(log10(tb_MR_master$MMR))
hist(log10(tb_MR_master$MMR))

qqnorm(log10(tb_MR_master$AAS))
qqline(log10(tb_MR_master$AAS))
hist(log10(tb_MR_master$AAS))


#### SMR Linear Model ####

model_SMR <- lm(
  formula = log10(SMR)~ log_mass*Temp_class*O2,
  data = tb_MR_master
)

model_SMR %>% summary()
model_SMR %>% confint()
model_SMR %>% sjPlot::plot_model()
model_SMR %>% residuals() %>% qqnorm()
model_SMR %>% plot()
model_SMR %>% Anova(contrasts=list(log_mass=contr.sum, Temp_class=contr.sum, O2=contr.sum), type=3)
MuMIn::r.squaredGLMM(model_SMR)
hist(model_SMR$residuals)


#### SMR Slope and Intercept Comparison ####

SMR_slope_compare <- lstrends (model_SMR, ~ Temp_class*O2, var = "log_mass")

SMR_slope <- cld(SMR_slope_compare)

SMR_int_compare_grand_mean <- emmeans(model_SMR, ~ Temp_class*O2, at = list(log_mass = log_grand_mean))

SMR_int_grand_mean <- cld(SMR_int_compare_grand_mean)

SMR_compare_factor <- emmeans(model_SMR, ~ Temp_class*O2, at = list(log_mass = 0))

SMR_factor <- cld(SMR_compare_factor)


#### Extract SMR table values ####

tb_SMR_slope <- as.data.frame(SMR_slope)  %>% 
  mutate(
    Treatment = paste0(Temp_class, "°C, ", O2),
    SMR_slope = log_mass.trend,
    SMR_slope_SE = SE,
    SMR_slope_df = df,
    SMR_slope_lower_CL = lower.CL,
    SMR_slope_upper_CL = upper.CL,
    SMR_slope_group = .group,
    Temp_class = NULL,
    O2 = NULL,
    log_mass.trend = NULL,
    SE = NULL,
    df = NULL,
    lower.CL = NULL,
    upper.CL = NULL,
    .group = NULL
  ) %>% 
  filter(
    Treatment != "21°C, hyperoxia"
  )

tb_SMR_factor <- as.data.frame(SMR_factor)  %>% 
  mutate(
    Treatment = paste0(Temp_class, "°C, ", O2),
    SMR_factor = 10^emmean,
    #SMR_int_grand_mean_SE = 10^SE,
    SMR_factor_df = df,
    SMR_factor_lower_CL = 10^lower.CL,
    SMR_factor_upper_CL = 10^upper.CL,
    SMR_factor_group = .group,
    Temp_class = NULL,
    O2 = NULL,
    emmean = NULL,
    SE = NULL,
    df = NULL,
    lower.CL = NULL,
    upper.CL = NULL,
    .group = NULL
  ) %>% 
  filter(
    Treatment != "21°C, hyperoxia"
  )

tb_SMR_int_grand_mean <- as.data.frame(SMR_int_grand_mean)  %>% 
  mutate(
    Treatment = paste0(Temp_class, "°C, ", O2),
    SMR_int_grand_mean = 10^emmean,
    #SMR_int_grand_mean_SE = 10^SE,
    SMR_int_grand_mean_df = df,
    SMR_int_grand_mean_lower_CL = 10^lower.CL,
    SMR_int_grand_mean_upper_CL = 10^upper.CL,
    SMR_int_grand_mean_group = .group,
    Temp_class = NULL,
    O2 = NULL,
    emmean = NULL,
    SE = NULL,
    df = NULL,
    lower.CL = NULL,
    upper.CL = NULL,
    .group = NULL
  ) %>% 
  filter(
    Treatment != "21°C, hyperoxia"
    )


#### MMR Linear Model ####

model_MMR <- lm(
  formula = log10(MMR)~ log_mass*Temp_class*O2,
  data = tb_MR_master
)

model_MMR %>% summary()
model_MMR %>% confint()
model_MMR %>% sjPlot::plot_model()
model_MMR %>% residuals() %>% qqnorm()
model_MMR %>% plot()
model_MMR %>% Anova(contrasts=list(log_mass=contr.sum, Temp_class=contr.sum, O2=contr.sum), type=3)
MuMIn::r.squaredGLMM(model_MMR)
hist(model_MMR$residuals)


#### MMR Slope and Intercept Comparison ####

MMR_compare <- lstrends (model_MMR_1, ~ Temp_class*O2, var = "log_mass")

MMR_slope <- cld(MMR_compare)

MMR_int_compare_grand_mean <- emmeans(model_MMR_1, ~ Temp_class*O2, at = list(log_mass = log_grand_mean))

MMR_int_grand_mean <- cld(MMR_int_compare_grand_mean)

MMR_compare_factor <- emmeans(model_MMR_1, ~ Temp_class*O2, at = list(log_mass = 0))

MMR_factor <- cld(MMR_compare_factor)

MMR_compare_factor <- emmeans(model_MMR_1, ~ Temp_class*O2, at = list(log_mass = 0))

MMR_factor <- cld(MMR_compare_factor)


#### Extract MMR table values ####

tb_MMR_slope <- as.data.frame(MMR_slope)  %>% 
  mutate(
    Treatment = paste0(Temp_class, "°C, ", O2),
    MMR_slope = log_mass.trend,
    MMR_slope_SE = SE,
    MMR_slope_df = df,
    MMR_slope_lower_CL = lower.CL,
    MMR_slope_upper_CL = upper.CL,
    MMR_slope_group = .group,
    Temp_class = NULL,
    O2 = NULL,
    log_mass.trend = NULL,
    SE = NULL,
    df = NULL,
    lower.CL = NULL,
    upper.CL = NULL,
    .group = NULL
  ) %>% 
  filter(
    Treatment != "21°C, hyperoxia"
  )

tb_MMR_factor <- as.data.frame(MMR_factor)  %>% 
  mutate(
    Treatment = paste0(Temp_class, "°C, ", O2),
    MMR_factor = 10^emmean,
    #SMR_int_grand_mean_SE = 10^SE,
    MMR_factor_df = df,
    MMR_factor_lower_CL = 10^lower.CL,
    MMR_factor_upper_CL = 10^upper.CL,
    MMR_factor_group = .group,
    Temp_class = NULL,
    O2 = NULL,
    emmean = NULL,
    SE = NULL,
    df = NULL,
    lower.CL = NULL,
    upper.CL = NULL,
    .group = NULL
  ) %>% 
  filter(
    Treatment != "21°C, hyperoxia"
  )

tb_MMR_int_grand_mean <- as.data.frame(MMR_int_grand_mean)  %>% 
  mutate(
    Treatment = paste0(Temp_class, "°C, ", O2),
    MMR_int_grand_mean = 10^emmean,
    #MMR_int_grand_mean_SE = 10^SE,
    MMR_int_grand_mean_df = df,
    MMR_int_grand_mean_lower_CL = 10^lower.CL,
    MMR_int_grand_mean_upper_CL = 10^upper.CL,
    MMR_int_grand_mean_group = .group,
    Temp_class = NULL,
    O2 = NULL,
    emmean = NULL,
    SE = NULL,
    df = NULL,
    lower.CL = NULL,
    upper.CL = NULL,
    .group = NULL
  ) %>% 
  filter(
    Treatment != "21°C, hyperoxia"
  )

#### AAS Linear Model ####

model_AAS <- lm(
  formula = log10(AAS)~ log_mass*Temp_class*O2,
  data = tb_MR_master
)

model_AAS %>% summary()
model_AAS %>% confint()
model_AAS %>% sjPlot::plot_model()
model_AAS %>% residuals() %>% qqnorm()
model_AAS %>% plot()
model_AAS %>% Anova(contrasts=list(log_mass=contr.sum, Temp_class=contr.sum, O2=contr.sum), type=3)
MuMIn::r.squaredGLMM(model_AAS)
hist(model_AAS$residuals)


#### AAS Slope and Intercept Comparison ####

AAS_compare <- lstrends (model_AAS_1, ~ Temp_class*O2, var = "log_mass")

AAS_slope <- cld(AAS_compare)

AAS_int_compare_grand_mean <- emmeans(model_AAS_1, ~ Temp_class*O2, at = list(log_mass = log_grand_mean))

AAS_int_grand_mean <- cld(AAS_int_compare_grand_mean)

AAS_compare_factor <- emmeans(model_AAS_1, ~ Temp_class*O2, at = list(log_mass = 0))

AAS_factor <- cld(AAS_compare_factor)


#### Extract AAS table values ####

tb_AAS_slope <- as.data.frame(AAS_slope)  %>% 
  mutate(
    Treatment = paste0(Temp_class, "°C, ", O2),
    AAS_slope = log_mass.trend,
    AAS_slope_SE = SE,
    AAS_slope_df = df,
    AAS_slope_lower_CL = lower.CL,
    AAS_slope_upper_CL = upper.CL,
    AAS_slope_group = .group,
    Temp_class = NULL,
    O2 = NULL,
    log_mass.trend = NULL,
    SE = NULL,
    df = NULL,
    lower.CL = NULL,
    upper.CL = NULL,
    .group = NULL
  ) %>% 
  filter(
    Treatment != "21°C, hyperoxia"
  )

tb_AAS_factor <- as.data.frame(AAS_factor)  %>% 
  mutate(
    Treatment = paste0(Temp_class, "°C, ", O2),
    AAS_factor = 10^emmean,
    #SMR_int_grand_mean_SE = 10^SE,
    AAS_factor_df = df,
    AAS_factor_lower_CL = 10^lower.CL,
    AAS_factor_upper_CL = 10^upper.CL,
    AAS_factor_group = .group,
    Temp_class = NULL,
    O2 = NULL,
    emmean = NULL,
    SE = NULL,
    df = NULL,
    lower.CL = NULL,
    upper.CL = NULL,
    .group = NULL
  ) %>% 
  filter(
    Treatment != "21°C, hyperoxia"
  )

tb_AAS_int_grand_mean <- as.data.frame(AAS_int_grand_mean)  %>% 
  mutate(
    Treatment = paste0(Temp_class, "°C, ", O2),
    AAS_int_grand_mean = 10^emmean,
    #AAS_int_grand_mean_SE = 10^SE,
    AAS_int_grand_mean_df = df,
    AAS_int_grand_mean_lower_CL = 10^lower.CL,
    AAS_int_grand_mean_upper_CL = 10^upper.CL,
    AAS_int_grand_mean_group = .group,
    Temp_class = NULL,
    O2 = NULL,
    emmean = NULL,
    SE = NULL,
    df = NULL,
    lower.CL = NULL,
    upper.CL = NULL,
    .group = NULL
  ) %>% 
  filter(
    Treatment != "21°C, hyperoxia"
  )


#### Combine table values ####

tb_data <- tb_SMR_slope %>%
  full_join( tb_SMR_int_grand_mean, by = "Treatment") %>%
  full_join( tb_MMR_slope, by = "Treatment") %>%
  full_join( tb_MMR_int_grand_mean, by = "Treatment") %>%
  full_join( tb_AAS_slope, by = "Treatment") %>%
  full_join( tb_AAS_int_grand_mean, by = "Treatment") %>%
  full_join( tb_SMR_factor, by = "Treatment") %>%
  full_join( tb_MMR_factor, by = "Treatment") %>%
  full_join( tb_AAS_factor, by = "Treatment") 
