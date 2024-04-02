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


#### Functions ####

calcSMR = function(Y, q=c(0.1,0.15,0.2,0.25,0.3), G=1:4){
  u = sort(Y)
  the.Mclust <- Mclust(Y,G=G,verbose=F)
  cl <- the.Mclust$classification
  cl2 <- as.data.frame(table(cl))
  cl2$cl <- as.numeric(levels(cl2$cl))
  valid <- cl2$Freq>=0.1*length(time)
  the.cl <- min(cl2$cl[valid])
  left.distr <- Y[the.Mclust$classification==the.cl]
  mlnd = the.Mclust$parameters$mean[the.cl]
  tenpc <- round(0.1*length(u))
  SD10pc <- sd(u[1:tenpc])
  low10pc = mean(u[(which((u > (mean(u[1:tenpc])-SD10pc)))):(tenpc+which((u > (mean(u[1:tenpc])-SD10pc-1))))])
  return(list(mlnd=mlnd, low10pc=low10pc))
}

#### Import Raw Data ####

tb_respirometry_morts <- 
  read_excel("Morts.xlsx") %>%
  mutate(
    dateTime       = as.POSIXct(ymd_hms(dateTime)),
    Time           = as_hms(ymd_hms(dateTime)),
    Date           = as.Date(dateTime, format = "%Y/%m/%d"),
    Temp_class_dup = Temp_class,
    O2_dup         = O2,
  ) %>%
  unite(
    Temp_class_dup, O2_dup, col = "Treatment", sep = "°C, "
  )


#### Processing ####

tb_respirometry_morts <- tb_respirometry_morts %>%
  mutate(
    volume_net = volume_ch - mass,
    MR_wBR     = abs(slope_wBR*(volume_net/1000)*60), 
    BR         = BRSlope*(volume_ch/1000)*60, 
    MR         = MR_wBR + BR,
  )

tb_respirometry_n25m <- tb_respirometry_morts[which(tb_respirometry_morts$Treatment=="25°C, normoxia"),] %>%
  arrange(desc(mass))

tb_respirometry_h25m <- tb_respirometry_morts[which(tb_respirometry_morts$Treatment=="25°C, hyperoxia"),] %>%
  arrange(desc(mass))


#### Mortality graphs ####

p_n25_morts <- ggplot(data=tb_respirometry_n25m) +
  geom_point(aes(ch_time, MR)) +
  facet_wrap(~factor(mass),  scales = "free_y") +
  labs(x="Chamber time (h)", y = expression(Metabolic~rate~(mg~O[2]~min^-1)), title = "25°C, normoxia - Mortality" ) +
  theme_classic()

plot(p_n25_morts)

p_h25_morts <- ggplot(data=tb_respirometry_h25m) +
  geom_point(aes(ch_time, MR)) +
  facet_wrap(~factor(mass),  scales = "free_y") +
  labs(x="Chamber time (h)", y = expression(Metabolic~rate~(mg~O[2]~min^-1)), title = "25°C, hyperoxia - Mortality" ) +
  theme_classic()

plot(p_h25_morts)
