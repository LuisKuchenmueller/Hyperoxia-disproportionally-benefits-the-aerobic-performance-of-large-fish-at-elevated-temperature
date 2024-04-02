#### Workspace ####

rm(list=ls())
setwd("XX") 


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


#### Functions ####

# calcSMR: computes standard metabolic rate (SMR)
# based on: Chabot, D., Steffensen, J. F. and Farrell, A. P. (2016). The determination of standard metabolic rate in fishes: measuring smr in fishes. J Fish Biol 88, 81-121.
# 
# mlnd:       'mean of the lowest normal distribution'
# tenpc:      number of 10th percentile of sorted oxygen slopes
# SD10pc:     standard deviation of the lowest 10% of oxygen slopes
# low10pc:    average oxygen slope of the lowest 10% of oxygen slopes; numerical handling of outliers, and subsequent shift of 'tenpc'
# low10pc_SD: as above; outliers defined by being larger than 2 SD from 'tenpc', and subsequent shift of 'tenpc'; note warning of multiple numerical expressions used expression is validated in script '6_rainbow_diagnostic_graphs.R'

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
  low10pc = mean(u[6:(5 + round(0.1*(length(u)-5)))])
  low10pc_SD = mean(u[(which((u > (mean(u[1:tenpc])- 2*SD10pc)))):(tenpc+which((u > (mean(u[1:tenpc])- 2*SD10pc-1))))])
  return(list(mlnd=mlnd, low10pc=low10pc, low10pc_SD))
}


#### Import respirometry data ####

tb_respirometry <- 
  read_excel("SMR.xlsx") %>%
  mutate(
    dateTime       = as.POSIXct(ymd_hms(dateTime)),
    Time           = as_hms(ymd_hms(dateTime)),
    Date           = as.Date(dateTime, format = "%Y/%m/%d"),
  ) %>%
  mutate(
    Temp_class_dup = Temp_class,
    O2_dup         = O2,
  ) %>%
  unite(
    Temp_class_dup, O2_dup, col = "Treatment", sep = "°C, "
  )

tb_mmr <- 
  read_excel("MMR.xlsx") %>%
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


#### Metabolic rate calculations ####

tb_respirometry <- tb_respirometry %>%
  mutate(
    volume_net = volume_ch - mass,
    MR_wBR     = abs(slope_wBR*(volume_net/1000)*60), 
    BR         = BRSlope*(volume_ch/1000)*60, 
    MR         = MR_wBR + BR,
  )

tb_smr <-
  tb_respirometry %>%
  group_by(
    ID_fish, mass, volume_net, ID_chamber, Treatment, Trial, O2, Temp_class
  ) %>% 
  summarise(
    SMR        = calcSMR(MR)$low10pc %>% unname(),
    SMR_mlnd    = calcSMR(MR)$mlnd %>% unname(),
    time_start = dateTime %>% min(),
    time_end   = dateTime %>% max(),
  ) 

tb_mmr <- tb_mmr %>%
  mutate(
    volume_net = volume_ch - mass,
    MMR_wBR    = abs(slope_wBR)*(volume_net/1000)*60, 
    BR         = BRSlope*(volume_ch/1000)*60,
    MMR        = MMR_wBR + BR
  ) 


#### Export data ####

  tb_MR_master <- full_join(tb_smr, tb_mmr %>% select(ID_fish, MMR, Maturity), by = c("ID_fish"), copy = FALSE,
                            keep = FALSE,
                            na_matches = "na")
  
  tb_MR_master <- tb_MR_master %>%
    add_column(
      log_mass = log10(tb_MR_master$mass),
      log_SMR  = log10(tb_MR_master$SMR),
      log_MMR  = log10(tb_MR_master$MMR),
      AAS      = tb_MR_master$MMR-tb_MR_master$SMR
      )
  
write.csv(tb_respirometry, file = "tb_respirometry.csv", row.names = FALSE) 
write.csv(tb_smr, file = "tb_smr.csv", row.names = FALSE)
write.csv(tb_mmr, file = "tb_mmr.csv", row.names = FALSE)
write.csv(tb_MR_master, file = "tb_MR_master.csv", row.names = FALSE)
