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
tb_MR_master <- read.csv("tb_MR_master.csv") 
tb_data <- read.csv("tb_data.csv") 


#### Data Preparation ####

tb_MR_master$Treatment <-as.factor(tb_MR_master$Treatment)
tb_MR_master$Maturity <-as.factor(tb_MR_master$Maturity)

tb_MR_master$Treatment <- factor(tb_MR_master$Treatment, levels = c("17°C, normoxia", "17°C, hyperoxia", "21°C, normoxia", "25°C, normoxia", "25°C, hyperoxia"))


#### Metabolic allometry graph ####

p_SMR_17<- ggplot(data=tb_MR_master, aes(mass, SMR, shape = Maturity,  colour = Treatment, fill = Treatment )) +
  geom_point(data = subset(tb_MR_master, Treatment == "17°C, normoxia"), size=3, alpha =0.5) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "17°C, normoxia"), aes(group = Treatment), method = lm) +
  geom_point(data = subset(tb_MR_master, Treatment == "17°C, hyperoxia"), size=3, alpha =0.5) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "17°C, hyperoxia"), aes(group = Treatment), method = lm) + 
  geom_point(data = subset(tb_MR_master, Treatment == "21°C, normoxia"), size=3, alpha =0.3) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "21°C, normoxia"), aes(group = Treatment), alpha= 0.3, method = lm) +
  scale_x_log10(limits = c(10, 359.2), breaks = c(10, 20, 30, 40, 50, 100, 200, 300, 400, 500)) +
  scale_y_log10(limits = c(0.02, 1.5), breaks = c(0.02,0.03,0.04, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 1, 1.5)) +
  labs(x=NULL, y = expression(SMR~(mg~O[2]~min^-1))) +  
  theme_classic()+
  scale_color_manual(values = c("slateblue1", "blue","grey")) +
  scale_fill_manual(values = c("slateblue1", "blue","grey")) +
  theme(axis.text = element_text(size = 14, colour = "black"), 
        axis.title = element_text(size =20, colour = "black"))

# Add the equation text
eq_SMR_norm_17 <- expression(SMR[100] ~ "=" ~ 0.0047 ~ M[B] ~ phantom()^0.76)

eq_SMR_hyper_17 <- expression(SMR[150] ~ "=" ~ 0.0043 ~ M[B] ~ phantom()^0.77)

p_SMR_17 <- p_SMR_17 +
  annotate("text", x = 25, y = 1.5, label = eq_SMR_norm_17, size = 5) +
  annotate("text", x = 25, y = 1.1, label = eq_SMR_hyper_17, size = 5)


p_SMR_25 <- ggplot(data=tb_MR_master, aes(mass, SMR, shape = Maturity, colour = Treatment, fill = Treatment )) +
  geom_point(data = subset(tb_MR_master, Treatment == "25°C, normoxia"), size=3, alpha =0.5) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "25°C, normoxia"), aes(group = Treatment), method = lm) +
  geom_point(data = subset(tb_MR_master, Treatment == "25°C, hyperoxia"), size=3, alpha =0.5) +
  geom_point(data = subset(tb_MR_master, Treatment == "21°C, normoxia"), size=3, alpha =0.3) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "21°C, normoxia"), aes(group = Treatment), alpha= 0.3, method = lm) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "25°C, hyperoxia"), aes(group = Treatment), method = lm) +
  scale_x_log10(limits = c(10, 359.2), breaks = c(10, 20, 30, 40, 50, 100, 200, 300, 400, 500)) +
  scale_y_log10(limits = c(0.02, 1.5), breaks = c(0.02,0.03,0.04, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 1, 1.5)) +
  labs(x=NULL, y=NULL) + 
  theme_classic() +
  scale_color_manual(values = c( "darkorange2", "red3","grey")) +
  scale_fill_manual(values = c("darkorange2", "red3","grey")) +
  theme(axis.text = element_text(size = 14, colour = "black"), 
        axis.title = element_text(size =20, colour = "black"))

eq_SMR_norm_25 <- expression(SMR[100] ~ "=" ~ 0.0058 ~ M[B] ~ phantom()^"0.90")

eq_SMR_hyper_25 <- expression(SMR[150] ~ "=" ~ 0.0091 ~ M[B] ~ phantom()^0.78)

p_SMR_25 <- p_SMR_25 +
  annotate("text", x = 25, y = 1.5, label = eq_SMR_norm_25, size = 5) +
  annotate("text", x = 25, y = 1.1, label = eq_SMR_hyper_25, size = 5)

p_MMR_17<- ggplot(data=tb_MR_master, aes(mass, MMR,shape = Maturity, colour = Treatment, fill = Treatment )) +
  geom_point(data = subset(tb_MR_master, Treatment == "17°C, normoxia"), size=3, alpha =0.5) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "17°C, normoxia"), aes(group = Treatment), method = lm) +
  geom_point(data = subset(tb_MR_master, Treatment == "17°C, hyperoxia"), size=3, alpha =0.5) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "17°C, hyperoxia"), aes(group = Treatment), method = lm) +
  geom_point(data = subset(tb_MR_master, Treatment == "21°C, normoxia"), size=3, alpha =0.3) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "21°C, normoxia"), aes(group = Treatment), alpha= 0.3, method = lm) +
  scale_x_log10(limits = c(10, 359.2), breaks = c(10, 20, 30, 40, 50, 100, 200, 300, 400, 500)) +
  scale_y_log10(limits = c(0.14, 5), breaks = c(0.1, 0.2,0.3, 0.4, 0.5, 1, 2, 3, 5)) +
  labs(x=NULL, y = expression(MMR~(mg~O[2]~min^-1))) +  
  theme_classic() +
  scale_color_manual(values = c("slateblue1", "blue","grey")) +
  scale_fill_manual(values = c("slateblue1", "blue","grey")) +
  theme(axis.text = element_text(size = 14, colour = "black"), 
        axis.title = element_text(size =20, colour = "black"))

eq_MMR_norm_17 <- expression(MMR[100] ~ "=" ~ 0.0222 ~ M[B] ~ phantom()^0.86)

eq_MMR_hyper_17 <- expression(MMR[150] ~ "=" ~ 0.0257 ~ M[B] ~ phantom()^0.86)

p_MMR_17 <- p_MMR_17 +
  annotate("text", x = 25, y = 5.0, label = eq_MMR_norm_17, size = 5) +
  annotate("text", x = 25, y = 3.8, label = eq_MMR_hyper_17, size = 5)

p_MMR_25 <- ggplot(data = tb_MR_master, aes(mass, MMR, shape = Maturity, colour = Treatment, fill = Treatment)) +
  geom_point(data = subset(tb_MR_master, Treatment == "25°C, normoxia"), size=3, alpha =0.5) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "25°C, normoxia"), aes(group = Treatment), method = lm) +
  geom_point(data = subset(tb_MR_master, Treatment == "25°C, hyperoxia"), size=3, alpha =0.5) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "25°C, hyperoxia"), aes(group = Treatment), method = lm) +
  geom_point(data = subset(tb_MR_master, Treatment == "21°C, normoxia"), size=3, alpha =0.3) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "21°C, normoxia"), aes(group = Treatment), alpha= 0.3, method = lm) +
  scale_x_log10(limits = c(10, 359.2), breaks = c(10, 20, 30, 40, 50, 100, 200, 300, 400, 500)) +
  scale_y_log10(limits = c(0.14, 5), breaks = c(0.1, 0.2,0.3, 0.4, 0.5, 1, 2, 3, 5)) +
  labs(x = NULL, y = NULL) + 
  theme_classic() +
  scale_color_manual(values = c("darkorange2", "red3","grey")) +
  scale_fill_manual(values = c("darkorange2", "red3","grey")) +
  theme(axis.text = element_text(size = 14, colour = "black"), 
        axis.title = element_text(size =20, colour = "black"))

eq_MMR_norm_25 <- expression(MMR[100] ~ "=" ~ 0.0302 ~ M[B] ~ phantom()^0.79)

eq_MMR_hyper_25 <- expression(MMR[150] ~ "=" ~ "0.0320" ~ M[B] ~ phantom()^0.84)

p_MMR_25 <- p_MMR_25 +
  annotate("text", x = 25, y = 5.0, label = eq_MMR_norm_25, size = 5) +
  annotate("text", x = 25, y = 3.8, label = eq_MMR_hyper_25, size = 5)

p_AAS_17 <- ggplot(data=tb_MR_master, aes(mass, AAS, shape = Maturity, colour = Treatment, fill = Treatment)) +
  geom_point(data = subset(tb_MR_master, Treatment == "17°C, normoxia"), size=3, alpha =0.5) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "17°C, normoxia"), aes(group = Treatment), method = lm) +
  geom_point(data = subset(tb_MR_master, Treatment == "17°C, hyperoxia"), size=3, alpha =0.5) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "17°C, hyperoxia"), aes(group = Treatment), method = lm) +
  geom_point(data = subset(tb_MR_master, Treatment == "21°C, normoxia"), size=3, alpha =0.3) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "21°C, normoxia"), aes(group = Treatment), alpha= 0.3, method = lm) +
  scale_x_log10(limits = c(10, 359.2), breaks = c(10, 20, 30, 40, 50, 100, 200, 300, 400, 500)) +
  scale_y_log10(limits = c(0.14, 5), breaks = c(0.1, 0.2,0.3, 0.4, 0.5, 1, 2, 3, 5)) +
  labs(x=NULL, y = expression(AAS~(mg~O[2]~min^-1))) +  
  theme_classic() +
  scale_color_manual(values = c("slateblue1", "blue", "grey")) +
  scale_fill_manual(values = c("slateblue1", "blue", "grey")) +
  theme(axis.text = element_text(size = 14, colour = "black"), 
        axis.title = element_text(size =20, colour = "black"))

eq_AAS_norm_17 <- expression(AAS[100] ~ "=" ~ 0.0177 ~ M[B] ~ phantom()^0.87)

eq_AAS_hyper_17 <- expression(AAS[150] ~ "=" ~ 0.0217 ~ M[B] ~ phantom()^0.87)

p_AAS_17 <- p_AAS_17 +
  annotate("text", x = 25, y = 5.0, label = eq_AAS_norm_17, size = 5) +
  annotate("text", x = 25, y = 3.8, label = eq_AAS_hyper_17, size = 5)

p_AAS_25 <- ggplot(data = tb_MR_master, aes(mass, AAS, shape = Maturity, colour = as.factor(Treatment), fill = as.factor(Treatment))) +
  geom_point(data = subset(tb_MR_master, Treatment == "25°C, normoxia"), size=3, alpha =0.5) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "25°C, normoxia"), aes(group = Treatment), method = lm) +
  geom_point(data = subset(tb_MR_master, Treatment == "25°C, hyperoxia"), size=3, alpha =0.5) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "25°C, hyperoxia"), aes(group = Treatment), method = lm) +
  geom_point(data = subset(tb_MR_master, Treatment == "21°C, normoxia"), size=3, alpha =0.3) +
  geom_smooth(data = subset(tb_MR_master, Treatment == "21°C, normoxia"), aes(group = Treatment), alpha= 0.3, method = lm) +
  scale_x_log10(limits = c(10, 359.2), breaks = c(10, 20, 30, 40, 50, 100, 200, 300, 400, 500)) +
  scale_y_log10(limits = c(0.14, 5), breaks = c(0.1, 0.2,0.3, 0.4, 0.5, 1, 2, 3, 5)) +
  labs(x = NULL, y = NULL) +  
  theme_classic() +
  scale_color_manual(values = c("darkorange2", "red3","grey50")) +
  scale_fill_manual(values = c("darkorange2", "red3","grey50")) +
  theme(axis.text = element_text(size = 14, colour = "black"), 
        axis.title = element_text(size =20, colour = "black"))

eq_AAS_norm_25 <- expression(AAS[100] ~ "=" ~ 0.0258 ~ M[B] ~ phantom()^0.74)

eq_AAS_hyper_25 <- expression(AAS[150] ~ "=" ~ 0.0227 ~ M[B] ~ phantom()^0.86)

p_AAS_25 <- p_AAS_25 +
  annotate("text", x = 25, y = 5.0, label = eq_AAS_norm_25, size = 5) +
  annotate("text", x = 25, y = 3.8, label = eq_AAS_hyper_25, size = 5)


#### Final Plot #### 
p_final <- ggarrange(p_SMR_17, 
                     p_SMR_25, 
                     p_MMR_17,
                     p_MMR_25, 
                     p_AAS_17, 
                     p_AAS_25, 
                     nrow = 3, 
                     ncol = 2,
                     align = "v", label.y = 1, label.x = 0.1, 
                     common.legend = TRUE 
                     )


p_final <- p_SMR_17 + p_SMR_25 + p_MMR_17 + p_MMR_25 + p_AAS_17 + p_AAS_25 + 
  plot_layout(ncol = 2, guides = "collect") 

plot(p_final)
