library(dplyr)
library(tidyr)
library(ggpubr)
library(corrplot)
library(mgcv)
library(DHARMa)
library(mgcViz)
library(gridExtra)
library(ROCR)
library(recdata)
library(predRecruit)
library(dplyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(corrplot)
library(car)
library(gratia)
library(ggpubr)

#### Read In Data ####
marginals<-readRDS("Output/Data/marginals.rds")
nulls<-marginals%>%
  select(RMSEnull,species)%>%
  unique()
ps <- readRDS("Output/Data/ps_model_fits.rds")
ps_results<-data.frame(ps[["LOO"]][["results"]])%>%
  select(-RMSE)%>%
  rename(rsq=rsq_full, RMSE=RMSE_loo)%>%
  mutate(RMSE_type="LOO")%>%
add_row(data.frame(ps[["LFO5"]][["results"]])
        %>%select(-n_pred)%>%
          mutate(RMSE_type="LFO 5"))%>%
add_row(data.frame(ps[["LFO10"]][["results"]])%>%
          select(-n_pred)%>%
          mutate(RMSE_type="LFO 10"))

sb <- readRDS("Output/Data/sb_model_fits.rds")
sb_results<-data.frame(sb[["LOO"]][["results"]])%>%select(-RMSE)%>%
rename(rsq=rsq_full, RMSE=RMSE_loo)%>%
  mutate(RMSE_type="LOO")%>%
  add_row(data.frame(sb[["LFO5"]][["results"]])
          %>%select(-n_pred)%>%
            mutate(RMSE_type="LFO 5"))%>%
  add_row(data.frame(sb[["LFO10"]][["results"]])%>%
            select(-n_pred)%>%
            mutate(RMSE_type="LFO 10"))

yt <- readRDS("Output/Data/yt_model_fits.rds")
yt_results<-data.frame(yt[["LOO"]][["results"]])%>%#select(-RMSE)%>%
rename(rsq=rsq_full)%>%#, RMSE=RMSE_loo)%>%
  mutate(RMSE_type="LOO")%>%
  add_row(data.frame(yt[["LFO5"]][["results"]])
          %>%select(-n_pred)%>%
            mutate(RMSE_type="LFO 5"))%>%
  add_row(data.frame(yt[["LFO10"]][["results"]])%>%
            select(-n_pred)%>%
            mutate(RMSE_type="LFO 10"))

results<- data.frame(yt_results%>%
  add_row(sb_results)%>%
  add_row(ps_results))%>%
  merge(nulls)%>%
  mutate(rmse_null=RMSE/RMSEnull)
class(results$Hit)
results2<-results%>%select(species, ModelID,rmse_null ,RMSE_type)%>%
  pivot_wider(names_from = RMSE_type, values_from = rmse_null)
results3<-results%>%select(species, ModelID,Hit ,RMSE_type)%>%
  pivot_wider(names_from = RMSE_type, values_from = Hit)
#### Hit Metric Plot ####

ggplot(results, aes(x =RMSE_type, y =rmse_null)) +
  facet_wrap(~species)+
  geom_jitter(height = 0, width = 0.1)+
  geom_violin() + 
  geom_boxplot(width = 0.1)+
   geom_hline(yintercept=1, col='red', lty=2)+
  labs(x = "", y = "RMSE/RMSEnull")+
#  scale_fill_gradient(low = "white", high = "Darkgreen") +
  theme_classic()

ggplot(results, aes(x =RMSE_type, y =Hit)) +
  facet_wrap(~species)+
  geom_jitter(height = 0, width = 0.1)+
  geom_violin() + 
  geom_boxplot(width = 0.1)+
  geom_hline(yintercept=0.5, col='red', lty=2)+
  labs(x = "", y = "RMSE/RMSEnull")+
  #  scale_fill_gradient(low = "white", high = "Darkgreen") +
  theme_classic()

ggplot(results, aes(x =rsq, y =rmse_null)) +
  facet_grid(RMSE_type~species)+
  geom_point()+
  #geom_violin() + 
 # geom_boxplot(width = 0.1)+
  #geom_hline(yintercept=1, col='red', lty=2)+
  labs(x = "R-squared", y = "RMSE/RMSEnull")+
  #  scale_fill_gradient(low = "white", high = "Darkgreen") +
  theme_classic()

ggplot(results2) +
  facet_wrap(~species)+
  geom_hline(yintercept=1, col='red', lty=2)+
  geom_jitter(width = 0.025, height = 0.025, alpha = 0.5,aes(y=`LFO 10`, x=`LFO 5`))+
  labs(x = "", y = "RMSE/RMSEnull")+
  ylim(c(0.5,1.5))+
  geom_abline(intercept = 0, slope = 1, linetype = 5, color = "blue")+
  
  #  scale_fill_gradient(low = "white", high = "Darkgreen") +
  theme_classic()

ggplot(results2) +
  facet_wrap(~species)+
  geom_hline(yintercept=1, col='red', lty=2)+
  geom_jitter( height = 0.01, alpha = 0.5,aes(y=`LFO 10`, x=LOO))+
  labs(x = "", y = "RMSE/RMSEnull")+
  ylim(c(0.5,1.5))+
  #  scale_fill_gradient(low = "white", high = "Darkgreen") +
  theme_classic()

ggplot(results3) +
  facet_wrap(~species)+
  geom_hline(yintercept=0.5, col='red', lty=2)+
  geom_jitter(width = 0.01, height = 0.01, alpha = 0.5,aes(y=`LFO 10`, x=`LFO 5`))+
  labs(x = "Hit Rate LFO 5", y = "Hit Rate LFO 10")+
  geom_abline(intercept = 0, slope = 1, linetype = 5, color = "blue")+
  ylim(c(0,1))+
  xlim(c(0,1))+
  #  scale_fill_gradient(low = "white", high = "Darkgreen") +
  theme_classic()

ggplot(results3) +
  facet_wrap(~species)+
  geom_hline(yintercept=0.5, col='red', lty=2)+
  geom_jitter(width = 0.01, height = 0.01, alpha = 0.5,aes(y=`LFO 10`, x=LOO))+
  labs(x = "Hit Rate LFO 5", y = "Hit Rate LFO 10")+
  geom_abline(intercept = 0, slope = 1, linetype = 5, color = "blue")+
  ylim(c(0,1))+
  xlim(c(0,1))+
  #  scale_fill_gradient(low = "white", high = "Darkgreen") +
  theme_classic()
