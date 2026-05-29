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
library(NatParksPalettes)
library(scales)

col_dv<-natparks.pals("DeathValley", 7)
names(NatParksPalettes)
#### Read In Data ####
ps <- readRDS("Output/Data/Analysis Part 1/ps_model_fits.rds")
ps_loo<-data.frame(ps[["LOO"]][["results"]])
ps_LFO5<-data.frame(ps[["LFO5"]][["results"]])
ps_LFO10<-data.frame(ps[["LFO10"]][["results"]])
ps_loo_pred<-data.frame(ps[["LOO"]][["predicted"]])
ps_LFO5_pred<-data.frame(ps[["LFO5"]][["predicted"]])
ps_LFO10_pred<-data.frame(ps[["LFO10"]][["predicted"]])
ps_dat <- data.frame(read.csv("Data/Petrale/DATA_Combined_glorys_petrale_STANDARDIZED.csv"))%>%
  select(-c(X, year.1))%>%
  # filter(type=="Main_RecrDev")%>%
  filter(year>1994&year<=2018)

sb <- readRDS("Output/Data/Analysis Part 1/sb_model_fits.rds")
sb_loo<-data.frame(sb[["LOO"]][["results"]])
sb_LFO5<-data.frame(sb[["LFO5"]][["results"]])
sb_LFO10<-data.frame(sb[["LFO10"]][["results"]])
sb_loo_pred<-data.frame(sb[["LOO"]][["predicted"]])
sb_LFO5_pred<-data.frame(sb[["LFO5"]][["predicted"]])
sb_LFO10_pred<-data.frame(sb[["LFO10"]][["predicted"]])
sb_dat <- data.frame(read.csv("Data/Sablefish/data-combined-glorys-sablefish_STANDARDIZED.csv"))%>%
  select(-c(X))%>%
  filter(type=="Main_RecrDev")%>%
  filter(year>1994&year<=2018)

yt <- readRDS("Output/Data/Analysis Part 1/yt_model_fits.rds")
yt_loo<-data.frame(yt[["LOO"]][["results"]])
yt_LFO5<-data.frame(yt[["LFO5"]][["results"]])
yt_LFO10<-data.frame(yt[["LFO10"]][["results"]])
yt_loo_pred<-data.frame(yt[["LOO"]][["predicted"]])
yt_LFO5_pred<-data.frame(yt[["LFO5"]][["predicted"]])
yt_LFO10_pred<-data.frame(yt[["LFO10"]][["predicted"]])
yt_dat <- data.frame(read.csv("Data/Yellowtail/yt_fulldataset_STANDARDIZED.csv"))%>%
  select(-c(X, hci2_pjuv,hci1_pjuv, lusi_annual, hci2_larv,hci1_larv))%>%
  filter(type=="Main")%>%
  filter(Datatreatment=="2025 Final"&year>1994&year<=2018)

subsethk<-readRDS("Output/Data/hakesubset.rds")
hk <- readRDS("Output/Data/Analysis Part 1/hk_model_fits.rds")
hk_loo<-data.frame(hk[["LOO"]][["results"]])
hk_LFO5<-data.frame(hk[["LFO5"]][["results"]])
hk_LFO10<-data.frame(hk[["LFO10"]][["results"]])
hk_loo_pred<-data.frame(hk[["LOO"]][["predicted"]])
hk_LFO5_pred<-data.frame(hk[["LFO5"]][["predicted"]])
hk_LFO10_pred<-data.frame(hk[["LFO10"]][["predicted"]])
hk_dat <- data.frame(read.csv("Data/Hake/DATA_Combined_glorys_hake_STANDARDIZED.csv"))%>%
  select(-X)%>%
  filter(type=="Main_RecrDev"&year>1993&year<=2018)%>%
  select(all_of(subsethk), Y_rec, year, type, sd)
z<- ggplot()+theme_void()

#### Generating a Combined Index Selection ####
loo_weight<-2
hk_combined<- hk_loo%>%select(ModelID, species, AIC, RMSE_loo, var1, var2, var3, var4)%>%
  left_join(hk_LFO5%>%select(ModelID, species, RMSE, var1, var2, var3, var4)%>%
  rename(RMSE_5=RMSE))%>%
  left_join(hk_LFO10%>%select(ModelID, species, RMSE, var1, var2, var3, var4)%>%
              rename(RMSE_10=RMSE))%>%
  mutate(RMSE_combined=RMSE_loo*loo_weight+RMSE_5)%>%arrange(RMSE_combined)

yt_combined<- yt_loo%>%select(ModelID, species, AIC, RMSE_loo, var1, var2, var3, var4)%>%
  left_join(yt_LFO5%>%select(ModelID, species, RMSE, var1, var2, var3, var4)%>%
              rename(RMSE_5=RMSE))%>%
  left_join(yt_LFO10%>%select(ModelID, species, RMSE, var1, var2, var3, var4)%>%
              rename(RMSE_10=RMSE))%>%
  mutate(RMSE_combined=RMSE_loo*loo_weight+RMSE_5)%>%arrange(RMSE_combined)

ps_combined<- ps_loo%>%select(ModelID, species, AIC, RMSE_loo, var1, var2, var3, var4)%>%
  left_join(ps_LFO5%>%select(ModelID, species, RMSE, var1, var2, var3, var4)%>%
              rename(RMSE_5=RMSE))%>%
  left_join(ps_LFO10%>%select(ModelID, species, RMSE, var1, var2, var3, var4)%>%
              rename(RMSE_10=RMSE))%>%
  mutate(RMSE_combined=RMSE_loo*loo_weight+RMSE_5)%>%arrange(RMSE_combined)

sb_combined<- sb_loo%>%select(ModelID, species, AIC, RMSE_loo, var1, var2, var3, var4)%>%
  left_join(sb_LFO5%>%select(ModelID, species, RMSE, var1, var2, var3, var4)%>%
              rename(RMSE_5=RMSE))%>%
  left_join(sb_LFO10%>%select(ModelID, species, RMSE, var1, var2, var3, var4)%>%
              rename(RMSE_10=RMSE))%>%
  mutate(RMSE_combined=RMSE_loo*loo_weight+RMSE_5)%>%arrange(RMSE_combined)

#### Generating combined index predictions ####
hk_combined_pred<-hk_loo_pred%>%rename(pred_loo=pred)%>%
  left_join(hk_LFO10_pred%>%rename(pred_10=pred))%>%
  left_join(hk_LFO5_pred%>%rename(pred_5=pred))

sb_combined_pred<-sb_loo_pred%>%rename(pred_loo=pred)%>%
  left_join(sb_LFO10_pred%>%rename(pred_10=pred))%>%
  left_join(sb_LFO5_pred%>%rename(pred_5=pred))

ps_combined_pred<-ps_loo_pred%>%rename(pred_loo=pred)%>%
  left_join(ps_LFO10_pred%>%rename(pred_10=pred))%>%
  left_join(ps_LFO5_pred%>%rename(pred_5=pred))

yt_combined_pred<-yt_loo_pred%>%rename(pred_loo=pred)%>%
  left_join(yt_LFO10_pred%>%rename(pred_10=pred))%>%
  left_join(yt_LFO5_pred%>%rename(pred_5=pred))
#### Time SeriesRMSE#### Time Series ####
##### Yellowtail #####
yt_baseline <-  null_RMSE(yt_dat)

best_loo<-yt_loo%>%arrange(RMSE_loo)%>%select(ModelID,RMSE_loo)%>%left_join(yt_loo_pred)%>%
  rename(RMSE=RMSE_loo)%>%
  mutate(Criteria="LOO")%>%
  filter(RMSE ==min(RMSE))
  
best_lfo5<- yt_LFO5%>%arrange(RMSE)%>%select(ModelID,RMSE)%>%left_join(yt_LFO5_pred)%>%
  mutate(Criteria="LFO5")%>%
  filter(RMSE ==min(RMSE))

best_lfo10<-yt_LFO10%>%arrange(RMSE)%>%select(ModelID,RMSE)%>%left_join(yt_LFO10_pred)%>%
  mutate(Criteria="LFO10")%>%
  filter(RMSE ==min(RMSE))

best_combined<-yt_combined%>%arrange(RMSE_combined)%>%select(ModelID,RMSE_combined)%>%left_join(yt_combined_pred)%>%
  mutate(Criteria="Combined")%>%
  filter(RMSE_combined ==min(RMSE_combined))

yt_lfo10_plot<-ggplot()+
  geom_line(data=yt_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=yt_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=yt_LFO10_pred%>%
               filter(ModelID %in% unique(best_loo$ModelID)), (aes(y=pred, x=year)), col=col_dv[1], pch=17, cex=2)+
  geom_point(data=yt_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo10$ModelID)), (aes(y=pred, x=year)), col=col_dv[5], pch=16, cex=2)+
  geom_point(data=yt_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo5$ModelID)), (aes(y=pred, x=year)), col=col_dv[3], pch=18, cex=2)+
  geom_point(data=best_combined%>%
               filter(ModelID %in% unique(best_combined$ModelID)), (aes(y=pred_10, x=year)), col=col_dv[7], pch=15, cex=2)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))

yt_loo_plot<-ggplot()+
  geom_point(data=yt_dat,aes( y=Y_rec, x=year), cex=2)+
  geom_ribbon(data=yt_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_line(data=yt_dat,aes( y=Y_rec, x=year), cex=0.25)+
  geom_line(data=yt_loo_pred%>%
               filter(ModelID %in% unique(best_loo$ModelID[1])), (aes(y=pred, x=year)), col=col_dv[1], cex=1.5)+
  geom_line(data=yt_loo_pred%>%
               filter(ModelID %in% unique(best_lfo10$ModelID[1])), (aes(y=pred, x=year)), col=col_dv[5], cex=1.5)+
  geom_line(data=yt_loo_pred%>%
               filter(ModelID %in% unique(best_lfo5$ModelID[1])), (aes(y=pred, x=year)), col=col_dv[3],  cex=1.5)+
  geom_line(data=yt_combined_pred%>%
              filter(ModelID %in% unique(best_combined$ModelID[1])), (aes(y=pred_loo, x=year)), col=col_dv[7],  cex=1.5)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))  

yt_ts<- annotate_figure(top="Yellowtail", bottom="Year",
                ggarrange(yt_lfo10_plot,yt_loo_plot,z, ncol = 3, nrow = 1,widths=c(1,1,0.6), labels=c("E.", "F.","")))


pred_yt_ensemble<- yt_LFO10_pred%>%
  filter(ModelID %in% unique(best_loo$ModelID)|ModelID %in% unique(best_lfo10$ModelID)|ModelID %in% unique(best_loo$ModelID))%>%
  group_by(year)%>%
  summarise(pred_ensemble=mean(pred))

yt_ensemble_plot<-ggplot()+
  geom_line(data=yt_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=yt_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=pred_yt_ensemble, aes(y=pred_ensemble, x=year), col=col_dv[1], pch=17, cex=2)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  ggtitle("Yellowtail")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))
##### Sablefish #####
sb_baseline <-  null_RMSE(sb_dat)

best_loo<-sb_loo%>%arrange(RMSE_loo)%>%select(ModelID,RMSE_loo)%>%left_join(sb_loo_pred)%>%
  rename(RMSE=RMSE_loo)%>%
  mutate(Criteria="LOO")%>%
  filter(RMSE ==min(RMSE))

best_lfo5<- sb_LFO5%>%arrange(RMSE)%>%select(ModelID,RMSE)%>%left_join(sb_LFO5_pred)%>%
  mutate(Criteria="LFO5")%>%
  filter(RMSE ==min(RMSE))

best_lfo10<-sb_LFO10%>%arrange(RMSE)%>%select(ModelID,RMSE)%>%left_join(sb_LFO10_pred)%>%
  mutate(Criteria="LFO10")%>%
  filter(RMSE ==min(RMSE))

best_combined<-sb_combined%>%arrange(RMSE_combined)%>%select(ModelID,RMSE_combined)%>%left_join(sb_combined_pred)%>%
  mutate(Criteria="Combined")%>%
  filter(RMSE_combined ==min(RMSE_combined))

sb_lfo10_plot<-ggplot()+
  geom_line(data=sb_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=sb_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=sb_LFO10_pred%>%
               filter(ModelID %in% unique(best_loo$ModelID)), (aes(y=pred, x=year)), col=col_dv[1], pch=17, cex=2)+
  geom_point(data=sb_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo10$ModelID)), (aes(y=pred, x=year)), col=col_dv[5], pch=16, cex=2)+
  geom_point(data=sb_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo5$ModelID)), (aes(y=pred, x=year)), col=col_dv[3], pch=18, cex=2)+
  geom_point(data=best_combined%>%
               filter(ModelID %in% unique(best_combined$ModelID)), (aes(y=pred_10, x=year)), col=col_dv[7], pch=15, cex=2)+

  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))

sb_loo_plot<-ggplot()+
  geom_point(data=sb_dat,aes( y=Y_rec, x=year), cex=2)+
  geom_ribbon(data=sb_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_line(data=sb_dat,aes( y=Y_rec, x=year,col="Recruitment \n Deviation"), cex=0.25)+
  geom_line(data=sb_loo_pred%>%
              filter(ModelID %in% unique(best_loo$ModelID[1])), aes(y=pred, x=year, col="LOO"),  cex=1.5)+
  geom_line(data=sb_loo_pred%>%
              filter(ModelID %in% unique(best_lfo10$ModelID[1])), aes(y=pred, x=year, col="LFO10"), cex=1.5)+
  geom_line(data=sb_loo_pred%>%
              filter(ModelID %in% unique(best_lfo5$ModelID[1])), aes(y=pred, x=year, col="LFO5"),  cex=1.5)+
  geom_line(data=sb_combined_pred%>%
            filter(ModelID %in% unique(best_combined$ModelID[1])), (aes(y=pred_loo, x=year, col="Combined")),   cex=1.5)+
  
   geom_hline(yintercept = 0, lty=2)+
  scale_color_manual(name = "Selection Criteria",
                    values = c("Recruitment \n Deviation"= "black",
                               "LOO" = col_dv[1],
                               "LFO5" = col_dv[3],
                               "LFO10" = col_dv[5],
                               "Combined"=col_dv[7]
                               
                               
                    ), 
                    labels = c("Combined","Leave-future-out 10 year","Leave-future-out 5 year","Leave-one-out","Recruitment \n Deviations"))+
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))  
sb_loo_plot
sb_ts<- annotate_figure(top="Sablefish", 
                        ggarrange(sb_lfo10_plot,sb_loo_plot, ncol = 2, nrow = 1, labels=c("C.", "D."), widths=c(1,1.75)))
sb_ts

pred_sb_ensemble<- sb_LFO10_pred%>%
  filter(ModelID %in% unique(best_loo$ModelID)|ModelID %in% unique(best_lfo10$ModelID)|ModelID %in% unique(best_loo$ModelID))%>%
  group_by(year)%>%
  summarise(pred_ensemble=mean(pred))

sb_ensemble_plot<-ggplot()+
  geom_line(data=sb_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=sb_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=pred_sb_ensemble, aes(y=pred_ensemble, x=year), col=col_dv[1], pch=17, cex=2)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  ggtitle("Sablefish")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))
##### Petrale #####
ps_baseline <-  null_RMSE(ps_dat)

best_loo<-ps_loo%>%arrange(RMSE_loo)%>%select(ModelID,RMSE_loo)%>%left_join(ps_loo_pred)%>%
  rename(RMSE=RMSE_loo)%>%
  mutate(Criteria="LOO")%>%
  filter(RMSE ==min(RMSE))

best_lfo5<- ps_LFO5%>%arrange(RMSE)%>%select(ModelID,RMSE)%>%left_join(ps_LFO5_pred)%>%
  mutate(Criteria="LFO5")%>%
  filter(RMSE ==min(RMSE))

best_lfo10<-ps_LFO10%>%arrange(RMSE)%>%select(ModelID,RMSE)%>%left_join(ps_LFO10_pred)%>%
  mutate(Criteria="LFO10")%>%
  filter(RMSE ==min(RMSE))

best_combined<-yt_combined%>%arrange(RMSE_combined)%>%select(ModelID,RMSE_combined)%>%left_join(yt_combined_pred)%>%
  mutate(Criteria="Combined")%>%
  filter(RMSE_combined ==min(RMSE_combined))

ps_lfo10_plot<-ggplot()+
  geom_line(data=ps_dat,aes( y=Y_rec, x=year))+
  geom_ribbon(data=ps_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=ps_LFO10_pred%>%
               filter(ModelID %in% unique(best_loo$ModelID)), (aes(y=pred, x=year)), col=col_dv[1], pch=17, cex=2)+
  geom_point(data=ps_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo10$ModelID)), (aes(y=pred, x=year)), col=col_dv[5], pch=16, cex=2)+
  geom_point(data=ps_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo5$ModelID)), (aes(y=pred, x=year)), col=col_dv[3], pch=18, cex=2)+
  geom_point(data=best_combined%>%
               filter(ModelID %in% unique(best_combined$ModelID)), (aes(y=pred_10, x=year)), col=col_dv[7], pch=15, cex=2)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  ggtitle("One-year-ahead predictions")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))

ps_loo_plot<-ggplot()+
  geom_point(data=ps_dat,aes( y=Y_rec, x=year), cex=2)+
  geom_ribbon(data=ps_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_line(data=ps_dat,aes( y=Y_rec, x=year), cex=0.25)+
  geom_line(data=ps_loo_pred%>%
              filter(ModelID %in% unique(best_loo$ModelID[1])), (aes(y=pred, x=year)), col=col_dv[1], cex=1.5)+
  geom_line(data=ps_loo_pred%>%
              filter(ModelID %in% unique(best_lfo10$ModelID[1])), (aes(y=pred, x=year)), col=col_dv[5], cex=1.5)+
  geom_line(data=ps_loo_pred%>%
              filter(ModelID %in% unique(best_lfo5$ModelID[1])), (aes(y=pred, x=year)), col=col_dv[3],  cex=1.5)+
  geom_line(data=ps_combined_pred%>%
              filter(ModelID %in% unique(best_combined$ModelID[1])), (aes(y=pred_loo, x=year)), col=col_dv[7],  cex=1.5)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  ggtitle("Leave-one-out predictions")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))  

ps_ts<- annotate_figure(top="Petrale Sole", 
                        ggarrange(ps_lfo10_plot,ps_loo_plot,z, ncol = 3, nrow = 1, widths=c(1,1,0.6),labels=c("A.", "B.","")))



pred_ps_ensemble<- ps_LFO10_pred%>%
  filter(ModelID %in% unique(best_loo$ModelID)|ModelID %in% unique(best_lfo10$ModelID)|ModelID %in% unique(best_loo$ModelID))%>%
  group_by(year)%>%
  summarise(pred_ensemble=mean(pred))

ps_ensemble_plot<-ggplot()+
  geom_line(data=ps_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=ps_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=pred_ps_ensemble, aes(y=pred_ensemble, x=year), col=col_dv[1], pch=17, cex=2)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  ggtitle("Petrale Sole")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))

##### Hake #####
hk_baseline <-  null_RMSE(hk_dat)

best_loo<-hk_loo%>%arrange(RMSE_loo)%>%select(ModelID,RMSE_loo)%>%left_join(hk_loo_pred)%>%
  rename(RMSE=RMSE_loo)%>%
  mutate(Criteria="LOO")%>%
  filter(RMSE ==min(RMSE))

best_lfo5<- hk_LFO5%>%arrange(RMSE)%>%select(ModelID,RMSE)%>%left_join(hk_LFO5_pred)%>%
  mutate(Criteria="LFO5")%>%
  filter(RMSE ==min(RMSE))

best_lfo10<-hk_LFO10%>%arrange(RMSE)%>%select(ModelID,RMSE)%>%left_join(hk_LFO10_pred)%>%
  mutate(Criteria="LFO10")%>%
  filter(RMSE ==min(RMSE))

best_combined<-hk_combined%>%arrange(RMSE_combined)%>%select(ModelID,RMSE_combined)%>%left_join(hk_combined_pred)%>%
  mutate(Criteria="Combined")%>%
  filter(RMSE_combined ==min(RMSE_combined))


hk_lfo10_plot<-ggplot()+
  geom_line(data=hk_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=hk_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=hk_LFO10_pred%>%
               filter(ModelID %in% unique(best_loo$ModelID)), (aes(y=pred, x=year)), col=col_dv[1], pch=17, cex=2)+
  geom_point(data=hk_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo10$ModelID)), (aes(y=pred, x=year)), col=col_dv[5], pch=16, cex=2)+
  geom_point(data=hk_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo5$ModelID)), (aes(y=pred, x=year)), col=col_dv[3], pch=18, cex=2)+
  geom_point(data=best_combined%>%
               filter(ModelID %in% unique(best_combined$ModelID)), (aes(y=pred_10, x=year)), col=col_dv[7], pch=15, cex=2)+
  
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))

hk_loo_plot<-ggplot()+
  geom_point(data=hk_dat,aes( y=Y_rec, x=year), cex=2)+
  geom_ribbon(data=hk_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_line(data=hk_dat,aes( y=Y_rec, x=year), cex=0.25)+
  geom_line(data=hk_loo_pred%>%
              filter(ModelID %in% unique(best_loo$ModelID[1])), (aes(y=pred, x=year)), col=col_dv[1], cex=1.5)+
  geom_line(data=hk_loo_pred%>%
              filter(ModelID %in% unique(best_lfo10$ModelID[1])), (aes(y=pred, x=year)), col=col_dv[5], cex=1.5)+
  geom_line(data=hk_loo_pred%>%
              filter(ModelID %in% unique(best_lfo5$ModelID[1])), (aes(y=pred, x=year)), col=col_dv[3],  cex=1.5)+
  geom_line(data=hk_combined_pred%>%
              filter(ModelID %in% unique(best_combined$ModelID[1])), (aes(y=pred_loo, x=year)), col=col_dv[7],  cex=1.5)+
  
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))  

hk_ts<- annotate_figure(top="Hake", bottom="Year",
                        ggarrange(hk_lfo10_plot,hk_loo_plot,z, ncol = 3, nrow = 1,widths=c(1,1,0.6), labels=c("G.", "H.","")))



pred_hk_ensemble<- hk_LFO10_pred%>%
  filter(ModelID %in% unique(best_loo$ModelID)|ModelID %in% unique(best_lfo10$ModelID)|ModelID %in% unique(best_loo$ModelID))%>%
    group_by(year)%>%
    summarise(pred_ensemble=mean(pred))

hk_ensemble_plot<-ggplot()+
  geom_line(data=hk_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=hk_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=pred_hk_ensemble, aes(y=pred_ensemble, x=year), col=col_dv[1], pch=17, cex=2)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  ggtitle("Hake")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))

#### Cobining plots ####
all_ts<- ggarrange(ps_ts, sb_ts,yt_ts, hk_ts,ncol = 1, nrow = 4)
all_ts
pdf(file = "Output/Figures/Combined/all_ts.pdf", width =9, height =12)
all_ts
dev.off()

png(file = "Output/Figures/Combined/all_ts.png",width = 900, height = 1200, res = 100)
all_ts
dev.off()


all_ensemble<- ggarrange(ps_ensemble_plot, sb_ensemble_plot,yt_ensemble_plot,hk_ensemble_plot,ncol = 2, nrow = 2)
all_ensemble
pdf(file = "Output/Figures/Combined/all_ensemble.pdf", width =9, height =7)
all_ts
dev.off()

png(file = "Output/Figures/Combined/all_ensemble.png",width = 900, height = 700, res = 100)
all_ensemble
dev.off()

#### Time Series Covariates ####
df$value[is.na(df$value)] <- ""

full_table_sb<- sb_loo%>%
 mutate(var2 = replace_na(var2, ""))%>%
  mutate(var3 = replace_na(var3, ""))%>%
  mutate(variables = paste(var1, var2, var3, sep=' '))%>%
  select(ModelID, species, RMSE_loo, dev.ex, AIC, Hit, variables)%>%
  left_join(sb_lfo5%>%
              rename(RMSE5=RMSE)%>%
              select(ModelID, RMSE5))%>%
  left_join(sb_lfo10%>%
    rename(RMSE10=RMSE)%>%
    select(ModelID, RMSE10))%>%
  left_join(sb_combined%>%
              select(ModelID,RMSE_combined))

full_table_ps<- ps_loo%>%
  mutate(var2 = replace_na(var2, ""))%>%
  mutate(var3 = replace_na(var3, ""))%>%
  mutate(variables = paste(var1, var2, var3, sep=' '))%>%
  select(ModelID, species, RMSE_loo, dev.ex, AIC, Hit, variables)%>%
  left_join(ps_lfo5%>%
              rename(RMSE5=RMSE)%>%
              select(ModelID, RMSE5))%>%
  left_join(ps_lfo10%>%
              rename(RMSE10=RMSE)%>%
              select(ModelID, RMSE10))%>%
  left_join(ps_combined%>%
              select(ModelID,RMSE_combined))

full_table_yt<- yt_loo%>%
  mutate(var2 = replace_na(var2, ""))%>%
  mutate(var3 = replace_na(var3, ""))%>%
  mutate(variables = paste(var1, var2, var3, sep=' '))%>%
  select(ModelID, species, RMSE_loo, dev.ex, AIC, Hit, variables)%>%
  left_join(yt_lfo5%>%
              rename(RMSE5=RMSE)%>%
              select(ModelID, RMSE5))%>%
  left_join(yt_lfo10%>%
              rename(RMSE10=RMSE)%>%
              select(ModelID, RMSE10))%>%
  left_join(yt_combined%>%
              select(ModelID,RMSE_combined))

full_table_hk<- hk_loo%>%
  mutate(var2 = replace_na(var2, ""))%>%
  mutate(var3 = replace_na(var3, ""))%>%
  mutate(variables = paste(var1, var2, var3, sep=' '))%>%
  select(ModelID, species, RMSE_loo, dev.ex, AIC, Hit, variables)%>%
  left_join(hk_lfo5%>%
              rename(RMSE5=RMSE)%>%
              select(ModelID, RMSE5))%>%
  left_join(hk_lfo10%>%
              rename(RMSE10=RMSE)%>%
              select(ModelID, RMSE10))%>%
  left_join(hk_combined%>%
              select(ModelID,RMSE_combined))
write.csv(full_table_yt, "Output/Data/Combined/yt_table.csv")
write.csv(full_table_ps, "Output/Data/Combined/ps_table.csv")
write.csv(full_table_sb, "Output/Data/Combined/sb_table.csv")
write.csv(full_table_hk, "Output/Data/Combined/hk_table.csv")


#### Best Combined ####



hk_lfo10_plot<-ggplot()+
  geom_line(data=hk_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=hk_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=hk_LFO10_pred%>%
               filter(ModelID %in% unique(best_loo$ModelID)), (aes(y=pred, x=year)), col=col_dv[1], pch=17, cex=2)+
  geom_point(data=hk_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo10$ModelID)), (aes(y=pred, x=year)), col=col_dv[5], pch=16, cex=2)+
  geom_point(data=hk_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo5$ModelID)), (aes(y=pred, x=year)), col=col_dv[3], pch=18, cex=2)+
  geom_point(data=best_combined%>%
               filter(ModelID %in% unique(best_combined$ModelID)), (aes(y=pred_10, x=year)), col=col_dv[7], pch=15, cex=2)+
  

#### Generatung ustable file ####
unstable <-read.csv("Output/Data/Analysis Part 1/ustablevariables_1.csv")%>%
  filter(nspecies>1)

total<-unstable%>%
  count(Type)

total_type<-unstable%>%
  group_by(Type)%>%
  summarise(mean=mean(TotalType))%>%
  left_join(total)%>%
  mutate(rate=n/mean)

count <-ggplot(total_type, aes(y = n, x=Type)) +
  geom_bar(stat = "identity")+
  theme_classic() +
  scale_x_discrete(labels = label_wrap(10))
rate<-ggplot(total_type, aes(y = rate, x=Type)) +
  geom_bar(stat = "identity")+
  theme_classic() +
  scale_x_discrete(labels = label_wrap(10))

countrate<- ggarrange(count, rate, ncol = 2)

pdf(file = "Output/Figures/Analysis Part 1/countrate.pdf", width =8, height = 2.5)
countrate
dev.off()

png(file = "Output/Figures/Analysis Part 1/countrate.png",width = 800, height = 250, res = 100)
countrate
dev.off()
