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

col_dv<-natparks.pals("DeathValley", 7)
names(NatParksPalettes)
#### Read In Data ####
ps <- readRDS("Output/Data/ShortenedTS/ps_model_fits.rds")
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

sb <- readRDS("Output/Data/ShortenedTS/sb_model_fits.rds")
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

yt <- readRDS("Output/Data/ShortenedTS/yt_model_fits.rds")
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
hk <- readRDS("Output/Data/ShortenedTS/hk_model_fits.rds")
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

#### Time Series ####
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


yt_lfo10_plot<-ggplot()+
  geom_line(data=yt_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=yt_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=yt_LFO10_pred%>%
               filter(ModelID %in% unique(best_loo$ModelID)), (aes(y=pred, x=year)), col=col_dv[1], pch=17, cex=2)+
  geom_point(data=yt_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo10$ModelID)), (aes(y=pred, x=year)), col=col_dv[5], pch=16, cex=2)+
  geom_point(data=yt_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo5$ModelID)), (aes(y=pred, x=year)), col=col_dv[3], pch=18, cex=2)+
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
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))  

yt_ts<- annotate_figure(top="Yellowtail", bottom="Year",
                ggarrange(yt_lfo10_plot,yt_loo_plot,z, ncol = 3, nrow = 1,widths=c(1,1,0.6), labels=c("E.", "F.","")))

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



sb_lfo10_plot<-ggplot()+
  geom_line(data=sb_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=sb_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=sb_LFO10_pred%>%
               filter(ModelID %in% unique(best_loo$ModelID)), (aes(y=pred, x=year)), col=col_dv[1], pch=17, cex=2)+
  geom_point(data=sb_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo10$ModelID)), (aes(y=pred, x=year)), col=col_dv[5], pch=16, cex=2)+
  geom_point(data=sb_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo5$ModelID)), (aes(y=pred, x=year)), col=col_dv[3], pch=18, cex=2)+
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
  geom_hline(yintercept = 0, lty=2)+
  scale_color_manual(name = "Selection Criteria",
                    values = c("Recruitment \n Deviation"= "black",
                               "LOO" = col_dv[1],
                               "LFO5" = col_dv[3],
                               "LFO10" = col_dv[5]
                    ), 
                    labels = c("Leave-future-out 10 year","Leave-future-out 5 year","Leave-one-out","Recruitment \n Deviations"))+
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))  
sb_loo_plot
sb_ts<- annotate_figure(top="Sablefish", 
                        ggarrange(sb_lfo10_plot,sb_loo_plot, ncol = 2, nrow = 1, labels=c("C.", "D."), widths=c(1,1.75)))
sb_ts
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


ps_lfo10_plot<-ggplot()+
  geom_line(data=ps_dat,aes( y=Y_rec, x=year))+
  geom_ribbon(data=ps_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=ps_LFO10_pred%>%
               filter(ModelID %in% unique(best_loo$ModelID)), (aes(y=pred, x=year)), col=col_dv[1], pch=17, cex=2)+
  geom_point(data=ps_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo10$ModelID)), (aes(y=pred, x=year)), col=col_dv[5], pch=16, cex=2)+
  geom_point(data=ps_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo5$ModelID)), (aes(y=pred, x=year)), col=col_dv[3], pch=18, cex=2)+
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
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  ggtitle("Leave-one-out predictions")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))  

ps_ts<- annotate_figure(top="Petrale Sole", 
                        ggarrange(ps_lfo10_plot,ps_loo_plot,z, ncol = 3, nrow = 1, widths=c(1,1,0.6),labels=c("A.", "B.","")))


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


hk_lfo10_plot<-ggplot()+
  geom_line(data=hk_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=hk_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=hk_LFO10_pred%>%
               filter(ModelID %in% unique(best_loo$ModelID)), (aes(y=pred, x=year)), col=col_dv[1], pch=17, cex=2)+
  geom_point(data=hk_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo10$ModelID)), (aes(y=pred, x=year)), col=col_dv[5], pch=16, cex=2)+
  geom_point(data=hk_LFO10_pred%>%
               filter(ModelID %in% unique(best_lfo5$ModelID)), (aes(y=pred, x=year)), col=col_dv[3], pch=18, cex=2)+
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
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))  

hk_ts<- annotate_figure(top="Hake", bottom="Year",
                        ggarrange(hk_lfo10_plot,hk_loo_plot,z, ncol = 3, nrow = 1,widths=c(1,1,0.6), labels=c("G.", "H.","")))


all_ts<- ggarrange(ps_ts, sb_ts,yt_ts, hk_ts,ncol = 1, nrow = 4)
all_ts
pdf(file = "Output/Figures/ShortenedTS/all_ts.pdf", width =9, height =12)
all_ts
dev.off()

png(file = "Output/Figures/ShortenedTS/all_ts.png",width = 900, height = 1200, res = 100)
all_ts
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
    select(ModelID, RMSE10))

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
              select(ModelID, RMSE10))

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
              select(ModelID, RMSE10))

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
              select(ModelID, RMSE10))
write.csv(full_table_yt, "Output/Data/ShortenedTS/yt_table.csv")
write.csv(full_table_ps, "Output/Data/ShortenedTS/ps_table.csv")
write.csv(full_table_sb, "Output/Data/ShortenedTS/sb_table.csv")
write.csv(full_table_hk, "Output/Data/ShortenedTS/hk_table.csv")

