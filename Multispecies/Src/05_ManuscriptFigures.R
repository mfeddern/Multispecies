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
library(tidytext)
library(ggpubr)

#### Figure 1 #####
all_marginals <- readRDS( "Output/Data/JuneUpdate/marginals.rds")
all_marginals$cov<-as.character(all_marginals$cov)
str(all_marginals)
all_marginals <-all_marginals%>%
  mutate(cov = fct_recode(cov, "ONIlarv"="oni_larv"))%>%
  mutate(cov = fct_recode(cov, "ONIpre"="oni_pre"))%>%
  mutate(cov = fct_recode(cov, "ONIpjuv"="oni_pjuv"))%>%
  mutate(cov = fct_recode(cov, "PDOpjuv"="pdo_pjuv"))%>%
  mutate(cov = fct_recode(cov, "PDOlarv"="pdo_larv"))
  
cols<- c('#dd4124',"#edd746",'#7cae00','#0f85a0')
xlim<-c(-.05,0.06)

yt_marg<-all_marginals%>%filter(species=="Yellowtail")%>%filter(RMSE=="LOO"|RMSE=="LFO 10")%>% 
  mutate(category = case_when(
      str_detect(cov, "CutiSTIpjuv")  ~ "Yes",
      str_detect(cov, "ONIpjuv") ~ "Yes",
      str_detect(cov, "LSTpjuv") ~ "Yes",
      TRUE                              ~ "No"       # Default catch-all
    ))

relevel_yt_loo<-yt_marg%>%
  filter(RMSE=="LOO")%>%
  arrange(total_rmse_st)
yt_marg$cov <- factor(yt_marg$cov, levels =relevel_yt_loo$cov)


marginal_yt <- ggplot(yt_marg, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y =cov,
  x = total_rmse_st,
  fill = RMSE
)) +
  xlim(c(xlim))+
  ggtitle("Yellowtail Rockfish")+
 # facet_grid(species ~ RMSE, scales = "free_y") + # Use free_y scale
  geom_bar(position="dodge",stat = "identity") +
  labs(x = "", y = "") +
  #scale_fill_gradient(low = "gray100", high = "Darkgreen",limits=c(-.1,0.1)) +
  theme_bw()+
  geom_text(
    aes(label = ifelse(category == "Yes", "*", "")), 
    nudge_y = -.4,x=0.035,  # Nudges the star slightly above the top of the bar
    size = 5         # Adjusts the size of the star
  ) +
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
marginal_yt


hk_marg<-all_marginals%>%filter(species=="Hake")%>%filter(RMSE=="LOO"|RMSE=="LFO 10")
relevel_hk_loo<-hk_marg%>%
  filter(RMSE=="LOO")%>%
  arrange(total_rmse_st)
hk_marg<-hk_marg%>% 
  mutate(category = case_when(
    str_detect(cov, "LSTyolk")  ~ "Yes",
    TRUE                              ~ "No"       # Default catch-all
  ))
hk_marg$cov <- factor(hk_marg$cov, levels =relevel_hk_loo$cov)

marginal_hk <- ggplot(hk_marg, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y =cov,
  x = total_rmse_st,
  fill = RMSE
)) +
  ggtitle("Hake")+
  xlim(c(xlim))+
  # facet_grid(species ~ RMSE, scales = "free_y") + # Use free_y scale
  geom_bar(position="dodge",stat = "identity") +
  labs(x = "", y = "") +
  #scale_fill_gradient(low = "gray100", high = "Darkgreen",limits=c(-.1,0.1)) +
  theme_bw()+
  geom_text(
    aes(label = ifelse(category == "Yes", "*", "")), 
    nudge_y = -.45,x=0.015,  # Nudges the star slightly above the top of the bar
    size = 5        # Adjusts the size of the star
  ) +
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
marginal_hk

ps_marg<-all_marginals%>%filter(species=="Petrale Sole")%>%filter(RMSE=="LOO"|RMSE=="LFO 10")
relevel_ps_loo<-ps_marg%>%
  filter(RMSE=="LOO")%>%
  arrange(total_rmse_st)


ps_marg<-ps_marg%>% 
  mutate(category = case_when(
    str_detect(cov, "DDpre")  ~ "Yes",
    str_detect(cov, "MLDegg") ~ "Yes",
    str_detect(cov, "CSTlarv") ~ "Yes",
    str_detect(cov, "CSTbjuv") ~ "Yes",
    str_detect(cov, "DDpjuv") ~ "Yes",
    str_detect(cov, "LSTlarv") ~ "Yes",
    TRUE                              ~ "No"       # Default catch-all
  ))
ps_marg$cov <- factor(ps_marg$cov, levels =relevel_ps_loo$cov)

marginal_ps <- ggplot(ps_marg, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y =cov,
  x = total_rmse_st,
  fill = RMSE
)) +
  ggtitle("Petrale Sole")+
  xlim(c(xlim))+
  # facet_grid(species ~ RMSE, scales = "free_y") + # Use free_y scale
  geom_bar(position="dodge",stat = "identity") +
  labs(x = "", y = "") +
  #scale_fill_gradient(low = "gray100", high = "Darkgreen",limits=c(-.1,0.1)) +
  theme_bw()+
  geom_text(
    aes(label = ifelse(category == "Yes", "*", "")), 
    nudge_y = -.3,x=0.06,  # Nudges the star slightly above the top of the bar
    size = 5        # Adjusts the size of the star
  ) +
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
marginal_ps

sb_marg<-all_marginals%>%filter(species=="Sablefish")%>%filter(RMSE=="LOO"|RMSE=="LFO 10")
relevel_sb_loo<-sb_marg%>%
  filter(RMSE=="LOO")%>%
  arrange(total_rmse_st)

sb_marg<-sb_marg%>% 
  mutate(category = case_when(
    str_detect(cov, "CSTegg")  ~ "Yes",
    str_detect(cov, "DDpre") ~ "Yes",
    str_detect(cov, "DDegg") ~ "Yes",
    str_detect(cov, "LSTyolk") ~ "Yes",
    str_detect(cov, "DDlarv") ~ "Yes",
    TRUE                              ~ "No"       # Default catch-all
  ))
sb_marg$cov <- factor(sb_marg$cov, levels =relevel_sb_loo$cov)

marginal_sb <- ggplot(sb_marg, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y =cov,
  x = total_rmse_st,
  fill = RMSE
)) +
  ggtitle("Sablefish")+
  xlim(c(xlim))+
  # facet_grid(species ~ RMSE, scales = "free_y") + # Use free_y scale
  geom_bar(position="dodge",stat = "identity") +
  labs(x = "", y = "") +
  geom_text(
    aes(label = ifelse(category == "Yes", "*", "")), 
    nudge_y = -.3,x=0.06,  # Nudges the star slightly above the top of the bar
    size = 5         # Adjusts the size of the star
  ) +
  #scale_fill_gradient(low = "gray100", high = "Darkgreen",limits=c(-.1,0.1)) +
  theme_bw()+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
marginal_sb


marginal_rmse<-ggarrange(marginal_hk,marginal_ps,marginal_yt,marginal_sb, ncol = 2, nrow = 2,     
          widths=c(10.,9), common.legend = TRUE,  
          legend = "bottom", labels = c("A", "B","C","D") )

marginal_rmse<-annotate_figure(marginal_rmse,
                bottom = text_grob("Marginal Mean Improvement RMSE"))


ggsave("Output/Figures/Manuscript/Figure1RMSE.pdf", plot = marginal_rmse, width = 8, height = 7)
ggsave("Output/Figures/Manuscript/Figure1RMSE.png", plot = marginal_rmse, width = 8, height = 7,bg = "white", dpi = 300)

#### Figure 2 ####
arrange_all_marginals<-all_marginals%>%
  filter(RMSE=="Combined")%>%
  arrange(total_rmse_st)
relevel_all_marginals$cov_species <- factor(relevel_all_marginals$cov_species, levels =arrange_all_marginals$cov_species)

marginal<- ggplot(relevel_all_marginals, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y =reorder_within(cov, total_rmse_st, species),
  x = total_rmse_st,
  fill = RMSE
)) +
  xlim(c(-0.05,0.07))+
  facet_wrap(~species, scales = "free_y") + # Use free_y scale
  geom_bar(position="dodge",stat = "identity") +
  labs(x = "", y = "") +
  scale_y_reordered() + 
  labs(x = "Mean Marginal Mean Improvement", y = "")+
  #scale_fill_gradient(low = "gray100", high = "Darkgreen",limits=c(-.1,0.1)) +
  theme_classic()+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
marginal

ggsave("Output/Figures/Manuscript/Figure2Combined.pdf", plot = marginal, width = 6, height = 6)
ggsave("Output/Figures/Manuscript/Figure2Combined.png", plot = marginal, width = 6, height = 6,bg = "white", dpi = 300)



#### Figure 4 ####

sb_loo_pred<-data.frame(sb[["LOO"]][["predicted"]])
sb_LFO5_pred<-data.frame(sb[["LFO5"]][["predicted"]])
sb_LFO10_pred<-data.frame(sb[["LFO10"]][["predicted"]])

ps_loo_pred<-data.frame(ps[["LOO"]][["predicted"]])
ps_LFO5_pred<-data.frame(ps[["LFO5"]][["predicted"]])
ps_LFO10_pred<-data.frame(ps[["LFO10"]][["predicted"]])

yt_loo_pred<-data.frame(yt[["LOO"]][["predicted"]])
yt_LFO5_pred<-data.frame(yt[["LFO5"]][["predicted"]])
yt_LFO10_pred<-data.frame(yt[["LFO10"]][["predicted"]])

hk_loo_pred<-data.frame(hk[["LOO"]][["predicted"]])
hk_LFO5_pred<-data.frame(hk[["LFO5"]][["predicted"]])
hk_LFO10_pred<-data.frame(hk[["LFO10"]][["predicted"]])

modselection<-readRDS("Output/Data/JuneUpdate/model_selection.rds")

ytselect<-modselection%>%filter(species=="Yellowtail" &Rank=="Best"&selection!="AIC"&selection!="LFO-5")
pred_yt_ensemble<- yt_LFO10_pred%>%
  filter(ModelID %in% unique(ytselect)$ModelID)%>%
  left_join(ytselect%>%select(ModelID,selection))
  #group_by(year)%>%
  #summarise(pred_ensemble=mean(pred))

yt_10pred<-ggplot()+
  geom_line(data=yt_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=yt_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=pred_yt_ensemble, aes(y=pred, x=year, col=selection),  pch=17, cex=2)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
  scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Yellowtail Rockfish")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
theme(legend.position = "none",plot.title = element_text(hjust = 0.5))



hkselect<-modselection%>%filter(species=="Hake" & Rank=="Best"&selection!="AIC"&selection!="LFO-5")
pred_hk_ensemble<- hk_LFO10_pred%>%
  filter(ModelID %in% unique(hkselect)$ModelID)%>%
  left_join(hkselect%>%select(ModelID,selection))
#group_by(year)%>%
#summarise(pred_ensemble=mean(pred))

hk_10pred<-ggplot()+
  geom_line(data=hk_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=hk_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=pred_hk_ensemble, aes(y=pred, x=year, col=selection),  pch=17, cex=2)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
  scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Hake")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
theme(legend.position = "none",plot.title = element_text(hjust = 0.5))


psselect<-modselection%>%filter(species=="Petrale Sole" & Rank=="Best"&selection!="AIC"&selection!="LFO-5")
pred_ps_ensemble<- ps_LFO10_pred%>%
  filter(ModelID %in% unique(psselect)$ModelID)%>%
  left_join(psselect%>%select(ModelID,selection))
#group_by(year)%>%
#summarise(pred_ensemble=mean(pred))

ps_10pred<-ggplot()+
  geom_line(data=ps_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=ps_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=pred_ps_ensemble, aes(y=pred, x=year, col=selection),  pch=17, cex=2,position = position_jitter(width = 0.125, height = 0))+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
  scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Petrale Sole")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
theme(legend.position = "none",plot.title = element_text(hjust = 0.5))



sbselect<-modselection%>%filter(species=="Sablefish" & Rank=="Best"&selection!="AIC"&selection!="LFO-5")
pred_sb_ensemble<- sb_LFO10_pred%>%
  filter(ModelID %in% unique(sbselect)$ModelID)%>%
  filter(ModelID !=8)%>%
  left_join(sbselect%>%select(ModelID,selection))
#group_by(year)%>%
#summarise(pred_ensemble=mean(pred))

sb_10pred<-ggplot()+
  geom_line(data=sb_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=sb_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=pred_sb_ensemble, aes(y=pred, x=year, col=selection),  pch=17, cex=2,position = position_jitter(width = 0.125, height = 0))+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
  scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Sablefish")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
theme(legend.position = "none",plot.title = element_text(hjust = 0.5))


oneyear<-ggarrange(hk_10pred,ps_10pred,yt_10pred,sb_10pred, ncol = 2, nrow = 2,     
                         widths=c(10.,9), common.legend = TRUE,  
                         legend = "right", labels = c("A", "B","C","D") )

oneyear<-annotate_figure(oneyear,left = text_grob("ln(recruitment deviations)",rot = 90),
                               bottom = text_grob("Year"))

ggsave("Output/Figures/Manuscript/Figure4oneyear.pdf", plot = oneyear, width = 7, height = 5.5)
ggsave("Output/Figures/Manuscript/Figure4oneyear.png", plot = oneyear, width = 7, height = 5.5,bg = "white", dpi = 300)



##### Figure 5 #####

ytselect<-modselection%>%filter(species=="Yellowtail" &Rank=="Best"&selection!="AIC"&selection!="LFO-5")
pred_yt_ensemble<- yt_loo_pred%>%
  filter(ModelID %in% unique(ytselect)$ModelID)%>%
  left_join(ytselect%>%select(ModelID,selection))
#group_by(year)%>%
#summarise(pred_ensemble=mean(pred))

yt_pred<-ggplot()+
  geom_line(data=yt_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=yt_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=pred_yt_ensemble, aes(y=pred, x=year, col=selection),  pch=17, cex=2)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
 # scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Yellowtail Rockfish")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))



hkselect<-modselection%>%filter(species=="Hake" & Rank=="Best"&selection!="AIC"&selection!="LFO-5")
pred_hk_ensemble<- hk_loo_pred%>%
  filter(ModelID %in% unique(hkselect)$ModelID)%>%
  left_join(hkselect%>%select(ModelID,selection))
#group_by(year)%>%
#summarise(pred_ensemble=mean(pred))

hk_pred<-ggplot()+
  geom_line(data=hk_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=hk_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=pred_hk_ensemble, aes(y=pred, x=year, col=selection),  pch=17, cex=2)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
  #scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Hake")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))


psselect<-modselection%>%filter(species=="Petrale Sole" & Rank=="Best"&selection!="AIC"&selection!="LFO-5")
pred_ps_ensemble<- ps_loo_pred%>%
  filter(ModelID %in% unique(psselect)$ModelID)%>%
  left_join(psselect%>%select(ModelID,selection))
#group_by(year)%>%
#summarise(pred_ensemble=mean(pred))

ps_pred<-ggplot()+
  geom_line(data=ps_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=ps_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=pred_ps_ensemble, aes(y=pred, x=year, col=selection),  pch=17, cex=2,position = position_jitter(width = 0.125, height = 0))+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
  #scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Petrale Sole")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))



sbselect<-modselection%>%filter(species=="Sablefish" & Rank=="Best"&selection!="AIC"&selection!="LFO-5")
pred_sb_ensemble<- sb_loo_pred%>%
  filter(ModelID %in% unique(sbselect)$ModelID)%>%
  filter(ModelID !=8)%>%
  left_join(sbselect%>%select(ModelID,selection))
#group_by(year)%>%
#summarise(pred_ensemble=mean(pred))

sb_pred<-ggplot()+
  geom_line(data=sb_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=sb_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=pred_sb_ensemble, aes(y=pred, x=year, col=selection),  pch=17, cex=2,position = position_jitter(width = 0.125, height = 0))+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
  #scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Sablefish")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))


outofsample<-ggarrange(hk_pred,ps_pred,yt_pred,sb_pred, ncol = 2, nrow = 2,     
                   widths=c(10.,9), common.legend = TRUE,  
                   legend = "right", labels = c("A", "B","C","D") )

outofsample<-annotate_figure(outofsample,left = text_grob("ln(recruitment deviations)",rot = 90),
                         bottom = text_grob("Year"))

ggsave("Output/Figures/Manuscript/Figure5outofsample.pdf", plot = outofsample, width = 8, height = 5.5)
ggsave("Output/Figures/Manuscript/Figure5outofsample.png", plot = outofsample, width = 8, height = 5.5,bg = "white", dpi = 300)

#### Supplemental Figures ####
##### previous models #####
ps1 <- readRDS("Output/Data/Manuscript/ps_model_best_og1.rds")
ps2 <- readRDS("Output/Data/Manuscript/ps_model_best_og2.rds")
sb <- readRDS("Output/Data/Manuscript/sb_model_best_og1.rds")
yt <- readRDS("Output/Data/Manuscript/yt_model_best_og1.rds")


sb_loo_pred<-data.frame(sb[["LOO"]][["predicted"]])
sb_LFO5_pred<-data.frame(sb[["LFO5"]][["predicted"]])
sb_LFO10_pred<-data.frame(sb[["LFO10"]][["predicted"]])

ps1_loo_pred<-data.frame(ps1[["LOO"]][["predicted"]])
ps1_LFO5_pred<-data.frame(ps1[["LFO5"]][["predicted"]])
ps1_LFO10_pred<-data.frame(ps1[["LFO10"]][["predicted"]])

ps2_loo_pred<-data.frame(ps2[["LOO"]][["predicted"]])
ps2_LFO5_pred<-data.frame(ps2[["LFO5"]][["predicted"]])
ps2_LFO10_pred<-data.frame(ps2[["LFO10"]][["predicted"]])


yt_loo_pred<-data.frame(yt[["LOO"]][["predicted"]])
yt_LFO5_pred<-data.frame(yt[["LFO5"]][["predicted"]])
yt_LFO10_pred<-data.frame(yt[["LFO10"]][["predicted"]])


sb_10pred<-ggplot()+
  geom_line(data=sb_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=sb_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=sb_LFO10_pred, aes(y=pred, x=year),  pch=17, cex=2,position = position_jitter(width = 0.125, height = 0))+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
  scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Sablefish")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
sb_10pred

ps1_10pred<-ggplot()+
  geom_line(data=ps_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=ps_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=ps1_LFO10_pred, aes(y=pred, x=year),  pch=17, cex=2,position = position_jitter(width = 0.125, height = 0))+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
  scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Petrale Sole 1")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
ps1_10pred

ps2_10pred<-ggplot()+
  geom_line(data=ps_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=ps_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=ps2_LFO10_pred, aes(y=pred, x=year),  pch=17, cex=2,position = position_jitter(width = 0.125, height = 0))+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
  scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Petrale Sole 2")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
ps2_10pred

yt_10pred<-ggplot()+
  geom_line(data=yt_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=yt_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=yt_LFO10_pred, aes(y=pred, x=year),  pch=17, cex=2,position = position_jitter(width = 0.125, height = 0))+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
  scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Yellowtail Rockfish")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
yt_10pred

oneyear<-ggarrange(ps1_10pred,ps2_10pred,yt_10pred,sb_10pred, ncol = 2, nrow = 2,     
                   widths=c(10.,9), common.legend = TRUE,  
                   legend = "right", labels = c("A", "B","C","D") )

oneyear<-annotate_figure(oneyear,left = text_grob("ln(recruitment deviations)",rot = 90),
                         bottom = text_grob("Year"))

ggsave("Output/Figures/Manuscript/FigureS3oneyear.pdf", plot = oneyear, width = 7, height = 5.5)
ggsave("Output/Figures/Manuscript/FigureS3oneyear.png", plot = oneyear, width = 7, height = 5.5,bg = "white", dpi = 300)




yt_pred<-ggplot()+
  geom_line(data=yt_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=yt_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=yt_loo_pred, aes(y=pred, x=year),  pch=17, cex=2)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
  # scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Yellowtail Rockfish")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))


sb_pred<-ggplot()+
  geom_line(data=sb_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=sb_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=sb_loo_pred, aes(y=pred, x=year),  pch=17, cex=2)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
  # scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Sablefish")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))


ps1_pred<-ggplot()+
  geom_line(data=ps_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=ps_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=ps1_loo_pred, aes(y=pred, x=year),  pch=17, cex=2)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
  # scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Petrale Sole 1")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))


ps2_pred<-ggplot()+
  geom_line(data=ps_dat,aes( y=Y_rec, x=year), lwd=0.25)+
  geom_ribbon(data=ps_dat,aes(ymin=Y_rec-1.96*sd, ymax=Y_rec+1.96*sd, x=year), alpha=0.1, col='lightgray')+
  geom_point(data=ps2_loo_pred, aes(y=pred, x=year),  pch=17, cex=2)+
  geom_hline(yintercept = 0, lty=2)+
  ylab("")+
  xlab("")+
  #ylim(c(-4,4))+
  # scale_x_continuous( limits = c(2009, 2018),breaks = seq(2008,2018, by = 2)) +
  ggtitle("Petrale Sole 2")+
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))


oneyear<-ggarrange(ps1_pred,ps2_pred,yt_pred,sb_pred, ncol = 2, nrow = 2,     
                   widths=c(10.,9), common.legend = TRUE,  
                   legend = "right", labels = c("A", "B","C","D") )

oneyear<-annotate_figure(oneyear,left = text_grob("ln(recruitment deviations)",rot = 90),
                         bottom = text_grob("Year"))

ggsave("Output/Figures/Manuscript/FigureS2loo.pdf", plot = oneyear, width = 7, height = 5.5)
ggsave("Output/Figures/Manuscript/FigureS2loo.png", plot = oneyear, width = 7, height = 5.5,bg = "white", dpi = 300)
