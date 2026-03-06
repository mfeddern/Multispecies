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

#### Reading in the data ####
### Functions ###
RMSE_improvement <-function(results,baseline_rmse,gam_model, covariates){
  
  # we can calculate the marginal improvement for each covariate
  results$n_cov <- ifelse(!is.na(results$var1), 1, 0) + ifelse(!is.na(results$var2), 1, 0) +
    ifelse(!is.na(results$var3), 1, 0) + ifelse(!is.na(results$var4), 1, 0)
  marginals <- data.frame(cov = covariates, "rmse_01" = NA, "rmse_12" = NA, "rmse_23" = NA,
                          "aic_01" = NA, "aic_12" = NA, "aic_23" = NA)
  results<-results%>%filter(n_cov<4)%>%select(-var4)
  results%>%filter(n_cov==1)
  for(i in 1:length(covariates)) {
    sub <- dplyr::filter(results, n_cov == 1,
                         var1 == covariates[i])
    marginals$rmse_01[i] <- (baseline_rmse-sub$RMSE) / baseline_rmse
    marginals$aic_01[i] <- AIC(gam_model) - sub$AIC
    
    # next look at all values of models that have 2 covariates and include this model
    sub1 <- dplyr::filter(results, n_cov == 1)
    sub2 <- dplyr::filter(results, n_cov == 2) %>%
      dplyr::mutate(keep = ifelse(var1 == covariates[i],1,0) + ifelse(var2 == covariates[i],1,0)) %>%
      dplyr::filter(keep == 1) %>% dplyr::select(-keep)
    # loop over every variable in sub2, and find the simpler model in sub1 that just represents the single covariate
    sub2$rmse_diff <- 0
    sub2$AIC_diff <- 0
    for(j in 1:nrow(sub2)) {
      vars <- sub2[j,c("var1", "var2")]
      vars <- vars[which(vars != covariates[i])]
      indx <- which(sub1$var1 == as.character(vars))
      sub2$rmse_diff[j] <- (sub1$RMSE[indx]-sub2$RMSE[j]) / sub1$RMSE[indx]
      sub2$AIC_diff[j] <- sub1$AIC[indx] - sub2$AIC[j]
    }
    
    # Finally compare models with 3 covariates to models with 2 covariates
    sub2_all <- dplyr::filter(results, n_cov == 2)
    # Apply a function across the rows to sort the values in var1 and var2
    sorted_names <- t(apply(sub2_all[, c("var1", "var2")], 1, function(x) sort(x)))
    # Replace the original columns with the sorted data
    sub2_all$var1 <- sorted_names[, 1]
    sub2_all$var2 <- sorted_names[, 2]
    
    sub3 <- dplyr::filter(results, n_cov == 3) %>%
      dplyr::mutate(keep = ifelse(var1 == covariates[i],1,0) + ifelse(var2 == covariates[i],1,0) + ifelse(var3 == covariates[i],1,0)) %>%
      dplyr::filter(keep == 1) %>% dplyr::select(-keep)
    sub3$rmse_diff <- 0
    sub3$AIC_diff <- 0
    for(j in 1:nrow(sub3)) {
      vars <- sub3[j,c("var1", "var2","var3")]
      vars <- sort(as.character(vars[which(vars != covariates[i])]))
      # find the same in sub2
      indx <- which(paste(sub2_all$var1, sub2_all$var2) == paste(vars, collapse=" "))
      sub3$rmse_diff[j] <- (sub2_all$RMSE[indx]-sub3$RMSE[j]) / sub2_all$RMSE[indx]
      sub3$AIC_diff[j] <- sub2_all$AIC[indx] - sub3$AIC[j]
    }
    
    # Fill in summary stats
    marginals$rmse_12[i] <- mean(sub2$rmse_diff)
    marginals$aic_12[i] <- mean(sub2$AIC_diff)
    marginals$rmse_23[i] <- mean(sub3$rmse_diff)
    marginals$aic_23[i] <- mean(sub3$AIC_diff)
  }
  return(marginals)
}
null_RMSE <-function(dat){
  null<-dat%>%mutate(sr_null = 0)
  sqerror<-(null$sr_null-null$Y_rec)^2
  rmse_sr_full <- sqrt(mean(sqerror, na.rm=T))
  return(rmse_sr_full)
}

#### Read In Data ####

ps <- readRDS("Output/Data/ps_model_fits.rds")
ps_loo<-data.frame(ps[["LOO"]][["results"]])
ps_LFO5<-data.frame(ps[["LFO5"]][["results"]])
ps_LFO10<-data.frame(ps[["LFO10"]][["results"]])

sb <- readRDS("Output/Data/sb_model_fits.rds")
sb_loo<-data.frame(sb[["LOO"]][["results"]])
sb_LFO5<-data.frame(sb[["LFO5"]][["results"]])
sb_LFO10<-data.frame(sb[["LFO10"]][["results"]])

yt <- readRDS("Output/Data/yt_model_fits.rds")

yt_loo<-data.frame(yt[["LOO"]][["results"]])
yt_LFO5<-data.frame(yt[["LFO5"]][["results"]])
yt_LFO10<-data.frame(yt[["LFO10"]][["results"]])
yt_dat <- data.frame(read.csv("Data/Yellowtail/yt_fulldataset_STANDARDIZED.csv"))%>%
  select(-c(X, hci2_pjuv,hci1_pjuv, lusi_annual, hci2_larv,hci1_larv))%>%
  filter(type=="Main")%>%
  filter(Datatreatment=="2025 Final"&year>1993)

#### Yellowtail ####
yt_covariates<-colnames(yt_dat%>%select(-c(Y_rec, sd, type, Datatreatment, year)))
yt_baseline <-  null_RMSE(yt_dat)
yt_model <- gam(Y_rec~1,data=yt_dat)
yt_marginals <- RMSE_improvement(yt_loo,yt_baseline,yt_model,yt_covariates)%>%
  mutate(species=yt_loo$species[1], RMSE="LOO")%>%
  add_row(RMSE_improvement(yt_LFO5,yt_baseline,yt_model,yt_covariates)%>%
            mutate(species=yt_LFO5$species[1],RMSE="LFO 5"))%>%
  add_row(RMSE_improvement(yt_LFO10,yt_baseline,yt_model,yt_covariates)%>%
            mutate(species=yt_LFO10$species[1],RMSE="LFO 10"))
yt_marginals$total_rmse <- apply(yt_marginals[,c("rmse_01","rmse_12","rmse_23")], 1, mean)
yt_marginals$total_aic <- apply(yt_marginals[,c("aic_12", "aic_23")], 1, mean)

yt_marginals<- yt_marginals%>%mutate(RMSEnull=yt_baseline)%>%
  dplyr::mutate(total_rmse_st=total_rmse/3)
relevel_yt_loo<-yt_marginals%>%
  filter(RMSE=="LOO")%>%
  arrange(total_rmse_st)
  
yt_marginals$cov <- factor(yt_marginals$cov, levels =relevel_yt_loo$cov)
#### Sablefish ####
sb_covariates<-colnames(sb_dat%>%select(-c(Y_rec, sd, type,  year)))

sb_baseline <-  null_RMSE(sb_dat)
sb_model <- gam(Y_rec~1,data=sb_dat)
sb_marginals <- RMSE_improvement(sb_loo,sb_baseline,sb_model,sb_covariates)%>%
  mutate(species=sb_loo$species[1], RMSE="LOO")%>%
  add_row(RMSE_improvement(sb_LFO5,sb_baseline,sb_model,sb_covariates)%>%
    mutate(species=sb_LFO5$species[1],RMSE="LFO 5"))%>%
  add_row(RMSE_improvement(sb_LFO10,sb_baseline,sb_model,sb_covariates)%>%
    mutate(species=sb_LFO10$species[1],RMSE="LFO 10"))

sb_marginals$total_rmse <- apply(sb_marginals[,c("rmse_01","rmse_12","rmse_23")], 1, mean)
sb_marginals$total_aic <- apply(sb_marginals[,c("aic_12", "aic_23")], 1, mean)

sb_marginals<- sb_marginals%>%mutate(RMSEnull=sb_baseline)%>%
  dplyr::mutate(total_rmse_st=total_rmse/3)

relevel_sb_loo<-sb_marginals%>%
  filter(RMSE=="LOO")%>%
  arrange(total_rmse_st)
sb_marginals$cov <- factor(sb_marginals$cov, levels =relevel_sb_loo$cov)



#### Petrale Sole ####
ps_covariates<-colnames(ps_dat%>%select(-c(Y_rec, sd, type,  year)))
ps_baseline <-  null_RMSE(ps_dat)
ps_model <- gam(Y_rec~1,data=ps_dat)
ps_marginals <- RMSE_improvement(ps_loo,ps_baseline,ps_model,ps_covariates)%>%
  mutate(species=ps_loo$species[1], RMSE="LOO")%>%
  add_row(RMSE_improvement(ps_LFO5,ps_baseline,ps_model,ps_covariates)%>%
            mutate(species=ps_LFO5$species[1],RMSE="LFO 5"))%>%
  add_row(RMSE_improvement(ps_LFO10,ps_baseline,ps_model,ps_covariates)%>%
            mutate(species=ps_LFO10$species[1],RMSE="LFO 10"))
ps_marginals$total_rmse <- apply(ps_marginals[,c("rmse_01","rmse_12","rmse_23")], 1, mean)
ps_marginals$total_aic <- apply(ps_marginals[,c("aic_01","aic_12", "aic_23")], 1, mean)

ps_marginals<- ps_marginals%>%mutate(RMSEnull=ps_baseline)%>%
  dplyr::mutate(total_rmse_st=total_rmse/3)

relevel_ps_loo<-ps_marginals%>%
  filter(RMSE=="LOO")%>%
  arrange(total_rmse_st)
ps_marginals$cov <- factor(ps_marginals$cov, levels =relevel_ps_loo$cov)

#### Single dataset ####
marginals<-ps_marginals%>%
  add_row(sb_marginals)%>%
  add_row(yt_marginals)
write_rds(marginals, "Output/Data/marginals.rds")

#### Figures #### 
cols<- c('#dd4124',"#edd746",'#7cae00','#0f85a0')

marginal_yt <- ggplot(yt_marginals, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y =cov,
  x = total_rmse_st,
  fill = total_rmse_st
)) +
  xlim(c(-.1,0.075))+
 facet_grid(species ~ RMSE, scales = "free_y") + # Use free_y scale
  geom_bar(stat = "identity") +
  labs(x = "Mean Marginal Improvement RMSE", y = "") +
  scale_fill_gradient(low = "gray100", high = "Darkgreen",limits=c(-.1,0.075)) +
  theme_classic()+
  theme(legend.position = "none")
marginal_yt

marginal_sb <- ggplot(sb_marginals, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y =cov,
  x = total_rmse_st,
  fill = total_rmse_st
)) +
  xlim(c(-.1,0.075))+
  facet_grid(species ~ RMSE, scales = "free_y") + # Use free_y scale
  geom_bar(stat = "identity") +
  labs(x = "", y = "Predictor") +
  scale_fill_gradient(low = "gray100", high = "Darkgreen", ,limits=c(-.1,0.075)) +
  theme_classic()
marginal_sb

marginal_ps <- ggplot(ps_marginals, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y =cov,
  x = total_rmse_st,
  fill = total_rmse_st
)) +
  xlim(c(-.1,0.075))+
  facet_grid(species ~ RMSE, scales = "free_y") + # Use free_y scale
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  scale_fill_gradient(low = "gray100", high = "Darkgreen",limits=c(-.1,0.075)) +
  theme_classic()+
  theme(legend.position = "none")
marginal_ps 

z<- ggplot()+theme_void()

marginal<- ggarrange(ggarrange(marginal_ps,z, widths = c(4.25,0.75), ncol=2, nrow=1),
                     marginal_sb, 
                     ggarrange(marginal_yt,z, widths = c(4.25,0.75), ncol=2, nrow=1), 
                     ncol = 1, nrow = 3)
marginal


pdf(file = "Output/Figures/MarginalMeanRMSE.pdf", width = 8, height = 11)
marginal 
dev.off()

png(file = "Output/Figures/MarginalMeanRMSE.png",width = 800, height = 1100, res = 100)
marginal 
dev.off()

#### Rolling Window ####

##### yellowtail #####
yt_RW<-data.frame(yt[["RW"]][["results"]])
length_window<-15
firstyr<-unique(yt_RW$firstyear)
lastyr<-unique(yt_RW$lastyear)
dat<-yt_RW
results <- data.frame()
results_temp <- data.frame()
for(k in 1:length(firstyr)){
  datwindow <- dat%>%filter(firstyear==firstyr[k])
  results_temp<-  RMSE_improvement(datwindow,yt_baseline,yt_model,yt_covariates)%>%
    mutate(range=paste(firstyr[k],"-",lastyr[k]), FirstYear=firstyr[k], LastYear=lastyr[k])
 
  results<- rbind(results,results_temp)   
}
rolling_yt<-results
rolling_yt$total_rmse <- apply(rolling_yt[,c("rmse_12","rmse_23")], 1, mean)/3
rolling_yt$total_aic <- apply(rolling_yt[,c("aic_12", "aic_23")], 1, mean)/3

rollingplot_yt<-ggplot(rolling_yt,aes(x=as.factor(LastYear), y=cov, fill= total_rmse )) + 
  xlab("First Year")+
  scale_fill_gradient(low = "white", high = "Darkgreen") +
  ylab("Oceanographic Conditions")+
  ggtitle("Yellowtail")+
  geom_tile()+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8),plot.title = element_text(hjust = 0.5))
rollingplot_yt
##### sablefish #####
sb_RW<-data.frame(sb[["RW"]][["results"]])
length_window<-15
firstyr<-unique(sb_RW$firstyear)
lastyr<-unique(sb_RW$lastyear)
dat<-sb_RW
results <- data.frame()
results_temp <- data.frame()
for(k in 1:length(firstyr)){
  datwindow <- dat%>%filter(firstyear==firstyr[k])
  results_temp<-  RMSE_improvement(datwindow,sb_baseline,sb_model,sb_covariates)%>%
    mutate(range=paste(firstyr[k],"-",lastyr[k]), FirstYear=firstyr[k], LastYear=lastyr[k])
  
  results<- rbind(results,results_temp)   
}

rolling_sb<-results
rolling_sb$total_rmse <- apply(rolling_sb[,c("rmse_12","rmse_23")], 1, mean)/3
rolling_sb$total_aic <- apply(rolling_sb[,c("aic_12", "aic_23")], 1, mean)/3

rollingplot_sb<-ggplot(rolling_sb,aes(x=as.factor(LastYear), y=cov, fill= total_rmse)) + 
  xlab("First Year")+
  scale_fill_gradient(low = "white", high = "Darkgreen") +
  ylab("Oceanographic Conditions")+
  ggtitle("Sablefish")+
  geom_tile()+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8),plot.title = element_text(hjust = 0.5))


##### sablefish #####
ps_RW<-data.frame(ps[["RW"]][["results"]])
firstyr<-unique(ps_RW$firstyear)
lastyr<-unique(ps_RW$lastyear)
dat<-ps_RW
results <- data.frame()
results_temp <- data.frame()
for(k in 1:length(firstyr)){
  datwindow <- dat%>%filter(firstyear==firstyr[k])
  results_temp<-  RMSE_improvement(datwindow,ps_baseline,ps_model,ps_covariates)%>%
    mutate(range=paste(firstyr[k],"-",lastyr[k]),FirstYear=firstyr[k],LastYear=lastyr[k])
  
  results<- rbind(results,results_temp)   
}
rolling_ps<-results
rolling_ps$total_rmse <- apply(rolling_ps[,c("rmse_12","rmse_23")], 1, mean)/3
rolling_ps$total_aic <- apply(rolling_ps[,c("aic_12", "aic_23")], 1, mean)/3

rollingplot_ps<-ggplot(rolling_ps,aes(x=as.factor(LastYear), y=cov, fill= total_rmse)) + 
  xlab("First Year")+
  scale_fill_gradient(low = "white", high = "Darkgreen") +
  ylab("Oceanographic Conditions")+
  ggtitle("Petrale Sole")+
  geom_tile()+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8),plot.title = element_text(hjust = 0.5))

##### Figures #####
rollingplot <- ggarrange(rollingplot_ps, rollingplot_sb, rollingplot_yt, ncol = 2, nrow = 2)

pdf(file = "Output/Figures/rollingplot.pdf", width =12, height = 8)
rollingplot
dev.off()

png(file = "Output/Figures/rollingplot.png",width = 1200, height = 800, res = 100)
rollingplot
dev.off()
