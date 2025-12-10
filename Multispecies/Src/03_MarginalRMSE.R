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

#### Yellowtail ####

yt_baseline <-  null_RMSE(yt_dat)
yt_model <- gam(Y_rec~1,data=yt_dat)
yt_marginals <- RMSE_improvement(yt_loo,yt_baseline,yt_model,yt_covariates)%>%
  mutate(species=yt_loo$species[1], RMSE="LOO")%>%
  add_row(RMSE_improvement(yt_LFO5,yt_baseline,yt_model,yt_covariates)%>%
            mutate(species=yt_LFO5$species[1],RMSE="LFO 5"))%>%
  add_row(RMSE_improvement(yt_LFO10,yt_baseline,yt_model,yt_covariates)%>%
            mutate(species=yt_LFO10$species[1],RMSE="LFO 10"))

yt_marginals$total_rmse <- apply(yt_marginals[,c("rmse_12","rmse_23")], 1, mean)
yt_marginals$total_aic <- apply(yt_marginals[,c("aic_12", "aic_23")], 1, mean)

yt_marginals<- yt_marginals%>%mutate(RMSEnull=yt_baseline)
#### Sablefish ####

sb_baseline <-  null_RMSE(sb_dat)
sb_model <- gam(Y_rec~1,data=sb_dat)
sb_marginals <- RMSE_improvement(sb_loo,sb_baseline,sb_model,sb_covariates)%>%
  mutate(species=sb_loo$species[1], RMSE="LOO")%>%
  add_row(RMSE_improvement(sb_LFO5,sb_baseline,sb_model,sb_covariates)%>%
    mutate(species=sb_LFO5$species[1],RMSE="LFO 5"))%>%
  add_row(RMSE_improvement(sb_LFO10,sb_baseline,sb_model,sb_covariates)%>%
    mutate(species=sb_LFO10$species[1],RMSE="LFO 10"))

sb_marginals$total_rmse <- apply(sb_marginals[,c("rmse_12","rmse_23")], 1, mean)
sb_marginals$total_aic <- apply(sb_marginals[,c("aic_12", "aic_23")], 1, mean)

sb_marginals<- sb_marginals%>%mutate(RMSEnull=sb_baseline)
#### Petrale Sole ####

ps_baseline <-  null_RMSE(ps_dat)
ps_model <- gam(Y_rec~1,data=ps_dat)
ps_marginals <- RMSE_improvement(ps_loo,ps_baseline,ps_model,ps_covariates)%>%
  mutate(species=ps_loo$species[1], RMSE="LOO")%>%
  add_row(RMSE_improvement(ps_LFO5,ps_baseline,ps_model,ps_covariates)%>%
            mutate(species=ps_LFO5$species[1],RMSE="LFO 5"))%>%
  add_row(RMSE_improvement(ps_LFO10,ps_baseline,ps_model,ps_covariates)%>%
            mutate(species=ps_LFO10$species[1],RMSE="LFO 10"))
ps_marginals$total_rmse <- apply(ps_marginals[,c("rmse_12","rmse_23")], 1, mean)
ps_marginals$total_aic <- apply(ps_marginals[,c("aic_12", "aic_23")], 1, mean)

ps_marginals<- ps_marginals%>%mutate(RMSEnull=ps_baseline)
#### Single dataset ####
marginals<-ps_marginals%>%
  add_row(sb_marginals)%>%
  add_row(yt_marginals)
write_rds(marginals, "Output/Data/marginals.rds")

#### Figures #### 

gam_loo_table <- marginals%>%
  group_by(species, RMSE)%>%
  dplyr::arrange(marginals, rmse_23)%>%
  dplyr::select(cov, rmse_23, total_rmse)
gam_loo_table[is.na(gam_loo_table)] <-"No"
cols<- c('#dd4124',"#edd746",'#7cae00','#0f85a0')

marginal <- ggplot(gam_loo_table, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y = reorder_within(cov, total_rmse, species),
  x = total_rmse,
  fill = total_rmse
)) +
  # Note that we swapped x and y in aes() because of coord_flip()
  facet_grid(species ~ RMSE, scales = "free_y") + # Use free_y scale
  geom_bar(stat = "identity") +
  # Add the scale_y_reordered function to clean up the axis labels
 scale_y_reordered() +
  #coord_flip() +
  labs(x = "Mean Marginal Improvement RMSE", y = "Predictor") +
  scale_fill_gradient(low = "white", high = "Darkgreen") +
  theme_classic()
marginal
pdf(file = "Output/Figures/MarginalMeanRMSE.pdf", width = 11, height = 8)
marginal 
dev.off()

png(file = "Output/Figures/MarginalMeanRMSE.png",width = 1100, height = 800, res = 100)
marginal 
dev.off()

