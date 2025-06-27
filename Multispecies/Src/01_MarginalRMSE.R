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

#### Reading in the data ####
yt_dat <- data.frame(read.csv("Data/Yellowtail/yt_fulldataset_STANDARDIZED.csv"))%>%
  select(-X)%>%
  filter(Datatreatment=="Expanded PacFIN"&type=="Main_RecrDev"&year>1993)
  
hk_dat <- data.frame(read.csv("Data/Hake/DATA_Combined_glorys_hake_STANDARDIZED.csv"))%>%
  select(-X)%>%
  filter(type=="Main_RecrDev"&year>1993)

ps_dat <- data.frame(read.csv("Data/Petrale/DATA_Combined_glorys_petrale_STANDARDIZED.csv"))%>%
  select(-X)%>%
  filter(type=="Main_RecrDev"&year>1993)

sb_dat <- data.frame(read.csv("Data/Sablefish/data-combined-glorys-sablefish_STANDARDIZED.csv"))%>%
  select(-X)%>%
  filter(type=="Main_RecrDev"&year>1993)

#### Functions ####
combinations <- function(data){
threshold <- 0.3 #0.95#0.3 #assiging a threshold of correlation to filter out 
yt_env <- data%>%select(-c(Y_rec, sd, type))
yt_env <- yt_env[complete.cases(yt_env), ] %>% 
  dplyr::select(!any_of( c('year'))) # getting environmental data
M = data.frame(cor(yt_env)) # generating a correlation matrix
M <- tibble::rownames_to_column(M, "xvar") #putting the rownames in their own column
M<-M%>%pivot_longer(!xvar, names_to = 'yvar', values_to = 'corr') #pivoting to a longer table
uncorr<-M%>%filter(abs(corr)<threshold) #generating the uncorrelated thresholds 
corrvar<-M%>%filter(abs(corr)>threshold)%>%filter(xvar!=yvar) #generating the correlated thresholds
corrvar2<-corrvar%>%select(-corr)
combinations_to_omit<-list() #setting up an empty list to fill with pairs to omit
for(i in 1:nrow(corrvar2)){
  combinations_to_omit[[i]]<-c(as.character(corrvar2[i, ]))
}

##### Models to Fit #####
covariates <- names(yt_env)
maxcovar <- 4 #max number of covars for a given model
combinations <- lapply(1:maxcovar, function(i) {
  combn(covariates, i, simplify = FALSE)
})

#setting up all possible combinations
combinations <- unlist(combinations, recursive = FALSE) 
combinations <- unique(lapply(combinations, function(x) sort(x)))
length(combinations) #all possible combinations 15275

# Function to check if a combination is a partial match of any combination to omit
is_partial_match <- function(comb, omit_list) {
  any(sapply(omit_list, function(omit) all(omit %in% comb)))
}

# Remove combinations that are partial matches (but not exact matches) to the combinations to omit
combinations <- combinations[!sapply(combinations, is_partial_match, omit_list = combinations_to_omit)]
# Check the length of remaining combinations
length(combinations) # how many left? only 887 (best model does not change compared to threshold of 0.4 but reduces model by 1/2)
combs <- list("combinations" = combinations, "covariates" = covariates)
return(combs)

}
model_fit <- function(combinations, dat){
  models <- list()
  results <- data.frame()
  predicted <- data.frame()
  jstart<-1
  for (i in seq_along(combinations)) {
    # k represent the number of parameters / knots estimating function at, should be small
    #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
    smooth_terms <- paste("s(", combinations[[i]], ", k = 3)", collapse = " + ")
    formula_str <- paste("Y_rec ~ ", smooth_terms)
    predictions <- numeric(nrow(dat))
    n_year <- length(unique(dat$year))
    # Loop over each observation
    for (j in jstart:n_year) {
      train_index <- setdiff(jstart:n_year, j)  # All indices except the j-th
      test_index <- j                 # The j-th index
      
      # Fit model on n-1 observations
      gam_model <- gam(as.formula(formula_str),
                       # weights = number_cwt_estimated,
                       data = dat[which(dat$year != unique(dat$year)[j]), ])
      
      # Predict the excluded observation
      predictions[which(dat$year == unique(dat$year)[j])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[j]), ])
    }
    
    # re-fit the model
    gam_model <- gam(as.formula(formula_str),
                     data = dat)
    predict_mod <- predict(gam_model)
    # keep in mind RMSE is weird for binomial things
    rmse_loo <- sqrt(mean((dat$Y_rec - predictions)^2, na.rm=T))
    rmse <- sqrt(mean((dat$Y_rec -  predict_mod)^2, na.rm=T))
    r2<-summary(gam_model)$r.sq
    dev.expl<-summary(gam_model)$dev.expl
    # Extract variable names
    var_names <- gsub("s\\(([^,]+),.*", "\\1", combinations[[i]])
    # Store results with variable names padded to ensure there are always 3 columns
    padded_vars <- c(var_names, rep(NA, 4 - length(var_names)))
    
    # Store results
    models[[i]] <- gam_model
    results <- rbind(results, data.frame(
      ModelID = i,
      AIC = AIC(gam_model),
      RMSE_loo = round(rmse_loo,3),
      RMSE = round(rmse,3),
      rsq_full=round(r2,2),
      dev.ex=round(dev.expl,4),
      #rmse_imp=(rmse_sr_full-rmse)/rmse_sr_full, 
      #rmse_ratio=rmse/rmse_sr_full,
      #AUC = auc,
      #direction = direction,
      var1 = padded_vars[1],
      var2 = padded_vars[2],
      var3 = padded_vars[3],
      var4 = padded_vars[4]
      
      
    ))
    
    predicted <- rbind(predicted, data.frame(
      ModelID = i,
      pred=predictions[jstart:n_year],
      year=unique(dat$year)[jstart:n_year],
      var1 = padded_vars[1],
      var2 = padded_vars[2],
      var3 = padded_vars[3],
      var3 = padded_vars[4]
      
    ))
    #saving the one step ahead predictions
    print(i)
  }
  results_output <- list("results" = results, "predicted" = predicted)
  
  return(results_output)
}
RMSE_improvement <-function(results,baseline_rmse,gam_model){
  baseline_rmse <-0.4
  
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

#### Yellowtail ####

yt_combinations_results <- combinations(yt_dat%>%select(-Datatreatment))
yt_combinations<- yt_combinations_results$combinations
yt_covariates<- yt_combinations_results$covariates
yt_results<- model_fit(yt_combinations, yt_dat)
yt_selection <- data.frame(yt_results$results)
yt_predicted <- yt_results$predicted
yt_baseline <- 0.4
yt_model <- gam(Y_rec~1,data=yt_dat)
yt_marginals <- RMSE_improvement(yt_selection,yt_baseline,yt_model)

yt_marginals$total_rmse <- apply(yt_marginals[,c("rmse_12","rmse_23")], 1, mean)
yt_marginals$total_aic <- apply(marginals[,c("aic_12", "aic_23")], 1, mean)

gam_loo_table <- dplyr::arrange(yt_marginals, rmse_23)%>%
  dplyr::select(cov, rmse_23, total_aic)%>%
  mutate(cv="LOO",model="GAM")
gam_loo_table[is.na(gam_loo_table)] <-"No"
cols<- c('#dd4124',"#edd746",'#7cae00','#0f85a0')

marginal <- ggplot(gam_loo_table, aes(x = reorder(cov,rmse_23), y = rmse_23)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip the axes to make a horizontal bar graph
  labs(y = "Marginal Improvement RMSE", x = "Predictor")+
  scale_fill_manual(values=c('grey',cols[3]))+
  theme_classic()
marginal 

### Calculating Relative Variable Importance ###

# Create a baseline model with only the intercept

hk_combinations <- combinations(hk_dat)
hk_results<-model_fit(hk_combinations[1:20], hk_dat)