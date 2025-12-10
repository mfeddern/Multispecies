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
  select(-c(X, hci2_pjuv,hci1_pjuv, lusi_annual, hci2_larv,hci1_larv))%>%
    filter(type=="Main")%>%
  filter(Datatreatment=="2025 Final"&year>1993)
#yt_loo <- readRDS("Output/Data/yt_selection.rds")

hk_dat <- data.frame(read.csv("Data/Hake/DATA_Combined_glorys_hake_STANDARDIZED.csv"))%>%
  select(-X)%>%
  filter(type=="Main_RecrDev"&year>1993)

ps_dat <- data.frame(read.csv("Data/Petrale/DATA_Combined_glorys_petrale_STANDARDIZED.csv"))%>%
  select(-c(X, year.1))%>%
  filter(type=="Main_RecrDev"&year>1993)

sb_dat <- data.frame(read.csv("Data/Sablefish/data-combined-glorys-sablefish_STANDARDIZED.csv"))%>%
  select(-c(X))%>%
  filter(type=="Main_RecrDev"&year>1993)

#### Functions ####
threshold<-3
combinations <- function(data,threshold,maxcovar){
  #threshold <- corr_threshold #0.95#0.3 #assiging a threshold of correlation to filter out 
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
  # <- 4 #max number of covars for a given model
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
    TP= sum(dat$Y_rec >0 & predictions >0)
    TN= sum(dat$Y_rec <0 & predictions <0)
    FP= sum(dat$Y_rec >0 & predictions <0)
    FN= sum(dat$Y_rec <0 & predictions >0)
    Hit <- (TP + TN)/((TP+TN)+(FP+FN))
    # Store results
    models[[i]] <- gam_model
    results <- rbind(results, data.frame(
      ModelID = i,
      AIC = AIC(gam_model),
      RMSE_loo = round(rmse_loo,3),
      RMSE = round(rmse,3),
      rsq_full=round(r2,2),
      dev.ex=round(dev.expl,4),
      Hit=Hit,
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
rw_model_fit <- function(combinations, dataset,yrlast, window,species_name){
  models <- list()
  results <- data.frame()
  predicted<-data.frame()
  jstart<-1
  #n_year <-lengthwindow
  length_window<-window
  firstyear<-1994:(yrlast-length_window)
  for(k in jstart:length(firstyear)){
    lastyear=firstyear[k]+length_window
    dat <- dataset%>%filter(year>=firstyear[k]&year<=lastyear)
    jstart<-1
    results_temp <- data.frame()
    predicted_temp <- data.frame()
    models_temp <- list()
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
      #models_temp[[i]] <- gam_model
      results_temp<- rbind(results_temp,data.frame(
        ModelID = i,
        species=species_name,
        AIC = AIC(gam_model),
        RMSE_loo = round(rmse_loo,3),
        RMSE = round(rmse,3),
        rsq_full=round(r2,2),
        dev.ex=round(dev.expl,4),
        firstyear=firstyear[k],
        lastyear=lastyear,
        var1 = padded_vars[1],
        var2 = padded_vars[2],
        var3 = padded_vars[3],
        var4 = padded_vars[4]
        
        
      ))
      
      predicted_temp <- rbind(predicted_temp, data.frame(
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
    results<- rbind(results,results_temp)
    predicted<- rbind(predicted,predicted_temp)
    #models <- list(models, models_temp)
  }
  results_output <- list("results" = results, "predicted" = predicted, "model"=models)
  return(results_output)
  
}
null_RMSE <-function(dat){
  null<-dat%>%mutate(sr_null = 0)
  sqerror<-(null$sr_null-null$Y_rec)^2
  rmse_sr_full <- sqrt(mean(sqerror, na.rm=T))
  return(rmse_sr_full)
}
LFO <- function(dat,combinations,n_pred, species_name){
  n_year <-length(unique(dat$year))
  n_train <- n_year-n_pred 
  kstart<-n_train+1
  models<-list()
  results<-data.frame()
  predicted<-data.frame()
  for (i in seq_along(combinations)) {
    # ps here represents a P-spline / penalized regression spline
    # k represent the number of parameters / knots estimating function at, should be small
    
    # smooth on total release isn't needed when we're not working with jack rate
    #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
    smooth_terms <- paste("s(", combinations[[i]], ", k = 3)")
    formula_str <- paste("Y_rec ~ ", smooth_terms)
    predictions <-  numeric(n_pred )
    
    # Loop over each observation
    for (k in kstart:n_year) {
      train_index <- setdiff(train_start:n_train, k)  # All indices except the k-th
      test_index <- k                 # The k-th index
      
      # Fit model on n-k observations
      gam_model <- gam(as.formula(formula_str),
                       data = dat[which(dat$year %notin% unique(dat$year)[k:n_year]), ])
      
      # Predict the excluded observation
      predictions[which(dat$year == unique(dat$year)[k])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[k]), ])
      # re-fit the model
      
    }
    # re-fit the model
    gam_model <- gam(as.formula(formula_str),
                     #data = dat)
                     data = dat[which(dat$year %notin% unique(dat$year)[kstart:n_year]), ])
    gam_model2 <- gam(as.formula(formula_str),
                      data = dat)
    # keep in mind RMSE is weird for binomial things
    rmse <- sqrt(mean((dat$Y_rec[kstart:n_year] - predictions[kstart:n_year])^2, na.rm=T))
    r2<-summary(gam_model)$r.sq
    dev.expl<-summary(gam_model)$dev.expl
    # Extract variable names
    var_names <- gsub("s\\(([^,]+),.*", "\\1", combinations[[i]])
    padded_vars <- c(var_names, rep(NA, 4 - length(var_names)))
    TP= sum(dat$Y_rec[kstart:n_year] >0 & predictions[kstart:n_year] >0)
    TN= sum(dat$Y_rec[kstart:n_year] <0 & predictions[kstart:n_year] <0)
    FP= sum(dat$Y_rec[kstart:n_year] >0 & predictions[kstart:n_year] <0)
    FN= sum(dat$Y_rec[kstart:n_year] <0 & predictions[kstart:n_year] >0)
    Hit <- (TP + TN)/((TP+TN)+(FP+FN))
    
    # Store results
    models[[i]] <- gam_model
    results <- rbind(results, data.frame(
      ModelID = i,
      species=species_name,
      Hit=Hit,
      AIC = AIC(gam_model),
      RMSE = round(rmse,3),
      rsq=round(r2,2),
      dev.ex=round(dev.expl,4),
      n_pred=n_pred,
      #AUC = auc,
      #direction = direction,
      var1 = padded_vars[1],
      var2= padded_vars[2],
      var3 = padded_vars[3],
      var4 = padded_vars[4]
      
    ))
    predicted <- rbind(predicted, data.frame(
      ModelID = i,
      species=species_name,
      pred=predictions[kstart:n_year],
      var1 = padded_vars[1],
      var2 = padded_vars[2],
      var3 = padded_vars[3],
      var4 = padded_vars[4]
      
    ))
    
    print(i)
  }
  results_output <- list("results" = results, "predicted" = predicted)
  
  return(results_output)
}
model_fit <- function(combinations, dat,species_name){
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
    TP= sum(dat$Y_rec >0 & predictions >0)
    TN= sum(dat$Y_rec <0 & predictions <0)
    FP= sum(dat$Y_rec >0 & predictions <0)
    FN= sum(dat$Y_rec <0 & predictions >0)
    Hit <- (TP + TN)/((TP+TN)+(FP+FN))
    # Store results
    models[[i]] <- gam_model
    results <- rbind(results, data.frame(
      ModelID = i,
      species=species_name,
      AIC = AIC(gam_model),
      RMSE_loo = round(rmse_loo,3),
      RMSE = round(rmse,3),
      rsq_full=round(r2,2),
      dev.ex=round(dev.expl,4),
      Hit=Hit,
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

train_start<-1
`%notin%` <- Negate(`%in%`)

#### Yellowtail ####

yt_combinations_results <- combinations(yt_dat%>%select(-Datatreatment), 0.5,3)
yt_combinations<- yt_combinations_results$combinations
yt_covariates<- yt_combinations_results$covariates
yt_results_5<- LFO(yt_dat,yt_combinations, 5,"Yellowtail")
yt_results_10<- LFO(yt_dat,yt_combinations, 10,"Yellowtail")
yt_results_loo<- model_fit(yt_combinations, yt_dat,"Yellowtail")
yt_results_rw<- rw_model_fit(yt_combinations, yt_dat,2020, 15,"Yellowtail")
yt_results <- list("LFO5" = yt_results_5, "LFO10" = yt_results_10,"LOO" = yt_results_loo,"RW" = yt_results_rw)

write_rds(yt_results, "Output/Data/yt_model_fits.rds")

#### Sablefish ####

sb_combinations_results <- combinations(sb_dat,0.5,3)
sb_combinations<- sb_combinations_results$combinations
sb_covariates<- sb_combinations_results$covariates
sb_results_5<- LFO(sb_dat,sb_combinations, 5,"Sablefish")
sb_results_10<- LFO(sb_dat,sb_combinations, 10,"Sablefish")
sb_results_loo<- model_fit(sb_combinations, sb_dat,"Sablefish")
sb_results_rw<- rw_model_fit(sb_combinations, sb_dat,2020, 15,"Sablefish")
sb_results <- list("LFO5" = sb_results_5, "LFO10" = sb_results_10,"LOO" = sb_results_loo,"RW" = sb_results_rw)

write_rds(sb_results, "Output/Data/sb_model_fits.rds")


#### Petrale Sole ####

ps_combinations_results <- combinations(ps_dat,0.5,3)
ps_combinations<- ps_combinations_results$combinations
ps_covariates<- ps_combinations_results$covariates
ps_results_5<- LFO(ps_dat,ps_combinations, 5,"Petrale Sole")
ps_results_10<- LFO(ps_dat,ps_combinations, 10,"Petrale Sole")
ps_results_loo<- model_fit(ps_combinations, ps_dat,"Petrale Sole")
ps_results_rw<- rw_model_fit(ps_combinations, ps_dat,2020, 15,"Petrale Sole")
ps_results <- list("LFO5" = ps_results_5, "LFO10" = ps_results_10,"LOO" = ps_results_loo,"RW" = ps_results_rw)

write_rds(ps_results, "Output/Data/ps_model_fits.rds")

