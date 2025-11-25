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
#  filter(type=="Main")%>%
  filter(Datatreatment=="2025 Final"&year>1993)
  
yt_loo <- readRDS("Output/Data/yt_selection.rds")

### Functions ###

FunctionalRelationships<- function(data,peel){
for(i in 1:peel){
  trainingmin <-nyears - peel
  peeled<-trainingmin+ i
  trainingdata<-  data[1:(nyears-i),]
  mod<-gam(as.formula(formula_str), data= trainingdata)
  x<- plot(mod)
  ts_smooths_temp<-data.frame()
  temp_edfs<- summary(mod)$s.table[, "edf"]
for(j in 1:numvar){ #loop over var
#this section is getting functional relationship for each var 
  results_temp<- data.frame(c(x[[j]]$fit),
  x[[j]]$x,
  x[[j]]$se,
  paste(x[[j]]$xlab)
  )
  #colnames(results_temp)<-c(paste(x[[j]]$xlab,"fit"), paste(x[[j]]$xlab,"x"),paste(x[[j]]$xlab,"se"))
  colnames(results_temp)<-c("fit", "x","se","var")
  ts_smooths_temp<-rbind(results_temp, ts_smooths_temp)
#this section is getting the calculations for Mohns

}
  xyy<-data.frame(mohns=temp_edfs-full_edfs, termyr=data$year[(nyears-i)])
  mohns_temp<- cbind(xyy, row.names(xyy))
  mohns<-rbind(mohns_temp,mohns)
  #ts_smooths[[i]]<-ts_smooths_temp
  ts_smooths_temp<-data.frame(termyr=data$year[(nyears-i)],ts_smooths_temp)
  ts_smooths<- rbind(ts_smooths,ts_smooths_temp)
  print(data$year[(nyears-i)])
}
  results<-list(ts_smooths, mohns)
  return(results)
}

### Yellowtail ###
modtemp <- arrange(yt_loo,RMSE_loo)[1,]%>%#select best model
  select(var1, var2, var3, var4) #select variables. Way to do this with all 4 and filter out 4th when it isn't used?
numvar<-sum(!is.na(modtemp))
mod<- modtemp [1:numvar]
smooth_terms <- paste("s(",mod, ", k = 3)", collapse = " + ")
formula_str <- paste("Y_rec ~ ", gsub(" ", "",smooth_terms))
model<-gam(as.formula(formula_str), data=yt_dat)
full_edfs<-summary(model)$s.table[, "edf"]
nyears<-length(yt_dat$year)

ts_smooths <- data.frame()
mohns <- data.frame()
results <- FunctionalRelationships(yt_dat, 15)


ggplot(data=ts_smooths,aes(x=x, y=fit))+ #you could add observations onto this to show contrast
  facet_wrap(~var)+
  geom_line(aes(group=termyr,col=termyr))
