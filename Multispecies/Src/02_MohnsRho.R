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
yt <- readRDS("Output/Data/yt_model_fits.rds")
yt_loo<-data.frame(yt[["LOO"]][["results"]])
yt_lfo5<-data.frame(yt[["LFO5"]][["results"]])
yt_lfo10<-data.frame(yt[["LFO10"]][["results"]])

sb_dat <- data.frame(read.csv("Data/Sablefish/data-combined-glorys-sablefish_STANDARDIZED.csv"))%>%
  select(-c(X))%>%
  filter(type=="Main_RecrDev"&year>1993)
sb <- readRDS("Output/Data/sb_model_fits.rds")
sb_loo<-data.frame(sb[["LOO"]][["results"]])
sb_lfo5<-data.frame(sb[["LFO5"]][["results"]])
sb_lfo10<-data.frame(sb[["LFO10"]][["results"]])

ps_dat <- data.frame(read.csv("Data/Petrale/DATA_Combined_glorys_petrale_STANDARDIZED.csv"))%>%
  select(-c(X, year.1))%>%
  filter(type=="Main_RecrDev"&year>1993)
ps <- readRDS("Output/Data/ps_model_fits.rds")
ps_loo<-data.frame(ps[["LOO"]][["results"]])
ps_lfo5<-data.frame(ps[["LFO5"]][["results"]])
ps_lfo10<-data.frame(ps[["LFO10"]][["results"]])

### Functions ###
par(mfrow = c(2, 2))

FunctionalRelationships<- function(data,peel,form,numvar){
  model<-gam(as.formula(form), data= data)
  full_edfs<-summary(model)$s.table[, "edf"]
  mohns<-data.frame()
  ts_smooths<-data.frame()
  nyears<-length(data)
for(i in 1:peel){

  trainingmin <-nyears - peel
  peeled<-trainingmin+ i
  trainingdata<-  data[1:(nyears-i),]
  mod<-gam(as.formula(form), data= trainingdata)
  x<- plot(mod,ask="N")
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
  xyy<-data.frame(mohns=((temp_edfs-full_edfs)/full_edfs), termyr=data$year[(nyears-i)])
  ts_smooths_temp<-data.frame(termyr=data$year[(nyears-i)],ts_smooths_temp)
  mohns_temp<- cbind(xyy, row.names(xyy))
  mohns<-rbind(mohns_temp,mohns)
  #ts_smooths[[i]]<-ts_smooths_temp
  ts_smooths<- rbind(ts_smooths,ts_smooths_temp)
  print(data$year[(nyears-i)])
}
  results<-list(ts_smooths, mohns)
  return(results)
}

single_covs<- function(models, data){
  nmodels<-length(unique(models$var1))
  nyears<-length(data$year)

  for(i in 1:nmodels){
    modtemp <- models[i,]%>%#select best model
      select(var1, var2, var3, var4) #select variables. Way to do this with all 4 and filter out 4th when it isn't used?
    numvar<-sum(!is.na(modtemp))
    mod2<- modtemp [1:numvar]
    smooth_terms <- paste("s(",mod2, ", k = 3)", collapse = " + ")
    formula_str <- paste("Y_rec ~ ", gsub(" ", "",smooth_terms))
    formulas<- rbind(formulas, formula_str)

  }
  
  for(j in 1:nmodels){
    results<-FunctionalRelationships(data, npeel,formulas[[j]],1)
    res_temp<-results[[1]]
    mohn_temp<-results[[2]]%>%
      add_column(variable =models$var1[j])
    res<-rbind(res,res_temp)
    mohn<-rbind(mohn,mohn_temp)
  }
  single_results<-list(res, mohn)
 return(single_results)
}

#### Best Models ####
##### Yellowtail #####
modtemp <- arrange(yt_lfo5,RMSE)[1,]%>%#select best model
  select(var1, var2, var3, var4)
#modtemp <- arrange(yt_loo,RMSE_loo)[1,]%>%#select best model
#  select(var1, var2, var3, var4) #select variables. Way to do this with all 4 and filter out 4th when it isn't used?
numvar<-sum(!is.na(modtemp))
mod<- modtemp [1:numvar]
smooth_terms <- paste("s(",mod, ", k = 3)", collapse = " + ")
formula_str <- paste("Y_rec ~ ", gsub(" ", "",smooth_terms))
ts_smooths <- data.frame()
mohns <- data.frame()
results_yt <- FunctionalRelationships(yt_dat, 15,formula_str,1)
model<-gam(as.formula(formula_str), data=yt_dat)
full_edfs<-summary(model)$s.table[, "edf"]
nyears<-length(yt_dat$year)

mohns_yt<- results_yt[[2]]%>%
  rename(variable='row.names(xyy)')%>%
  mutate(Species="Yellowtail")

plot_dat_yt<-yt_dat%>%select(year,Y_rec,unique(results_yt[[1]]$var))%>%
  pivot_longer(cols=c(unique(results_yt[[1]]$var)))%>%
  rename(var=name)%>%
  mutate(Species="Yellowtail")

##### Sablefish #####
npeel<-15
formulas<-list()
res<-data.frame()
mohn<-data.frame()
#rec<-single_covs(sb_loo,sb_dat)

modtemp <- arrange(sb_lfo5,RMSE)[1,]%>%#select best model
  select(var1, var2, var3, var4) #select variables. Way to do this with all 4 and filter out 4th when it isn't used?
numvar<-sum(!is.na(modtemp))
mod<- modtemp [1:numvar]
smooth_terms <- paste("s(",mod, ", k = 3)", collapse = " + ")
formula_str <- paste("Y_rec ~ ", gsub(" ", "",smooth_terms))
model<-gam(as.formula(formula_str), data=sb_dat)
full_edfs<-summary(model)$s.table[, "edf"]
nyears<-length(sb_dat$year)

ts_smooths <- data.frame()
mohns <- data.frame()
results_sb <- FunctionalRelationships(sb_dat, 15,formula_str,1)

mohns_sb<- results_sb[[2]]%>%
  rename(variable='row.names(xyy)')%>%
  mutate(Species="Sablefish")

plot_dat_sb<-sb_dat%>%select(year,Y_rec,unique(results_sb[[1]]$var))%>%
  pivot_longer(cols=c(unique(results_sb[[1]]$var)))%>%
  rename(var=name)%>%
  mutate(Species="Sablefish")

##### Petrale Sole #####
modtemp <- arrange(ps_lfo5,RMSE)[1,]%>%#select best model
  select(var1, var2, var3, var4) #select variables. Way to do this with all 4 and filter out 4th when it isn't used?
numvar<-sum(!is.na(modtemp))
mod<- modtemp [1:numvar]
smooth_terms <- paste("s(",mod, ", k = 3)", collapse = " + ")
formula_str <- paste("Y_rec ~ ", gsub(" ", "",smooth_terms))
model<-gam(as.formula(formula_str), data=ps_dat)
full_edfs<-summary(model)$s.table[, "edf"]
nyears<-length(ps_dat$year)

ts_smooths <- data.frame()
mohns <- data.frame()
results_ps<- FunctionalRelationships(ps_dat, 15,formula_str,1)

mohns_ps<- results_ps[[2]]%>%
  rename(variable='row.names(xyy)')%>%
  mutate(Species="Petrale Sole")

plot_dat_ps<-ps_dat%>%select(year,Y_rec,unique(results_ps[[1]]$var))%>%
  pivot_longer(cols=c(unique(results_ps[[1]]$var)))%>%
  rename(var=name)%>%
  mutate(Species="Petrale Sole")

##### Best Models Plot #####

mohns_all <- mohns_ps%>%
  bind_rows(mohns_yt)%>%
  bind_rows(mohns_sb)

plotdat_all <- plot_dat_yt%>%
  bind_rows(plot_dat_ps)%>%
  bind_rows(plot_dat_sb)

mohns_yt_plot<-ggplot(data=mohns_yt,aes(x=termyr, y=mohns))+ #you could add observations onto this to show contrast
  #facet_wrap(~Species, nrow=3)+
  geom_point(aes(col=variable))+
  ylab("Mohn's Rho")+
  xlab("Terminal Year")+
  ylim(c(-0.1, 0.8))+
 # ylim(c(-0.5, 0.2))+
  ggtitle("Yellowtail")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

mohns_sb_plot<-ggplot(data=mohns_sb,aes(x=termyr, y=mohns))+ #you could add observations onto this to show contrast
  geom_point(aes(col=variable))+
  ylab("Mohn's Rho")+
  xlab("Terminal Year")+
  ylim(c(-0.1, 0.8))+
  #ylim(c(-0.5, 0.2))+
  ggtitle("Sablefish")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

mohns_ps_plot<-ggplot(data=mohns_ps,aes(x=termyr, y=mohns))+ #you could add observations onto this to show contrast
  geom_point(aes(col=variable))+
  ylab("Mohn's Rho")+
  xlab("Terminal Year")+
  ylim(c(-0.1, 0.8))+
  ggtitle("Petrale Sole")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))


yt_var<- ggplot()+ #you could add observations onto this to show contrast
  facet_wrap(~var)+
  geom_point(data=plot_dat_yt,alpha=0.6, aes( label=year,y=Y_rec, x=value, col=year))+
  #geom_text(data=plot_dat_yt, aes( label=year,y=Y_rec, x=value,col=year),alpha=0.6,nudge_y = 0.05)+
  xlim(c(-3,3))+
  ylab("Recruitment Deviations")+
  xlab("Standardized Value")+
  geom_line(data=data.frame(results_yt[[1]]),aes(x=x, y=fit,group=termyr,col=termyr))+
  theme_classic()

ps_var<-ggplot()+ #you could add observations onto this to show contrast
  facet_wrap(~var)+
  geom_point(data=plot_dat_ps,alpha=0.6, aes( label=year,y=Y_rec, x=value, col=year))+
  #geom_text(data=plot_dat_ps, aes( label=year,y=Y_rec, x=value,col=year),alpha=0.6,nudge_y = 0.05)+
  xlim(c(-3,3))+
  ylab("Recruitment Deviations")+
  xlab("Standardized Value")+
  geom_line(data=data.frame(results_ps[[1]]),aes(x=x, y=fit,group=termyr,col=termyr))+
  theme_classic()

sb_var<-ggplot()+ #you could add observations onto this to show contrast
  facet_wrap(~var)+
  geom_point(data=plot_dat_sb,alpha=0.6, aes( label=year,y=Y_rec, x=value, col=year))+
  #geom_text(data=plot_dat_sb, aes( label=year,y=Y_rec, x=value,col=year),alpha=0.6,nudge_y = 0.05)+
  xlim(c(-3,3))+
  xlab("Standardized Value")+
  ylab("Recruitment Deviations")+
  geom_line(data=data.frame(results_sb[[1]]),aes(x=x, y=fit,group=termyr,col=termyr))+
  theme_classic()

yt_best <- ggarrange(mohns_yt_plot,yt_var, ncol = 2, nrow = 1)
sb_best <- ggarrange(mohns_sb_plot,sb_var, ncol = 2, nrow = 1)
ps_best <- ggarrange(mohns_ps_plot,ps_var, ncol = 2, nrow = 1)

all_best_lfo5 <- ggarrange(yt_best, sb_best, ps_best, ncol = 1, nrow = 3)

pdf(file = "Output/Figures/all_best_lfo5.pdf", width = 11, height = 8)
all_best_lfo5
dev.off()

png(file = "Output/Figures/all_best_lfo5.png",width = 1100, height = 800, res = 100)
all_best_lfo5
dev.off()


#### Single Variables ####

##### Yellowtail #####
yt_single <- readRDS("Output/Data/yt_model_single_fits.rds")
yt_single_loo<-data.frame(yt_single[["LOO"]][["results"]])
res<-data.frame()
mohn<-data.frame()
results_single_yt<- single_covs(yt_single_loo,yt_dat)


mohns_yt<- results_single_yt[[2]]%>%
  mutate(Species="Yellowtail")
npeel<-10
rec_yt<-single_covs(yt_single_loo,yt_dat)
plot.dat<-yt_dat%>%select(year,Y_rec,unique(rec_yt[[1]]$var))%>%
  pivot_longer(cols=c(unique(rec_yt[[1]]$var)))%>%
  rename(var=name)

mohns10yt<- mohns_yt%>%
  group_by(variable)%>%
  summarize(mohns=sum(mohns)/length(unique(termyr)))%>%
  mutate(Species="Yellowtail")


##### Sablefish #####

sb_single <- readRDS("Output/Data/sb_model_single_fits.rds")
sb_single_loo<-data.frame(sb_single[["LOO"]][["results"]])
res<-data.frame()
mohn<-data.frame()
npeel<-10
results_single_sb<- single_covs(sb_single_loo,sb_dat)
mohns_sb<- results_single_sb[[2]]%>%
  mutate(Species="Sablefish")

plot.dat<-sb_dat%>%select(year,Y_rec,unique(rec_sb[[1]]$var))%>%
  pivot_longer(cols=c(unique(rec_sb[[1]]$var)))%>%
  rename(var=name)

mohns10sb<- mohns_sb%>%
  group_by(variable)%>%
  summarize(mohns=sum(mohns)/length(unique(termyr)))%>%
  mutate(Species="Sablefish")

##### Petrale Sole #####

ps_single <- readRDS("Output/Data/ps_model_single_fits.rds")
ps_single_loo<-data.frame(ps_single[["LOO"]][["results"]])
res<-data.frame()
mohn<-data.frame()
npeel<-10
results_single_ps<- single_covs(ps_single_loo,ps_dat)
mohns_ps<- results_single_ps[[2]]%>%
  mutate(Species="Petrale Sole")

plot.dat<-ps_dat%>%select(year,Y_rec,unique(results_single_ps[[1]]$var))%>%
  pivot_longer(cols=c(unique(results_single_ps[[1]]$var)))%>%
  rename(var=name)

mohns10ps<- mohns_ps%>%
  group_by(variable)%>%
  summarize(mohns=sum(mohns)/length(unique(termyr)))%>%
  mutate(Species="Petrale Sole")

##### Figures #####

mohns10<- mohns10ps%>%
  add_row(mohns10sb)%>%
  add_row(mohns10yt)
  
mohns10yt_plot<-ggplot(mohns10yt, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y =  variable,
  x =   mohns,
  fill=mohns
)) +
  # Note that we swapped x and y in aes() because of coord_flip()
  # facet_grid(species ~ RMSE, scales = "free_y") + # Use free_y scale
  geom_bar(stat = "identity", fill="Darkgreen") +
  # Add the scale_y_reordered function to clean up the axis labels
  scale_y_reordered() +
  #coord_flip() +
  ggtitle("Yellowtail")+
  xlim(c(-0.3,0.3))+
  labs(x = "Mohn's Rho", y = "Predictor")+
#  scale_fill_gradient2(low ="Darkgreen", mid="white", high = "Darkgreen") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

mohns10sb_plot<-ggplot(mohns10sb, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y =  variable,
  x =   mohns,
  fill=mohns
)) +
  ggtitle("Sablefish")+
  # Note that we swapped x and y in aes() because of coord_flip()
  # facet_grid(species ~ RMSE, scales = "free_y") + # Use free_y scale
  geom_bar(stat = "identity", fill="Darkgreen") +
  # Add the scale_y_reordered function to clean up the axis labels
  scale_y_reordered() +
  #coord_flip() +
  xlim(c(-0.3,0.3))+
  labs(x = "Mohn's Rho", y = "Predictor")+
  #  scale_fill_gradient2(low ="Darkgreen", mid="white", high = "Darkgreen") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))


mohns10ps_plot<-ggplot(mohns10ps, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y =  variable,
  x =   mohns,
  fill=mohns
)) +
  ggtitle("Petrale Sole")+
  # Note that we swapped x and y in aes() because of coord_flip()
  # facet_grid(species ~ RMSE, scales = "free_y") + # Use free_y scale
  geom_bar(stat = "identity", fill="Darkgreen") +
  # Add the scale_y_reordered function to clean up the axis labels
  scale_y_reordered() +
  #coord_flip() +
  xlim(c(-0.3,0.3))+
  labs(x = "Mohn's Rho", y = "Predictor")+
  #  scale_fill_gradient2(low ="Darkgreen", mid="white", high = "Darkgreen") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))



mohns_single_plot <- ggarrange(mohns10yt_plot, mohns10sb_plot, mohns10ps_plot, ncol = 3, nrow = 1)

pdf(file = "Output/Figures/mohns_single_plot.pdf", width = 8, height = 4)
mohns_single_plot
dev.off()

png(file = "Output/Figures/mohns_single_plot.png",width = 800, height = 400, res = 100)
mohns_single_plot
dev.off()
