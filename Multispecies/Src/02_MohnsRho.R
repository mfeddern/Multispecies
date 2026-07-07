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

yrfirst<- 1993
yrlast<- 2018 #2010
npeel <- 15 #8 
#### Reading in the data ####
yt_dat <- data.frame(read.csv("Data/Yellowtail/yt_fulldataset_STANDARDIZED.csv"))%>%
  select(-c(X, HCI2pjuv,HCI1pjuv, LUSIannual, HCI2larv,HCI1larv))%>%
  filter(type=="Main")%>%
  #select(-c(data.frame(unstable%>%filter(Species == "Yellowtail"))$variable))%>% #turn off when using all variables
  filter(Datatreatment=="2025 Final"&year>yrfirst&year<=yrlast)
#yt <- readRDS("Output/Data/yt_model_fits.rds")
yt <- readRDS("Output/Data/Manuscript/yt_model_fits.rds")
yt_loo<-data.frame(yt[["LOO"]][["results"]])
yt_lfo5<-data.frame(yt[["LFO5"]][["results"]])
yt_lfo10<-data.frame(yt[["LFO10"]][["results"]])

sb_dat <- data.frame(read.csv("Data/Sablefish/data-combined-glorys-sablefish_STANDARDIZED.csv"))%>%
  select(-c(X))%>%
  filter(type=="Main_RecrDev"&year>yrfirst&year<=yrlast)
#sb <- readRDS("Output/Data/sb_model_fits.rds")
sb <- readRDS("Output/Data/Manuscript/sb_model_fits.rds")
sb_loo<-data.frame(sb[["LOO"]][["results"]])
sb_lfo5<-data.frame(sb[["LFO5"]][["results"]])
sb_lfo10<-data.frame(sb[["LFO10"]][["results"]])

ps_dat <- data.frame(read.csv("Data/Petrale/DATA_Combined_glorys_petrale_STANDARDIZED.csv"))%>%
  select(-c(X, year.1))%>%
  filter(type=="Main_RecrDev"&year>yrfirst&year<=yrlast)
#ps <- readRDS("Output/Data/ps_model_fits.rds")
ps <- readRDS("Output/Data/Manuscript/ps_model_fits.rds")
ps_loo<-data.frame(ps[["LOO"]][["results"]])
ps_lfo5<-data.frame(ps[["LFO5"]][["results"]])
ps_lfo10<-data.frame(ps[["LFO10"]][["results"]])

subsethk<-readRDS("Output/Data/AnalysisPart1/hakesubset.rds")
hk_dat <- data.frame(read.csv("Data/Hake/DATA_Combined_glorys_hake_STANDARDIZED.csv"))%>%
  select(-X)%>%
  filter(type=="Main_RecrDev"&year>yrfirst&year<=yrlast)%>%
  select(all_of(subsethk),LSTyolk, Y_rec, year, type, sd)
#hk <- readRDS("Output/Data/hk_model_fits.rds")
hk <- readRDS("Output/Data/Manuscript/hk_model_fits.rds")
hk_loo<-data.frame(hk[["LOO"]][["results"]])
hk_lfo5<-data.frame(hk[["LFO5"]][["results"]])
hk_lfo10<-data.frame(hk[["LFO10"]][["results"]])
### Functions ###
par(mfrow = c(2, 2))
nyears<-length(hk_dat$Y_rec)
FunctionalRelationships<- function(data,peel,form,numvar){
  model<-gam(as.formula(form), data= data)
  full_edfs<-summary(model)$s.table[, "edf"]
  mohns<-data.frame()
  ts_smooths<-data.frame()
  #nyears<-length(data$Y_rec)
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


#### Best Models LOO ####
##### Yellowtail #####
#modtemp <- arrange(yt_lfo5,S)[1,]%>%#select best model
#  select(var1, var2, var3, var4)
#modtemp <- arrange(yt_loo,AIC)[1,]%>%#select best model
#  select(var1, var2, var3, var4) #select variables. Way to do this with all 4 and filter out 4th when it isn't used?
modtemp <- arrange(yt_lfo5,CV)[1,]%>%#select best model
  select(var1, var2, var3, var4)
numvar<-sum(!is.na(modtemp))
mod<- modtemp [1:numvar]
smooth_terms <- paste("s(",mod, ", k = 3)", collapse = " + ")
formula_str <- paste("Y_rec ~ ", gsub(" ", "",smooth_terms))
ts_smooths <- data.frame()
mohns <- data.frame()
results_yt <- FunctionalRelationships(yt_dat, npeel,formula_str,numvar)
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

formulas<-list()
res<-data.frame()
mohn<-data.frame()
#rec<-single_covs(sb_loo,sb_dat)
#modtemp <- arrange(sb_loo,AIC)[1,]%>%#select best model
#  select(var1, var2, var3, var4) #select variables. Way to do this with all 4 and filter out 4th when it isn't used?
modtemp <- arrange(sb_lfo5,CV)[3,]%>%#select best model
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
results_sb <- FunctionalRelationships(sb_dat, npeel,formula_str,numvar)

mohns_sb<- results_sb[[2]]%>%
  rename(variable='row.names(xyy)')%>%
  mutate(Species="Sablefish")

plot_dat_sb<-sb_dat%>%select(year,Y_rec,unique(results_sb[[1]]$var))%>%
  pivot_longer(cols=c(unique(results_sb[[1]]$var)))%>%
  rename(var=name)%>%
  mutate(Species="Sablefish")

##### Petrale Sole #####
modtemp <- arrange(ps_lfo5,CV)[1,]%>%#select best model
  select(var1, var2, var3, var4) #select variables. Way to do this with all 4 and filter out 4th when it isn't used?
#modtemp <- arrange(ps_lfo10,S)[1,]%>%#select best model
#  select(var1, var2, var3, var4)
numvar<-sum(!is.na(modtemp))
mod<- modtemp [1:numvar]
smooth_terms <- paste("s(",mod, ", k = 3)", collapse = " + ")
formula_str <- paste("Y_rec ~ ", gsub(" ", "",smooth_terms))
model<-gam(as.formula(formula_str), data=ps_dat)
full_edfs<-summary(model)$s.table[, "edf"]
nyears<-length(ps_dat$year)

ts_smooths <- data.frame()
mohns <- data.frame()
results_ps<- FunctionalRelationships(ps_dat, npeel,formula_str,numvar)

mohns_ps<- results_ps[[2]]%>%
  rename(variable='row.names(xyy)')%>%
  mutate(Species="Petrale Sole")

plot_dat_ps<-ps_dat%>%select(year,Y_rec,unique(results_ps[[1]]$var))%>%
  pivot_longer(cols=c(unique(results_ps[[1]]$var)))%>%
  rename(var=name)%>%
  mutate(Species="Petrale Sole")

##### Hake #####
modtemp <- arrange(hk_lfo5,CV)[1,]%>%#select best model
  select(var1, var2, var3, var4) #select variables. Way to do this with all 4 and filter out 4th when it isn't used?
#modtemp <- arrange(hk_lfo10,S)[1,]%>%#select best model
#  select(var1, var2, var3, var4) 
numvar<-sum(!is.na(modtemp))
mod<- modtemp [1:numvar]
smooth_terms <- paste("s(",mod, ", k = 3)", collapse = " + ")
formula_str <- paste("Y_rec ~ ", gsub(" ", "",smooth_terms))
model<-gam(as.formula(formula_str), data=hk_dat)
full_edfs<-summary(model)$s.table[, "edf"]
nyears<-length(hk_dat$year)

ts_smooths <- data.frame()
mohns <- data.frame()
results_hk<- FunctionalRelationships(hk_dat, npeel,formula_str,numvar)

mohns_hk<- results_hk[[2]]%>%
  rename(variable='row.names(xyy)')%>%
  mutate(Species="Hake")

plot_dat_hk<-hk_dat%>%select(year,Y_rec,unique(results_hk[[1]]$var))%>%
  pivot_longer(cols=c(unique(results_hk[[1]]$var)))%>%
  rename(var=name)%>%
  mutate(Species="Hake")



##### Best Models Plot #####

mohns_all <- mohns_ps%>%
  bind_rows(mohns_yt)%>%
  bind_rows(mohns_sb)%>%
  bind_rows(mohns_hk)

plotdat_all <- plot_dat_yt%>%
  bind_rows(plot_dat_ps)%>%
  bind_rows(plot_dat_sb)%>%
  bind_rows(plot_dat_hk)

mohns_yt_plot<-ggplot(data=mohns_yt,aes(x=termyr, y=mohns))+ #you could add observations onto this to show contrast
  #facet_wrap(~Species, nrow=3)+
  geom_point(aes(col=variable))+
  ylab("Mohn's Rho")+
  xlab("Terminal Year")+
  #  ylim(c(-0.5, 0.8))+
 ylim(c(-0.5, 0.2))+
  ggtitle("Yellowtail")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

mohns_sb_plot<-ggplot(data=mohns_sb,aes(x=termyr, y=mohns))+ #you could add observations onto this to show contrast
  geom_point(aes(col=variable))+
  ylab("Mohn's Rho")+
  xlab("Terminal Year")+
  #ylim(c(-0.5, 0.8))+
  ylim(c(-0.5, 0.2))+
  ggtitle("Sablefish")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

mohns_ps_plot<-ggplot(data=mohns_ps,aes(x=termyr, y=mohns))+ #you could add observations onto this to show contrast
  geom_point(aes(col=variable))+
  ylab("Mohn's Rho")+
  xlab("Terminal Year")+
  ylim(c(-0.5, 0.2))+
  ggtitle("Petrale Sole")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

mohns_hk_plot<-ggplot(data=mohns_hk,aes(x=termyr, y=-mohns))+ #you could add observations onto this to show contrast
  geom_point(aes(col=variable))+
  ylab("Mohn's Rho")+
  xlab("Terminal Year")+
  ylim(c(-0.5, 0.2))+
  ggtitle("Hake")+
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

hk_var<-ggplot()+ #you could add observations onto this to show contrast
  facet_wrap(~var)+
  geom_point(data=plot_dat_hk,alpha=0.6, aes( label=year,y=Y_rec, x=value, col=year))+
  #geom_text(data=plot_dat_hk, aes( label=year,y=Y_rec, x=value,col=year),alpha=0.6,nudge_y = 0.05)+
  xlim(c(-3,3))+
  xlab("Standardized Value")+
  ylab("Recruitment Deviations")+
  geom_line(data=data.frame(results_hk[[1]]),aes(x=x, y=fit,group=termyr,col=termyr))+
  theme_classic()

yt_best <- ggarrange(mohns_yt_plot,yt_var, ncol = 2, nrow = 1)
sb_best <- ggarrange(mohns_sb_plot,sb_var, ncol = 2, nrow = 1)
ps_best <- ggarrange(mohns_ps_plot,ps_var, ncol = 2, nrow = 1)
hk_best <- ggarrange(mohns_hk_plot,hk_var, ncol = 2, nrow = 1)

all_best_loo <- ggarrange(yt_best, sb_best, ps_best,hk_best, ncol = 1, nrow = 4)
all_best_loo
pdf(file = "Output/Figures/Manuscript/all_best_lfo5.pdf", width = 11, height = 8)
all_best_loo
dev.off()

png(file = "Output/Figures/Manuscript/all_best_lf05.png",width = 1100, height = 800, res = 100)
all_best_loo
dev.off()


#### Single Variables ####

##### Yellowtail #####
yt_single <- readRDS("Output/Data/Manuscript/yt_model_single_fits.rds")
yt_single_loo<-data.frame(yt_single[["LOO"]][["results"]])
res<-data.frame()
mohn<-data.frame()
results_single_yt<- single_covs(yt_single_loo,yt_dat)


mohns_yt<- results_single_yt[[2]]%>%
  mutate(Species="Yellowtail")

rec_yt<-single_covs(yt_single_loo,yt_dat)
plot.dat<-yt_dat%>%select(year,Y_rec,unique(rec_yt[[1]]$var))%>%
  pivot_longer(cols=c(unique(rec_yt[[1]]$var)))%>%
  rename(var=name)

mohns10yt<- mohns_yt%>%
  group_by(variable)%>%
  summarize(mohns=sum(mohns)/length(unique(termyr)))%>%
  mutate(Species="Yellowtail Rockfish")%>%
  mutate(category = case_when(
    str_detect(variable, "CutiSTIpjuv")  ~ "Yes",
    str_detect(variable, "ONIpjuv") ~ "Yes",
    str_detect(variable, "LSTpjuv") ~ "Yes",
    TRUE                              ~ "No"       # Default catch-all
  ))


##### Sablefish #####

sb_single <- readRDS("Output/Data/Manuscript/sb_model_single_fits.rds")
sb_single_loo<-data.frame(sb_single[["LOO"]][["results"]])
res<-data.frame()
mohn<-data.frame()
results_single_sb<- single_covs(sb_single_loo,sb_dat)
mohns_sb<- results_single_sb[[2]]%>%
  mutate(Species="Sablefish")

plot.dat<-sb_dat%>%select(year,Y_rec,unique(results_single_sb[[1]]$var))%>%
  pivot_longer(cols=c(unique(results_single_sb[[1]]$var)))%>%
  rename(var=name)

mohns10sb<- mohns_sb%>%
  group_by(variable)%>%
  summarize(mohns=sum(mohns)/length(unique(termyr)))%>%
  mutate(Species="Sablefish")%>%
  mutate(category = case_when(
    str_detect(variable, "CSTegg")  ~ "Yes",
    str_detect(variable, "DDpre") ~ "Yes",
    str_detect(variable, "DDegg") ~ "Yes",
    str_detect(variable, "LSTyolk") ~ "Yes",
    str_detect(variable, "DDlarv") ~ "Yes",
    TRUE                              ~ "No"       # Default catch-all
  ))

##### Hake #####

hk_single <- readRDS("Output/Data/Manuscript/hk_model_single_fits.rds")
hk_single_loo<-data.frame(hk_single[["LOO"]][["results"]])
res<-data.frame()
mohn<-data.frame()
results_single_hk<- single_covs(hk_single_loo,hk_dat)
mohns_hk<- results_single_hk[[2]]%>%
  mutate(Species="Hake")

plot.dat<-hk_dat%>%select(year,Y_rec,unique(results_single_hk[[1]]$var))%>%
  pivot_longer(cols=c(unique(results_single_hk[[1]]$var)))%>%
  rename(var=name)

mohns10hk<- mohns_hk%>%
  group_by(variable)%>%
  summarize(mohns=sum(mohns)/length(unique(termyr)))%>%
  mutate(Species="Hake")%>%
  mutate(category = case_when(
    str_detect(variable, "LSTyolk")  ~ "Yes",
    TRUE                              ~ "No"       # Default catch-all
  ))

##### Petrale Sole #####

ps_single <- readRDS("Output/Data/Manuscript/ps_model_single_fits.rds")
ps_single_loo<-data.frame(ps_single[["LOO"]][["results"]])
res<-data.frame()
mohn<-data.frame()
#npeel<-15
results_single_ps<- single_covs(ps_single_loo,ps_dat)
mohns_ps<- results_single_ps[[2]]%>%
  mutate(Species="Petrale Sole")

plot.dat<-ps_dat%>%select(year,Y_rec,unique(results_single_ps[[1]]$var))%>%
  pivot_longer(cols=c(unique(results_single_ps[[1]]$var)))%>%
  rename(var=name)

mohns10ps<- mohns_ps%>%
  group_by(variable)%>%
  summarize(mohns=sum(mohns)/length(unique(termyr)))%>%
  mutate(Species="Petrale Sole")%>%
  mutate(category = case_when(
    str_detect(variable, "DDpre")  ~ "Yes",
    str_detect(variable, "MLDegg") ~ "Yes",
    str_detect(variable, "CSTlarv") ~ "Yes",
    str_detect(variable, "CSTbjuv") ~ "Yes",
    str_detect(variable, "DDpjuv") ~ "Yes",
    str_detect(variable, "LSTlarv") ~ "Yes",
    TRUE                              ~ "No"       # Default catch-all
  ))

##### Figures #####

mohns10<- mohns10ps%>%
  add_row(mohns10sb)%>%
  add_row(mohns10yt)%>%
  add_row(mohns10hk)
write_rds(mohns10, "Output/Data/Manuscript/mohns.rds")


mohns10yt_plot<-ggplot(mohns10yt, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y =  reorder( variable,abs(mohns), .desc = TRUE),
  x =   abs(mohns),
  fill=mohns
)) +
  # Note that we swapped x and y in aes() because of coord_flip()
  # facet_grid(species ~ RMSE, scales = "free_y") + # Use free_y scale
  geom_bar(stat = "identity", fill="Darkgreen") +
  # Add the scale_y_reordered function to clean up the axis labels
  scale_y_reordered() +
  #coord_flip() +
  ggtitle("Yellowtail Rockfish")+
  xlim(c(-0,0.7))+
  labs(x = "", y = "")+
  geom_text(
    aes(label = ifelse(category == "Yes", "*", "")), 
    nudge_y = -.25,x=-0.015,  # Nudges the star slightly above the top of the bar
    size = 5        # Adjusts the size of the star
  ) +
  geom_vline(xintercept=0.15, lty=2, col="grey")+
  
#  scale_fill_gradient2(low ="Darkgreen", mid="white", high = "Darkgreen") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
mohns10yt_plot

mohns10sb_plot<-ggplot(mohns10sb, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y =  reorder( variable,abs(mohns), .desc = TRUE),
  x =   abs(mohns),
  fill=mohns
)) +
  ggtitle("Sablefish")+
  # Note that we swapped x and y in aes() because of coord_flip()
  # facet_grid(species ~ RMSE, scales = "free_y") + # Use free_y scale
  geom_bar(stat = "identity", fill="Darkgreen") +
  # Add the scale_y_reordered function to clean up the axis labels
  scale_y_reordered() +
  #coord_flip() +
  xlim(c(-0,0.7))+
  labs(x = "", y = "")+
  geom_text(
    aes(label = ifelse(category == "Yes", "*", "")), 
    nudge_y = -.25,x=-0.015,  # Nudges the star slightly above the top of the bar
    size = 5        # Adjusts the size of the star
  ) +
  geom_vline(xintercept=0.15, lty=2, col="grey")+
  
  #  scale_fill_gradient2(low ="Darkgreen", mid="white", high = "Darkgreen") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
mohns10sb_plot

mohns10ps_plot<-ggplot(mohns10ps, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y =  reorder( variable,abs(mohns), .desc = TRUE),
  x =  abs(mohns),
  fill=mohns)) +
  ggtitle("Petrale Sole")+
  # Note that we swapped x and y in aes() because of coord_flip()
  # facet_grid(species ~ RMSE, scales = "free_y") + # Use free_y scale
  geom_bar(stat = "identity", fill="Darkgreen") +
  # Add the scale_y_reordered function to clean up the axis labels
  scale_y_reordered() +
  #coord_flip() +
  xlim(c(-0,0.7))+
  geom_text(
    aes(label = ifelse(category == "Yes", "*", "")), 
    nudge_y = -.25,x=-0.015,  # Nudges the star slightly above the top of the bar
    size = 5        # Adjusts the size of the star
  ) +
  geom_vline(xintercept=0.15, lty=2, col="grey")+
  
  labs(x = "", y = "")+
  #  scale_fill_gradient2(low ="Darkgreen", mid="white", high = "Darkgreen") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
mohns10ps_plot

mohns10hk_plot<-ggplot(mohns10hk, aes(
  # Use reorder_within, specifying 'cov', 'total_rmse', and the grouping variable 'species'
  y =  reorder( variable,abs(mohns), .desc = TRUE),
  x =  abs(mohns),
  fill=mohns)) +
  ggtitle("Hake")+
  # Note that we swapped x and y in aes() because of coord_flip()
  # facet_grid(species ~ RMSE, scales = "free_y") + # Use free_y scale
  geom_bar(stat = "identity", fill="Darkgreen") +
  # Add the scale_y_reordered function to clean up the axis labels
  scale_y_reordered() +
  #coord_flip() +
  xlim(c(-0,0.7))+
  geom_text(
    aes(label = ifelse(category == "Yes", "*", "")), 
    nudge_y = -.25,x=-0.015,  # Nudges the star slightly above the top of the bar
    size = 5        # Adjusts the size of the star
  ) +
  geom_vline(xintercept=0.15, lty=2, col="grey")+
  labs(x = "", y = "")+
  #  scale_fill_gradient2(low ="Darkgreen", mid="white", high = "Darkgreen") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
mohns10hk_plot

mohns_single_plot <- ggarrange( mohns10hk_plot,mohns10ps_plot, mohns10yt_plot, mohns10sb_plot,ncol = 2, nrow = 2, labels=c("A", "B", "C", "D"))
mohns_single_plot <-annotate_figure(mohns_single_plot, left= "Predictor",
                bottom = text_grob(expression(paste("|Mohn's", ~rho,"|"))))
ggsave("Output/Figures/Manuscript/mohns_single_plot.pdf", plot =mohns_single_plot, width = 8, height = 6.5)
ggsave("Output/Figures/Manuscript/mohns_single_plot.png", plot = mohns_single_plot, width = 8, height = 6.5,bg = "white", dpi = 300)



#### RMSE Comparison ####
results_single<-results_single_sb[[2]]%>%mutate(Species="Sablefish")%>%
  bind_rows(results_single_ps[[2]]%>%mutate(Species="Petrale Sole"))%>%
  bind_rows(results_single_yt[[2]]%>%mutate(Species="Yellowtail"))%>%
  bind_rows(results_single_hk[[2]]%>%mutate(Species="Hake"))
results_single_summ<-results_single%>%group_by(termyr, Species)%>%
  summarise(mean=mean(-abs(mohns)), sd=sd(-abs(mohns)))
ggplot(data=results_single, aes(x=termyr, y=-abs(mohns)))+
 # geom_point(se = FALSE, aes(col=variable), col = 'grey')+
  geom_smooth(aes(group=Species, col=Species))+
  theme_classic()
arrange(results_single%>%
  group_by(Species,variable)%>%
  summarise(mean=mean(mohns)), mean)

mohnstrend<- ggplot(data=results_single, aes(x=termyr, y=-abs(mohns)))+
  geom_point(data=results_single_summ,aes(y=mean, x=termyr, col=Species))+
  geom_smooth(aes(group=Species, col=Species))+#, method='gam')+
  ylab("Mean Mohns Rho")+
  xlab("Terminal Year")+
  theme_classic()

ggplot(data=results_single%>%filter(Species=='Yellowtail'), aes(x=termyr, y=-abs(mohns)))+
  geom_point(aes(col=variable))+
  geom_smooth(aes(col=variable),se = FALSE)+#, method='gam')+
  ylab("Mean Mohns Rho")+
  xlab("Terminal Year")+
  theme_classic()

pdf(file = "Output/Figures/Manuscript/mohnstrend.pdf", width = 5, height = 5)
mohnstrend
dev.off()

png(file = "Output/Figures/Manuscript/mohnstrend.png",width = 1500, height = 1500, res = 300)
mohnstrend
dev.off()


plot_dat_yt=yt_dat%>%select(year,Y_rec,unique(results_single_yt[[1]]$var))%>%
  pivot_longer(cols=c(unique(results_single_yt[[1]]$var)))%>%
  rename(var=name)%>%
  mutate(Species="Yellowtail")

 ggplot()+ #you could add observations onto this to show contrast
  facet_wrap(~var)+
  geom_point(data=plot_dat_yt,alpha=0.6, aes( label=year,y=Y_rec, x=value, col=year))+
  #geom_text(data=plot_dat_yt, aes( label=year,y=Y_rec, x=value,col=year),alpha=0.6,nudge_y = 0.05)+
  xlim(c(-3,3))+
  ylab("Recruitment Deviations")+
  xlab("Standardized Value")+
  geom_line(data=data.frame(results_yt[[1]]),aes(x=x, y=fit,group=termyr,col=termyr))+
  theme_classic()

 
##### small plot

 
 plot_dat_sb<-sb_dat%>%select(year,Y_rec,unique(results_single_sb[[1]]$var))%>%
   pivot_longer(cols=c(unique(results_single_sb[[1]]$var)))%>%
   rename(var=name)%>%
   mutate(Species="Sablefish")%>%
   left_join(mohns10sb%>%rename(var=variable, rho=mohns)%>%select(var,rho))
 
 ggplot()+ #you could add observations onto this to show contrast
   facet_wrap(~var)+
   geom_point(data = plot_dat_sb, alpha = 1, shape = 21, size = 2,
              aes(label = year, y = Y_rec, x = value, group = year, fill = year), col="white") +
   #geom_text(data=plot_dat_yt, aes( label=year,y=Y_rec, x=value,col=year),alpha=0.6,nudge_y = 0.05)+
   xlim(c(-3,3))+
   ylim(c(-2,3))+
    ggtitle("Sablefish")+
   ylab("ln(Recruitment Deviations)")+
   xlab("Standardized Predictor Value")+
   geom_text(data=plot_dat_sb,aes(label=round(rho,2), x=-2, y=2.8))+
   geom_line(data=data.frame(results_single_sb[[1]]),aes(x=x, y=fit,group=termyr,col=termyr))+
   theme_classic()+
   scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Point Year") +
   scale_color_gradient(low = "pink", high = "darkred", name = "Line Term Year") +
   theme(plot.title = element_text(hjust = 0.5))
   
plot_dat_sb<-sb_dat%>%select(year,Y_rec,unique(results_single_sb[[1]]$var))%>%
   pivot_longer(cols=c(unique(results_single_sb[[1]]$var)))%>%
   rename(var=name)%>%
   mutate(Species="Sablefish")%>%
   left_join(mohns10sb%>%rename(var=variable, rho=mohns)%>%select(var,rho))
 
functional_sb<-ggplot()+ #you could add observations onto this to show contrast
   facet_wrap(~var)+
   geom_point(data = plot_dat_sb, alpha = 1, shape = 21, size = 2,
              aes(label = year, y = Y_rec, x = value, group = year, fill = year), col="white") +
   #geom_text(data=plot_dat_yt, aes( label=year,y=Y_rec, x=value,col=year),alpha=0.6,nudge_y = 0.05)+
   xlim(c(-2,2))+
   ggtitle("Sablefish")+
   ylab("ln(Recruitment Deviations)")+
   xlab("Standardized Predictor Value")+
   geom_text(data=plot_dat_sb,aes(label=round(rho,2), x=-1.25, y=2.8))+
   geom_line(data=data.frame(results_single_sb[[1]]),aes(x=x, y=fit,group=termyr,col=termyr))+
   theme_classic()+
   ylim(c(-2,3))+
   scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Year") +
   scale_color_gradient(low = "pink", high = "darkred", name = "Terminal Year") +
   theme(plot.title = element_text(hjust = 0.5))
   
ggsave("Output/Figures/Manuscript/functional_sb.pdf", plot = functional_sb, width = 7, height = 7)
ggsave("Output/Figures/Manuscript/functional_sb.png", plot =  functional_sb, width = 7, height = 7,bg = "white", dpi = 300)


 plot_dat_yt<-yt_dat%>%select(year,Y_rec,unique(results_single_yt[[1]]$var))%>%
   pivot_longer(cols=c(unique(results_single_yt[[1]]$var)))%>%
   rename(var=name)%>%
   mutate(Species="Sablefish")%>%
   left_join(mohns10yt%>%rename(var=variable, rho=mohns)%>%select(var,rho))
 
 functional_yt<-ggplot()+ #you could add observations onto this to show contrast
   facet_wrap(~var)+
   geom_point(data = plot_dat_yt, alpha = 1, shape = 21, size = 2,
              aes(label = year, y = Y_rec, x = value, group = year, fill = year), col="white") +
   #geom_text(data=plot_dat_yt, aes( label=year,y=Y_rec, x=value,col=year),alpha=0.6,nudge_y = 0.05)+
   xlim(c(-2,2))+
   ggtitle("Yellowtail Rockfish")+
   ylab("ln(Recruitment Deviations)")+
   xlab("Standardized Predictor Value")+
   geom_text(data=plot_dat_yt,aes(label=round(rho,2), x=-1.25, y=2))+
   geom_line(data=data.frame(results_single_yt[[1]]),aes(x=x, y=fit,group=termyr,col=termyr))+
   theme_classic()+
   ylim(c(-2,2.8))+
   scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Year") +
   scale_color_gradient(low = "pink", high = "darkred", name = "Terminal Year") +
   theme(plot.title = element_text(hjust = 0.5))
 ggsave("Output/Figures/Manuscript/functional_yt.pdf", plot = functional_yt, width = 7, height = 7)
 ggsave("Output/Figures/Manuscript/functional_yt.png", plot =  functional_yt, width = 7, height = 7,bg = "white", dpi = 300)
 
 
 plot_dat_hk<-hk_dat%>%select(year,Y_rec,unique(results_single_hk[[1]]$var))%>%
   pivot_longer(cols=c(unique(results_single_hk[[1]]$var)))%>%
   rename(var=name)%>%
   mutate(Species="Sablefish")%>%
   left_join(mohns10hk%>%rename(var=variable, rho=mohns)%>%select(var,rho))
 
 functional_hk<-ggplot()+ #you could add observations onto this to show contrast
   facet_wrap(~var)+
   geom_point(data = plot_dat_hk, alpha = 1, shape = 21, size = 2,
              aes(label = year, y = Y_rec, x = value, group = year, fill = year), col="white") +
   #geom_text(data=plot_dat_hk, aes( label=year,y=Y_rec, x=value,col=year),alpha=0.6,nudge_y = 0.05)+
   xlim(c(-2,2))+
   ggtitle("Pacific Hake")+
   ylab("ln(Recruitment Deviations)")+
   xlab("Standardized Predictor Value")+   geom_text(data=plot_dat_hk,aes(label=round(rho,2), x=-1.25, y=2))+
   geom_line(data=data.frame(results_single_hk[[1]]),aes(x=x, y=fit,group=termyr,col=termyr))+
   theme_classic()+
   ylim(c(-2,2.8))+
   scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Year") +
   scale_color_gradient(low = "pink", high = "darkred", name = "Terminal Year") +
   theme(plot.title = element_text(hjust = 0.5))
 
 ggsave("Output/Figures/Manuscript/functional_hk.pdf", plot = functional_hk, width = 7, height = 7)
 ggsave("Output/Figures/Manuscript/functional_hk.png", plot =  functional_hk, width = 7, height = 7,bg = "white", dpi = 300)
 
 
 plot_dat_ps<-ps_dat%>%select(year,Y_rec,unique(results_single_ps[[1]]$var))%>%
   pivot_longer(cols=c(unique(results_single_ps[[1]]$var)))%>%
   rename(var=name)%>%
   mutate(Species="Sablefish")%>%
   left_join(mohns10ps%>%rename(var=variable, rho=mohns)%>%select(var,rho))
 
 functional_ps<-ggplot()+ #you could add observations onto this to show contrast
   facet_wrap(~var)+
   geom_point(data = plot_dat_ps, alpha = 1, shape = 21, size = 2,
              aes(label = year, y = Y_rec, x = value, group = year, fill = year), col="white") +
   #geom_text(data=plot_dat_ps, aes( label=year,y=Y_rec, x=value,col=year),alpha=0.6,nudge_y = 0.05)+
   xlim(c(-2,2))+
   ggtitle("Petrale Sole")+
   ylab("ln(Recruitment Deviations)")+
   xlab("Standardized Predictor Value")+
   geom_text(data=plot_dat_ps,aes(label=round(rho,2), x=-1.25, y=2))+
   geom_line(data=data.frame(results_single_ps[[1]]),aes(x=x, y=fit,group=termyr,col=termyr))+
   theme_classic()+
   ylim(c(-2,2.8))+
   scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Year") +
   scale_color_gradient(low = "pink", high = "darkred", name = "Terminal Year") +
   theme(plot.title = element_text(hjust = 0.5))

 ggsave("Output/Figures/Manuscript/functional_ps.pdf", plot = functional_ps, width = 7, height = 7)
 ggsave("Output/Figures/Manuscript/functional_ps.png", plot =  functional_ps, width = 7, height = 7,bg = "white", dpi = 300)
 