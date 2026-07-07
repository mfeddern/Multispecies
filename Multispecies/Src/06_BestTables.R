#### Reading in the data ####
yt_dat <- data.frame(read.csv("Data/Yellowtail/yt_fulldataset_STANDARDIZED.csv"))%>%
  select(-c(X, HCI2pjuv,HCI1pjuv, LUSIannual, HCI2larv,HCI1larv))%>%
  filter(type=="Main")%>%
  #select(-c(data.frame(unstable%>%filter(Species == "Yellowtail"))$variable))%>% #turn off when using all variables
  filter(Datatreatment=="2025 Final"&year>yrfirst&year<=yrlast)
#yt <- readRDS("Output/Data/yt_model_fits.rds")


#### Making a GIANT table ####
colnames(yt_lfo5)
yt <- readRDS("Output/Data/JuneUpdate/yt_model_fits.rds")
yt_loo<-data.frame(yt[["LOO"]][["results"]])
yt_lfo5<-data.frame(yt[["LFO5"]][["results"]])
yt_lfo10<-data.frame(yt[["LFO10"]][["results"]])


yt_AIC_select<-arrange(yt_loo,AIC)%>%#select best model
  mutate(AICdiff=AIC-min(AIC))%>%
   filter(AICdiff<2)%>%
  mutate(selection="AIC")%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)%>%
  mutate(Rank=ifelse(numvar==min(numvar),'Best', "Competing"))

yt_loo_select<-arrange(yt_loo,CV)%>%#select best model
  mutate(Sdiff=CV-min(CV))%>%
  filter(Sdiff<S[1])%>%
  mutate(selection="LOO")%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)%>%
  mutate(Rank=ifelse(numvar==min(numvar),'Best', "Competing"))

yt_lfo5_select<-arrange(yt_lfo5,CV)%>%#select best model
  mutate(Sdiff=CV-min(CV))%>%
  filter(Sdiff<S[1])%>%
  mutate(selection="LFO-5")%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)%>%
  mutate(Rank=ifelse(numvar==min(numvar),'Best', "Competing"))

yt_lfo10_select<-arrange(yt_lfo10,CV)%>%#select best model
  mutate(Sdiff=CV-min(CV))%>%
  filter(Sdiff<S[1])%>%
  mutate(selection="LFO-10")%>%
  mutate(Rank=ifelse(numvar==min(numvar),'Best', "Competing"))%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)

compdf_yt<- yt_AIC_select%>%
  bind_rows(yt_loo_select)%>%
  bind_rows(yt_lfo5_select)%>%
  bind_rows(yt_lfo10_select)

bestdf_yt<- compdf_yt%>%
  filter(Rank=="Best")%>%
  select(species, selection,Variables, Comp)


colnames(yt_lfo5)



#### sablefish ####

sb <- readRDS("Output/Data/JuneUpdate/sb_model_fits.rds")
sb_loo<-data.frame(sb[["LOO"]][["results"]])
sb_lfo5<-data.frame(sb[["LFO5"]][["results"]])
sb_lfo10<-data.frame(sb[["LFO10"]][["results"]])


sb_AIC_select<-arrange(sb_loo,AIC)%>%#select best model
  mutate(AICdiff=AIC-min(AIC))%>%
  filter(AICdiff<2)%>%
  mutate(selection="AIC")%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)%>%
  mutate(Rank=ifelse(numvar==min(numvar),'Best', "Competing"))

sb_loo_select<-arrange(sb_loo,CV)%>%#select best model
  mutate(Sdiff=CV-min(CV))%>%
  filter(Sdiff<S[1])%>%
  mutate(selection="LOO")%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)%>%
  mutate(Rank=ifelse(numvar==min(numvar),'Best', "Competing"))

sb_lfo5_select<-arrange(sb_lfo5,CV)%>%#select best model
  mutate(Sdiff=CV-min(CV))%>%
  filter(Sdiff<S[1])%>%
  mutate(selection="LFO-5")%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)%>%
  mutate(Rank=ifelse(numvar==min(numvar),'Best', "Competing"))

sb_lfo10_select<-arrange(sb_lfo10,CV)%>%#select best model
  mutate(Sdiff=CV-min(CV))%>%
  filter(Sdiff<S[1])%>%
  mutate(selection="LFO-10")%>%
  mutate(Rank=ifelse(numvar==min(numvar),'Best', "Competing"))%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)

compdf_sb<- sb_AIC_select%>%
  bind_rows(sb_loo_select)%>%
  bind_rows(sb_lfo5_select)%>%
  bind_rows(sb_lfo10_select)

bestdf_sb<- compdf_sb%>%
  filter(Rank=="Best")%>%
  select(species, selection,Variables, Comp)



#### hake ####

hk <- readRDS("Output/Data/JuneUpdate/hk_model_fits.rds")
hk_loo<-data.frame(hk[["LOO"]][["results"]])
hk_lfo5<-data.frame(hk[["LFO5"]][["results"]])
hk_lfo10<-data.frame(hk[["LFO10"]][["results"]])


hk_AIC_select<-arrange(hk_loo,AIC)%>%#select best model
  mutate(AICdiff=AIC-min(AIC))%>%
  filter(AICdiff<2)%>%
  mutate(selection="AIC")%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)%>%
  mutate(Rank=ifelse(numvar!=min(numvar), "Competing",ifelse(CV!=min(CV),"Competing", "Best")))

hk_loo_select<-arrange(hk_loo,CV)%>%#select best model
  mutate(Sdiff=CV-min(CV))%>%
  filter(Sdiff<S[1])%>%
  mutate(selection="LOO")%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)%>%
  mutate(Rank=ifelse(CV==min(CV),"Best", ifelse(numvar!=min(numvar), "Competing","Competing")))
  
hk_lfo5_select<-arrange(hk_lfo5,CV)%>%#select best model
  mutate(Sdiff=CV-min(CV))%>%
  filter(Sdiff<S[1])%>%
  mutate(selection="LFO-5")%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)%>%
  mutate(Rank=ifelse(numvar!=min(numvar), "Competing",ifelse(CV!=min(CV),"Competing", "Best")))
  
hk_lfo10_select<-arrange(hk_lfo10,CV)%>%#select best model
  mutate(Sdiff=CV-min(CV))%>%
  filter(Sdiff<S[1])%>%
  mutate(selection="LFO-10")%>%
  mutate(Rank=ifelse(numvar!=min(numvar), "Competing",ifelse(CV==min(CV), "Best","Competing")))%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)

compdf_hk<- hk_AIC_select%>%
  bind_rows(hk_loo_select)%>%
  bind_rows(hk_lfo5_select)%>%
  bind_rows(hk_lfo10_select)

bestdf_hk<- compdf_hk%>%
  filter(Rank=="Best")%>%
  select(species, selection,Variables, Comp, CV)


dat<-hk_loo[which(hk_loo$ModelID %in% unique(compdf_hk$ModelID)), ]

dat<-hk_lfo5[which(hk_lfo5$ModelID %in% unique(compdf_hk$ModelID)), ]%>%
  left_join(compdf_hk%>%select(ModelID, selection))

ggplot(data=dat, aes(y=rsq, x=selection))+
  geom_violin()

#### petrale sole ####

ps <- readRDS("Output/Data/JuneUpdate/ps_model_fits.rds")
ps_loo<-data.frame(ps[["LOO"]][["results"]])
ps_lfo5<-data.frame(ps[["LFO5"]][["results"]])
ps_lfo10<-data.frame(ps[["LFO10"]][["results"]])


ps_AIC_select<-arrange(ps_loo,AIC)%>%#select best model
  mutate(AICdiff=AIC-min(AIC))%>%
  filter(AICdiff<2)%>%
  mutate(selection="AIC")%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)%>%
  mutate(Rank=ifelse(numvar==min(numvar),'Best', "Competing"))

ps_loo_select<-arrange(ps_loo,CV)%>%#select best model
  mutate(Sdiff=CV-min(CV))%>%
  filter(Sdiff<S[1])%>%
  mutate(selection="LOO")%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)%>%
  mutate(Rank=ifelse(numvar!=min(numvar), "Competing",ifelse(CV!=min(CV),"Competing", "Best")))
  
ps_lfo5_select<-arrange(ps_lfo5,CV)%>%#select best model
  mutate(Sdiff=CV-min(CV))%>%
  filter(Sdiff<S[1])%>%
  mutate(selection="LFO-5")%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)%>%
  mutate(Rank=ifelse(numvar!=min(numvar), "Competing",ifelse(CV!=min(CV),"Competing", "Best")))
  
ps_lfo10_select<-arrange(ps_lfo10,CV)%>%#select best model
  mutate(Sdiff=CV-min(CV))%>%
  filter(Sdiff<S[1])%>%
  mutate(selection="LFO-10")%>%
  mutate(Rank=ifelse(numvar!=min(numvar), "Competing",ifelse(CV!=min(CV),"Competing", "Best")))%>%
  mutate(Comp=length(ModelID))%>%
  unite(col = "Variables", var1, var2, var3, var4,  sep = ", ", na.rm = TRUE)

compdf_ps<- ps_AIC_select%>%
  bind_rows(ps_loo_select)%>%
  bind_rows(ps_lfo5_select)%>%
  bind_rows(ps_lfo10_select)

bestdf_ps<- compdf_ps%>%
  filter(Rank=="Best")%>%
  select(species, selection,Variables, Comp)


model_selection<-compdf_ps%>%
  bind_rows(compdf_yt)%>%
  bind_rows(compdf_hk)%>%
  bind_rows(compdf_sb)
write_rds(model_selection, "Output/Data/JuneUpdate/model_selection.rds")  
  