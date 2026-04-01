library(dplyr)
library(reshape2)
library(bayesdfa)
library(MCMCvis)
library(ggplot2)
library(stringr)
library(ggpubr)

#### Data ####
yrlast<-2018
yrfirst<-1994
marg<-readRDS("Output/Data/marginals.rds")%>%
  mutate(names=paste(species,'_', cov, sep = ""))


mohns_dat<- readRDS("Output/mohns.rds")%>%
  mutate(mohns=abs(mohns))
ps_dat <- data.frame(read.csv("Data/Petrale/DATA_Combined_glorys_petrale_STANDARDIZED.csv"))%>%
  select(-c(X, year.1))%>%
  # filter(type=="Main_RecrDev")%>%
  filter(year>yrfirst&year<=yrlast)
sb_dat <- data.frame(read.csv("Data/Sablefish/data-combined-glorys-sablefish_STANDARDIZED.csv"))%>%
  select(-c(X))%>%
  filter(type=="Main_RecrDev")%>%
  filter(year>yrfirst&year<=yrlast)
yt_dat <- data.frame(read.csv("Data/Yellowtail/yt_fulldataset_STANDARDIZED.csv"))%>%
  select(-c(X, hci2_pjuv,hci1_pjuv, lusi_annual, hci2_larv,hci1_larv))%>%
  filter(type=="Main")%>%
  filter(Datatreatment=="2025 Final"&year>yrfirst&year<=yrlast)

subsethk<-readRDS("Output/Data/hakesubset.rds")
hk_dat <- data.frame(read.csv("Data/Hake/DATA_Combined_glorys_hake_STANDARDIZED.csv"))%>%
  select(-X)%>%
  filter(type=="Main_RecrDev"&year>yrfirst&year<=yrlast)%>%
  select(all_of(subsethk), year)


#### Wrangle Data ####
hk_combined <-pivot_longer(hk_dat,-year)%>%
  mutate(Species="Hake")%>%
  rename(variable=name)%>%
  left_join(mohns_dat)%>%
  mutate(names=paste("Hake_", variable, sep = ""))

ps_combined <-pivot_longer(ps_dat%>%select(c(-type, -Y_rec, -sd)), -year)%>%
  mutate(Species="Petrale Sole")%>%
  rename(variable=name)%>%
  left_join(mohns_dat)%>%
  mutate(names=paste("Petrale Sole_", variable, sep = ""))

sb_combined<-pivot_longer(sb_dat%>%select(c(-type, -Y_rec, -sd)), -year)%>%
  mutate(Species="Sablefish")%>%
  rename(variable=name)%>%
  left_join(mohns_dat)%>%
  mutate(names=paste("Sablefish_", variable, sep = ""))


yt_combined <-pivot_longer(yt_dat%>%select(c(-type, -Y_rec, -sd,-Datatreatment)), -year)%>%
  mutate(Species="Yellowtail")%>%
  rename(variable=name)%>%
  left_join(mohns_dat)%>%
  mutate(names=paste("Yellowtail_", variable, sep = ""))

combined<-hk_combined%>%
  bind_rows(ps_combined)%>%
  bind_rows(sb_combined)%>%
  bind_rows(yt_combined)

colnames(ps_dat)[which(colnames(ps_dat)%in%ps_combined$variable)] <- unique(ps_combined$names)
colnames(yt_dat)[which(colnames(yt_dat)%in%yt_combined$variable)] <- unique(yt_combined$names)
colnames(hk_dat)[which(colnames(hk_dat)%in%hk_combined$variable)] <- unique(hk_combined$names)
colnames(sb_dat)[which(colnames(sb_dat)%in%sb_combined$variable)] <- unique(sb_combined$names)

all_dat<-ps_dat%>%select(year,unique(ps_combined$names))%>%
  left_join(hk_dat%>%select(year,unique(hk_combined$names)))%>% 
  left_join(sb_dat%>%select(year,unique(sb_combined$names)))%>%
  left_join(yt_dat%>%select(year,unique(yt_combined$names)))
#### Contrast Comparison ####
summary<-combined%>%
  mutate(period=ifelse(year>=2012, "After", "Before"))%>%
  mutate(type=ifelse(mohns>=0.1, "Unstable",ifelse(mohns<0.01, "Stable", "Moderate")))%>%
  group_by(names, period, type,Species)%>%
  summarise(mean=mean(value), sd=sd(value))%>%
  mutate(grouping=paste(period, type))

summary2<-summary%>%pivot_wider(-c(grouping,sd),names_from=period,
                      values_from = mean)%>%
  mutate(difference=Before-After)%>%
  left_join(marg%>%filter(RMSE=="LFO 5")%>%select(rmse_01 , names))%>%
  distinct()

gam<-lm(rmse_01~difference*as.factor(Species), data=summary2)
summary(gam)
ggplot(data=summary2,aes(y=rmse_01,x=difference, col=as.factor(Species)))+
  geom_point()+
  geom_smooth(method = "lm")
  

sample_size =summary%>% group_by(period) %>% summarize(num=n())

ba<-ggplot(data = summary, aes(x = period, y =mean, fill=period)) +
  #facet_grid(~Species)+
 # geom_point()+
  geom_violin() +
  geom_boxplot(width=0.1, color="black", alpha=0.2) +
  #scale_fill_viridis(discrete = TRUE) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) 


pdf(file = "Output/Figures/ba.pdf", width =9, height =9)
ba
dev.off()

png(file = "Output/Figures/ba.png",width = 900, height = 900, res = 100)
ba
dev.off()


model <- aov(mean~period,data=summary)
summary(model)

model <- aov(sd~period,data=summary)
summary(model)
model <- lm(mean~period,data=summary)
summary(model)

summary<-combined%>%
  mutate(period=ifelse(year>=2012, "After", "Before"))%>%
  mutate(type=ifelse(mohns>=0.1, "Unstable",ifelse(mohns<0.01, "Stable", "Moderate")))%>%
  group_by(names, period, type,Species)%>%
#  summarise(mean=mean(value), sd=sd(value))%>%
  mutate(grouping=paste(period, type))


ggplot(data = summary, aes(x = value, y =period, col)) +
  geom_violin()+
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA) 
  
