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


#### Yellowtail Data ####
# developing a standardized dataset
yt_env <- data.frame(read.csv("Data/Yellowtail/2024Env-annual-yellowtail_GLORYS_UNSTANDARDIZED.csv"))%>%
  select(-c(Year,Year.1, X,ZOObenS,ZOOpjuvS,ZOObenN,ZOOpjuvN,CHLpjuv,PPpjuv,bakun_sti))
#yellowtail_recdev
yt_env_sl<- yt_env%>%
  select(-c(year))%>%
  scale()%>%
  cbind(year = yt_env$year)
write.csv(yt_env_sl,"Data/Yellowtail/2024Env-annual-yellowtail_GLORYS_STANDARDIZED.csv")

#combining standardized dataset with all rec devs
yt_recdev_final <- data.frame(read.csv("Data/Yellowtail/2025_yt_final.csv"))%>%
  mutate(Datatreatment="2025 Final")
yt_recdev_draft <- data.frame(read.csv("Data/Yellowtail/yt_RecruitmentDeviations2025draft.csv"))
full_recdev<-  yt_recdev_draft%>%bind_rows(yt_recdev_final)
yt_dat_std <-full_recdev%>%
  left_join(data.frame(yt_env_sl))
write.csv(yt_dat_std,"Data/Yellowtail/yt_fulldataset_STANDARDIZED.csv")

#combining all recdevs with unstandardized dataset
yt_dat_unst <-full_recdev%>%
  left_join(data.frame(yt_env))
write.csv(yt_dat_unst,"Data/Yellowtail/yt_fulldataset_UNSTANDARDIZED.csv")

### Hake Data ###
hk_recdev <- data.frame(read.csv("Data/Hake/DATA_Combined_glorys_hake_UNSTANDARDIZED.csv"))
hk_env_st <- hk_recdev%>%
  select(-c(year, type, Y_rec, sd))%>%
  scale()%>%
  cbind(year = hk_recdev$year)
hk_dat_st <-data.frame(hk_env_st)%>%
  left_join(data.frame(hk_recdev%>%select(type, year, Y_rec,sd)))
write.csv(hk_dat_st,"Data/Hake/DATA_Combined_glorys_hake_STANDARDIZED.csv")

### Petrale Sole Data ###

ps_recdev <- data.frame(read.csv("Data/Petrale/DATA_Combined_glorys_petrale_UNSTANDARDIZED.csv"))
ps_env_st <- ps_recdev%>%
  select(-c(year, type, Y_rec, sd))%>%
  scale()%>%
  cbind(year = ps_recdev$year)
ps_dat_st <-data.frame(ps_env_st)%>%
  left_join(data.frame(ps_recdev%>%select(type, year, Y_rec,sd)))
write.csv(ps_dat_st,"Data/Petrale/DATA_Combined_glorys_petrale_STANDARDIZED.csv")

### Sablefish Data ###

sb_recdev <- data.frame(read.csv("Data/Sablefish/data-combined-glorys-sablefish_UNSTANDARDIZED.csv"))
sb_env_st <-sb_recdev%>%
  select(-c(year, type, Y_rec, sd))%>%
  scale()%>%
  cbind(year = sb_recdev$year)
sb_dat_st <-data.frame(sb_env_st)%>%
  left_join(data.frame(sb_recdev%>%select(type, year, Y_rec,sd)))
write.csv(ps_dat_st,"Data/Sablefish/data-combined-glorys-sablefish_STANDARDIZED.csv")
