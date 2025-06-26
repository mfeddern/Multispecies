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
yt_env <- data.frame(read.csv("Data/Yellowtail/2024Env-annual-yellowtail_GLORYS_UNSTANDARDIZED.csv"))
#yellowtail_recdev
yt_env_sl<- yt_env%>%
  select(-c(year,Year,Year.1, X,ZOObenS,ZOOpjuvS,ZOObenN,ZOOpjuvN,CHLpjuv,PPpjuv,bakun_sti))%>%
  scale()%>%
  cbind(year = yt_env$year)

write.csv(yt_env_sl,"Data/Yellowtail/2024Env-annual-yellowtail_GLORYS_STANDARDIZED.csv")

yt_recdev_final <- data.frame(read.csv("Data/Yellowtail/2025_yt_final.csv"))%>%
  mutate(Datatreatment="2025 Final")
yt_recdev_draft <- data.frame(read.csv("Data/Yellowtail/yt_RecruitmentDeviations2025draft.csv"))
yt_dat <- yt_recdev_draft%>%bind_rows(yt_recdev_final)%>%
  left_join(data.frame(yt_env_sl))

write.csv(yt_dat,"Data/Yellowtail/yt_fulldataset.csv")

