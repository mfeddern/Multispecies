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

yt_env <- data.frame(read.csv("Data/2024Env-annual-yellowtail_GLORYS_UNSTANDARDIZED.csv"))
#yellowtail_recdev
yt_env_sl<- yt_env%>%
  select(-c(year,Year,Year.1, X,ZOObenS,ZOOpjuvS,ZOObenN,ZOOpjuvN,CHLpjuv,PPpjuv,bakun_sti))%>%
  scale()%>%
  cbind(year = yt_env$year)
