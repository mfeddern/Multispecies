library(dplyr)
library(reshape2)
library(bayesdfa)
library(MCMCvis)
library(ggplot2)
library(stringr)
library(ggpubr)

#### Data ####
yrlast<-2018
mohns_dat<- readRDS("Output/mohns.rds")%>%
  mutate(mohns=abs(mohns))
ps_dat <- data.frame(read.csv("Data/Petrale/DATA_Combined_glorys_petrale_STANDARDIZED.csv"))%>%
  select(-c(X, year.1))%>%
  # filter(type=="Main_RecrDev")%>%
  filter(year>1994&year<=yrlast)
sb_dat <- data.frame(read.csv("Data/Sablefish/data-combined-glorys-sablefish_STANDARDIZED.csv"))%>%
  select(-c(X))%>%
  filter(type=="Main_RecrDev")%>%
  filter(year>1994&year<=yrlast)
yt_dat <- data.frame(read.csv("Data/Yellowtail/yt_fulldataset_STANDARDIZED.csv"))%>%
  select(-c(X, hci2_pjuv,hci1_pjuv, lusi_annual, hci2_larv,hci1_larv))%>%
  filter(type=="Main")%>%
  filter(Datatreatment=="2025 Final"&year>1994&year<=yrlast)

subsethk<-readRDS("Output/Data/hakesubset.rds")
hk_dat <- data.frame(read.csv("Data/Hake/DATA_Combined_glorys_hake_STANDARDIZED.csv"))%>%
  select(-X)%>%
  filter(type=="Main_RecrDev"&year>1993&year<=2022)%>%
  select(all_of(subsethk), year)

ax<- 20
ti<-24
wid <- 28
plot_trends2 <- function(rotated_modelfit,
                         years = NULL,
                         highlight_outliers = FALSE,
                         threshold = 0.01) {
  rotated <- rotated_modelfit
  df <- dfa_trends(rotated, years = years)
  
  # make faceted ribbon plot of trends
  p1 <- ggplot(df, aes_string(x = "time", y = "estimate")) +
    geom_ribbon(aes_string(ymin = "lower", ymax = "upper"), alpha = 0.4) +
    geom_line() +
    facet_wrap("trend_number") +
    xlab("Time") +
    ylab("")+
    
    theme(axis.text=element_text(size=ax),
          axis.title=element_text(size=ti,face="bold"))+
    theme_bw()
  
  if (highlight_outliers) {
    swans <- find_swans(rotated, threshold = threshold)
    df$outliers <- swans$below_threshold
    p1 <- p1 + geom_point(data = df[which(df$outliers), ], color = "red")
  }
  
  p1
}


plot_loadings2 <- function(rotated_modelfit,
                           names = NULL,
                           facet = TRUE,
                           violin = TRUE,
                           conf_level = 0.95,
                           threshold = NULL) {
  v <- dfa_loadings(rotated_modelfit,
                    summary = FALSE,
                    names = names,
                    conf_level = conf_level
  )
  df <- dfa_loadings(rotated_modelfit,
                     summary = TRUE,
                     names = names,
                     conf_level = conf_level
  )
  
  # filter values below threshold
  if (!is.null(threshold)) {
    df <- df[df$prob_diff0 >= threshold, ]
    v <- v[v$prob_diff0 >= threshold, ]
  }
  
  if (!violin) {
    p1 <- ggplot(df, aes_string(
      x = "name", y = "median", col = "trend",
      alpha = "prob_diff0"
    )) +
      geom_point(size = 3, position = position_dodge(0.3)) +
      geom_errorbar(aes_string(ymin = "lower", ymax = "upper"),
                    position = position_dodge(0.3), width = 0
      ) +
      geom_hline(yintercept = 0, lty = 2) +
      coord_flip() +
      xlab("Time Series") +
      ylab("Loading")+
      guides(fill="none", alpha='none')+
      scale_x_discrete(labels = function(x) str_wrap(x, width = wid))+
      theme(legend.position="none")+
      theme_bw()
  }
  
  if (violin) {
    p1 <- ggplot(v, aes_string(
      x = "name", y = "loading", fill = "trend",
      alpha = "prob_diff0"
    )) +
      geom_violin(color = NA) +
      geom_hline(yintercept = 0, lty = 2) +
      coord_flip() +
      xlab("Time Series") +
      ylab("Loading")+
      theme(axis.text=element_text(size=ax),
            axis.title=element_text(size=ti,face="bold"))+
      guides(fill="none", alpha='none')+
      scale_x_discrete(labels = function(x) str_wrap(x, width = wid))+
      theme_bw()
  }
  
  if (facet) {
    p1 <- p1 + facet_wrap(~trend, scales = "free_x")
  }
  
  p1
}




#### Wrangle Data ####
hk_combined <-pivot_longer(hk_dat,-year)%>%
  mutate(Species="Hake")%>%
  rename(variable=name)%>%
  left_join(mohns)%>%
  mutate(names=paste("Hake_", variable))

ps_combined <-pivot_longer(ps_dat%>%select(c(-type, -Y_rec, -sd)), -year)%>%
  mutate(Species="Petrale Sole")%>%
  rename(variable=name)%>%
  left_join(mohns_dat)%>%
  mutate(names=paste("PetraleSole_", variable))

sb_combined<-pivot_longer(sb_dat%>%select(c(-type, -Y_rec, -sd)), -year)%>%
  mutate(Species="Sablefish")%>%
  rename(variable=name)%>%
  left_join(mohns)%>%
  mutate(names=paste("Sablefish_", variable))


yt_combined <-pivot_longer(yt_dat%>%select(c(-type, -Y_rec, -sd,-Datatreatment)), -year)%>%
  mutate(Species="Yellowtail")%>%
  rename(variable=name)%>%
  left_join(mohns)%>%
  mutate(names=paste("Yellowtail_", variable))

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
#### UNSTABLE VARIABLES ####
combined_unstable<-combined%>%
  filter(mohns>0.1)
n1 <-unique(combined_unstable$names)

dat.unstable<-all_dat[which(colnames(all_dat)%in%c(n1,"year"))]
  
remelt = melt(dat.unstable,id.vars = "year")
names(remelt)<-c("year","code","value")
Y <- dcast(remelt, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])


n_chains = 3
n_iter =8000

options(mc.cores = parallel::detectCores())

# load environmental covariate data, and average over spring months to have one value per covar per year
nspp=dim(Y)[1]
nyear=dim(Y)[2]
ntrend=1
n_row=nspp*nyear
n_row_pro=ntrend*nyear

model_df = expand.grid(estimate_trend_ma = FALSE,
                       estimate_trend_ar = TRUE, est_nu = TRUE, estimate_process_sigma = c(FALSE),
                       var_index = c("survey"), num_trends = 1,
                       elpd_loo = TRUE, se_elpd_loo=TRUE)
varIndx =c(rep(1,length(n1)))
##### one trend #####

fit.mod.unstable.v1 = fit_dfa(y = Y,
                          num_trends =1,
                          iter=n_iter,
                          varIndx = varIndx,
                          chains=n_chains, estimate_nu=TRUE,
                          estimate_trend_ma = FALSE,
                          estimate_trend_ar = TRUE,
                          estimate_process_sigma = FALSE,
                          seed=123)

pars = rstan::extract(fit.mod.unstable.v1$model)
r.unstable.v1 <- rotate_trends(fit.mod.unstable.v1)
p.unstable.v1 <- plot_trends2(r.unstable.v1,years =dat.unstable$year)
p.unstable.v1
l.unstable.v1 <- plot_loadings2(r.unstable.v1,names=n1)
l.unstable.v1
is_converged(fit.mod.unstable.v1)
summary(fit.mod.unstable.v1)

loo1 <- loo(fit.mod.unstable.v1)


##### two trends #####
fit.mod.unstable.v2 = fit_dfa(y = Y,
                              num_trends =2,
                              iter=n_iter,
                              varIndx = varIndx,
                              chains=n_chains, estimate_nu=model_df$est_nu[1],
                              estimate_trend_ma = model_df$estimate_trend_ma[1],
                              estimate_trend_ar = model_df$estimate_trend_ar[1],
                              estimate_process_sigma = model_df$estimate_process_sigma[1],
                              seed=123)

pars = rstan::extract(fit.mod.unstable.v2$model)
r.unstable.v2 <- rotate_trends(fit.mod.unstable.v2)
p.unstable.v2 <- plot_trends2(r.unstable.v2,years =dat.unstable$year)
p.unstable.v2
l.unstable.v2 <- plot_loadings2(r.unstable.v2,names=n1)
l.unstable.v2
is_converged(fit.mod.unstable.v2)
summary(fit.mod.unstable.v2)
loo2 <- loo(fit.mod.unstable.v2)
loo2$estimates
loo1$estimates

#### STABLE VARIABLES ####

combined_stable<-combined%>%
  filter(mohns>0.01)
n1 <-unique(combined_unstable$names)
dat.stable<-all_dat[which(colnames(all_dat)%in%c(n1,"year"))]

remelt = melt(dat.stable,id.vars = "year")
names(remelt)<-c("year","code","value")
Y <- dcast(remelt, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])


n_chains = 3
n_iter = 8000

options(mc.cores = parallel::detectCores())

# load environmental covariate data, and average over spring months to have one value per covar per year
nspp=dim(Y)[1]
nyear=dim(Y)[2]
ntrend=1
n_row=nspp*nyear
n_row_pro=ntrend*nyear


varIndx =c(rep(1,length(n1)))
fit.mod.stable.v1 = fit_dfa(y = Y,
                          num_trends =1,
                          iter=n_iter,
                          varIndx = varIndx,
                          chains=n_chains, estimate_nu=model_df$est_nu[1],
                          estimate_trend_ma = model_df$estimate_trend_ma[1],
                          estimate_trend_ar = model_df$estimate_trend_ar[1],
                          estimate_process_sigma = model_df$estimate_process_sigma[1],
                          seed=123)

pars = rstan::extract(fit.mod.stable.v1$model)
r.stable.v1 <- rotate_trends(fit.mod.stable.v1)
p.stable.v1 <- plot_trends2(r.stable.v1,years =dat.stable$year)
p.stable.v1
p.unstable.v1
l.stable.v1 <- plot_loadings2(r.stable.v1,names=n1)
l.stable.v1
is_converged(fit.mod.stable.v1)
summary(fit.mod.stable.v1)
loo1st <- loo(fit.mod.stable.v1)
##### two trends #####
fit.mod.stable.v2 = fit_dfa(y = Y,
                              num_trends =2,
                              iter=n_iter,
                              varIndx = varIndx,
                              chains=n_chains, estimate_nu=model_df$est_nu[1],
                              estimate_trend_ma = model_df$estimate_trend_ma[1],
                              estimate_trend_ar = model_df$estimate_trend_ar[1],
                              estimate_process_sigma = model_df$estimate_process_sigma[1],
                              seed=123)

pars = rstan::extract(fit.mod.stable.v2$model)
r.stable.v2 <- rotate_trends(fit.mod.stable.v2)
p.stable.v2 <- plot_trends2(r.stable.v2,years =dat.stable$year)
p.stable.v2
l.stable.v2 <- plot_loadings2(r.stable.v2,names=n1)
l.stable.v2
is_converged(fit.mod.stable.v2)
summary(fit.mod.stable.v2)
loo2 <- loo(fit.mod.stable.v2)
loo2$estimates
loo1$estimates


stable<-annotate_figure(ggarrange(p.stable.v1,l.stable.v1, labels = c("C.", "D.")), top ="Stable Indices")
unstable<-annotate_figure(ggarrange(p.unstable.v1,l.unstable.v1, labels = c("A.", "B.")), top ="Unstable Indices")


dfas<-ggarrange(unstable, stable, nrow=2, ncol=1)

pdf(file = "Output/Figures/dfas.pdf", width =12, height =9)
dfas
dev.off()

png(file = "Output/Figures/dfas.png",width = 1200, height = 900, res = 100)
dfas
dev.off()

