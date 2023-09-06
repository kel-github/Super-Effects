#########################################################################
# plot fx size/p vs dist measures
# K. Garner 
# for the supfx project, can we predict error in effect size estimates
# using higher order moment data
# note: If skewness is less than -1 or greater than 1, the distribution is highly skewed. 
# If skewness is between -1 and -0.5 or between 0.5 and 1, the distribution is moderately skewed. 
# If skewness is between -0.5 and 0.5, the distribution is approximately symmetric

rm(list=ls())

library(MASS)
library(tidyverse) # for data wrangling
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(moments)
library(stringr)
library(car)
library(GGally)
library(Hmisc)
library(lme4)
library(caret)
source("efilids_functions.R") # custom functions written for this project
source("moments_functions.R")
source("R_rainclouds.R") # functions for plotting

set.seed(42) # meaning of life
################################################################################
# Run things
################################################################################
not_new <- T
# define sub Ns
datpath <- "../data"
tasks <- c("SD", "AB", "SRT", "CC")
zipnms <- c("SD_wv.zip", "EPSAB_wv.zip", "SRT_wv.zip", "CC_wv.zip")
sub.Ns <- round(exp(seq(log(13), log(313), length.out = 20)))
maxN = 313
Ns <- list(c(42, 313), c(25, 313), c(36, 313), c(25, 313)) 
train <- trainControl(method = "cv", number = 10)

# first get a list, where each element is the data from one task
if (!not_new){ 
dat <- mapply(function(x,y,z) do.call(rbind, lapply(z, 
                                                  get_dat_4corrs, 
                                                  datpath=datpath, 
                                                  task=x, 
                                                  zipnm=y)),
                x=tasks,
                y = zipnms,
                z = Ns)
names(dat) <- tasks
save(dat, file = paste(datpath, "fx_by_dst_regression.RData", sep="/"))
} else {
  load(paste(datpath, "fx_by_dst_regression.RData", sep="/"))
}
##############
# now for each task, compute the appropriate distance variables
##############################################################
# SD
##############################################################
mu <- get_ev_non_norm(dat[["SD"]]$esz_ME[dat[["SD"]]$n == maxN])
# now compute the distance of each value of the fx size
dat[["SD"]]$esz_ME_dist <- dat[["SD"]]$esz_ME - mu
# now do the interaction effect
mu <- get_ev_non_norm(dat[["SD"]]$esz_int[dat[["SD"]]$n == maxN])
dat[["SD"]]$esz_int_dist <- dat[["SD"]]$esz_int - mu

##############################################################
# AB
##############################################################
mu <- get_ev_non_norm(dat[["AB"]]$esz[dat[["AB"]]$n == maxN])
# now compute the distance of each value of the fx size
dat[["AB"]]$esz_dist <- dat[["AB"]]$esz - mu

##############################################################
# SRT
##############################################################
mu <- get_ev_non_norm(dat[["SRT"]]$esz[dat[["SRT"]]$n == maxN])
mu
# now compute the distance of each value of the fx size
dat[["SRT"]]$esz_dist <- dat[["SRT"]]$esz - mu

##############################################################
# SD
##############################################################
mu <- get_ev_non_norm(dat[["CC"]]$esz_ME[dat[["CC"]]$n == maxN])
# now compute the distance of each value of the fx size
dat[["CC"]]$esz_ME_dist <- dat[["CC"]]$esz_ME - mu
# now do the interaction effect
mu <- get_ev_non_norm(dat[["CC"]]$esz_int[dat[["CC"]]$n == maxN])
dat[["CC"]]$esz_int_dist <- dat[["CC"]]$esz_int - mu

#############################################################
# Things I can say so far:
# for all tasks, skew, kurtosis, correlation between 
# variables and subject variability predicted error in effect 
# size estimates
#############################################################
tstN <- 25
full_model <- lm(esz_dist ~ AB.skew + AB.k + sigma_mu + skew_mu + k_mu, 
                 data = dat[['AB']] %>% filter(n == tstN))
summary(full_model)
# cross validate
AB_model <- train(esz_dist ~ AB.skew + AB.k + sigma_mu + skew_mu + k_mu, 
                    data = dat[['AB']] %>% filter(n == tstN),
                    method = "lm",
                    trControl = train)

AB_model
AB_model$finalModel

# now take the model and add group interaction terms
AB_dat <- dat[["AB"]]
AB_grp_int <- lm(esz_dist ~ AB.skew + AB.k + sigma_mu + k_mu +
                            AB.skew*n + AB.k*n + sigma_mu*n + k_mu*n ,
                        data = AB_dat) # no group interaction
plot(AB_grp_int, which=1)
summary(AB_grp_int)

### now compute tthe correlation between effect size error
# and skew
AB_skew_cor <- lm(esz_dist ~ AB.skew, data = dat[["AB"]] %>% filter(n == tstN))
summary(AB_skew_cor)

AB_skew_model <- train(esz_dist ~ AB.skew, 
                      data = dat[['AB']] %>% filter(n == tstN),
                      method = "lm",
                      trControl = train)

# and kurtosis
AB_k_cor <- lm(esz_dist ~ AB.k, data = dat[["AB"]] %>% filter(n == tstN))
summary(AB_k_cor)

AB_k_model <- train(esz_dist ~ AB.k, 
                    data = dat[['AB']] %>% filter(n == tstN),
                    method = "lm",
                    trControl = train)

# and within participant skew
AB_musku <- lm(esz_dist ~ skew_mu, data = dat[["AB"]] %>% filter(n == tstN))

# now predict danger for skew
# the effect size for N=313
AB_95 <- quantile(dat[["AB"]]$esz_dist[dat[["AB"]]$n == 313], c(.025, .975))
AB_95 # this gives the y-values

# the prediction is 
# pred_y = int + beta*X
# beta*X = pred_y - int
# X = (pred_y - int)/beta
ABX <- (AB_95 - AB_skew_cor$coefficients["(Intercept)"]) / 
               AB_skew_cor$coefficients["AB.skew"]

# now compute the standard error pf the estimate
# sy.x = sy*sqrt(1-r^2)
ABX_SE <- sd(dat[["AB"]]$esz_dist[dat[["AB"]]$n == tstN])*sqrt(1-0.09572)

ABmuskewX <- (AB_95 - AB_musku$coefficients["(Intercept)"]) /
                       AB_musku$coefficients["skew_mu"]

ABmuskewX_SE <- sd(dat[["AB"]]$esz_dist[dat[["AB"]]$n == tstN])*sqrt(1-0.1822)

AB_pred_list <- list(full=full_model, cv=AB_model, grp=AB_grp_int, 
                     skew=AB_skew_cor, k=AB_k_cor, mu_skew=AB_musku,
                     X=ABX, SE=ABX_SE)
#############################################################
# SRT
#############################################################
# SRT is a little more complex, got some work to do on the variables
# first
# 1) get EV of distribution
tstN <- 36
# 2) allocate data to a df and make variables behave
SRT_df <- dat[["SRT"]] %>% filter(n == tstN) %>%
                            select(esz_dist, SRT.skew,
                                    SRT.k, sigma_mu, skew_mu, k_mu)
ggpairs(SRT_df)
hist(SRT_df$esz_dist, breaks = 30) # positively skewed
apply(SRT_df %>% select(SRT.k, sigma_mu, skew_mu, k_mu), 2, min) # all +ve 
# first make all the values positive and log transform
x = 0.01 - min(SRT_df$esz_dist)
SRT_df <- SRT_df %>% mutate(esz_dist_t = log(esz_dist + x),
                            SRT.k_t = log(SRT.k), # same
                            sigma_mu_t = log(sigma_mu), # same
                            skew_mu_t = log(skew_mu),
                            k_mu_t = log(k_mu)) # same

SRT_df <- SRT_df %>% select(esz_dist_t, SRT.skew, SRT.k_t, 
                            sigma_mu_t, skew_mu_t, k_mu_t)

SRT_model <- train(esz_dist_t ~ SRT.skew + SRT.k_t + sigma_mu_t + skew_mu_t + k_mu_t, 
                          data = SRT_df,
                          method = "lm",
                          trControl = train)

SRT_model
SRT_model$finalModel

SRT_full_model <- lm(esz_dist_t ~ SRT.skew + SRT.k_t + sigma_mu_t + skew_mu_t + k_mu_t, 
                     data = SRT_df)
summary(SRT_full_model)

# now make group level dataframe
SRT_grp_df <- dat[["SRT"]] %>%
            select(esz_dist, n, SRT.skew,
              SRT.k, sigma_mu, skew_mu, k_mu)
apply(SRT_grp_df %>% select(SRT.k, sigma_mu, skew_mu, k_mu), 2, min) # all +ve 
# first make all the values positive and log transform
SRT_grp_df <- SRT_grp_df %>% mutate(esz_dist_t = log(esz_dist + x),
                                    SRT.k_t = log(SRT.k), # same
                                    sigma_mu_t = log(sigma_mu), # same
                                    skew_mu_t = log(skew_mu),
                                    k_mu_t = log(k_mu)) # same
SRT_grp_int <- lm(esz_dist ~ SRT.skew + SRT.k_t + sigma_mu_t + k_mu_t +
                                        SRT.skew*n + SRT.k_t*n + sigma_mu_t*n + 
                                        k_mu_t*n ,
                         data = SRT_grp_df) 
summary(SRT_grp_int)

### now compute tthe correlation between effect size error
# and skew
SRT_skew_cor <- lm(esz_dist_t ~ SRT.skew, data = SRT_df)
summary(SRT_skew_cor)

# and kurtosis
SRT_k_cor <- lm(esz_dist_t ~ SRT.k_t, data = SRT_df)
summary(SRT_k_cor)

SRT_mu_skew <- lm(esz_dist_t ~ skew_mu_t, data = SRT_df)

# now predict danger for skew
# the effect size for N=313
SRT_95 <- log(quantile(dat[["SRT"]]$esz_dist[dat[["SRT"]]$n == 313], c(.025, .975)) + x)
# now SRT quantiles are on the same scale as the regression model
SRT_95 # this gives the y-values

# the prediction is 
# pred_y = int + beta*X
# beta*X = pred_y - int
# X = (pred_y - int)/beta
SRTX <- (SRT_95 - SRT_skew_cor$coefficients["(Intercept)"]) / 
              SRT_skew_cor$coefficients["SRT.skew"]

# now compute the standard error pf the estimate
# sy.x = sy*sqrt(1-r^2)
SRTX_SE <- sd(log(dat[["SRT"]]$esz_dist[dat[["SRT"]]$n == tstN]+x))*sqrt(1-0.04297) # make sure you've changed R^2

SRT_pred_list <- list(full=SRT_full_model, cv=SRT_model, grp=SRT_grp_int,
                      skew=SRT_skew_cor, k=SRT_k_cor, mu_skew=SRT_mu_skew,
                      X=SRTX, SE=SRTX_SE)

#############################################################
# SD
#############################################################
# first look at data for ME and int effects, transform variables
# to make them behave well
############################################################
# SD ME
############################################################
tstN <- 42
SD_ME_df <- dat[["SD"]] %>% filter(n == tstN) %>%
              select(esz_ME_dist, ME_skew,
                     ME_k, sigma_mu, skew_mu, k_mu)

apply(SD_ME_df, 2, min) 
# k is clearly non-normal
SD_ME_mod_df <- SD_ME_df %>% select(esz_ME_dist, ME_skew,
                                    ME_k, sigma_mu, skew_mu, k_mu) %>%
                              mutate(ME_k_t = log(ME_k)) %>%
                              select(esz_ME_dist, ME_skew,
                                     ME_k_t, sigma_mu, skew_mu, k_mu)
ggpairs(SD_ME_mod_df)

SD_ME_model <- train(esz_ME_dist ~ ME_skew + ME_k_t + sigma_mu + skew_mu + k_mu, 
                      data = SD_ME_mod_df,
                      method = "lm",
                      trControl = train)
SD_ME_model
SD_ME_model$finalModel

SD_ME_full_model <- lm(esz_ME_dist ~ ME_skew + ME_k_t + sigma_mu + skew_mu + k_mu, 
                       data = SD_ME_mod_df)
summary(SD_ME_full_model)



# now do group comparison
SD_ME_grp_df <- dat[["SD"]] %>% select(esz_ME_dist, n, ME_skew,
                                        ME_k, sigma_mu, skew_mu, k_mu) %>% 
                                 mutate(ME_k_t = log(ME_k)) %>%
                                 select(esz_ME_dist, n, ME_skew,
                                        ME_k_t, sigma_mu, skew_mu, k_mu)
SD_ME_grp_int <- lm(esz_ME_dist ~ ME_skew + ME_k_t + sigma_mu + skew_mu + k_mu +
                                    ME_skew*n + ME_k_t*n + sigma_mu*n +
                                      skew_mu*n + k_mu*n,
                    data = SD_ME_grp_df)
plot(SD_ME_grp_int, which=1)
summary(SD_ME_grp_int)

### now compute tthe correlation between effect size error
# and skew
SD_ME_skew_cor <- lm(esz_ME_dist ~ ME_skew, data = SD_ME_mod_df)
summary(SD_ME_skew_cor)

# and kurtosis
SD_ME_k_cor <- lm(esz_ME_dist ~ ME_k_t, data = SD_ME_mod_df)
summary(SD_ME_k_cor)

SD_mu_skew <- lm(esz_ME_dist ~ skew_mu, data = SD_ME_mod_df)
summary(SD_mu_skew)

# now predict danger for skew
# the effect size for N=313
SD_ME_95 <- quantile(dat[["SD"]]$esz_ME_dist[dat[["SD"]]$n == 313], c(.025, .975))
# now SRT quantiles are on the same scale as the regression model
SD_ME_95 # this gives the y-values

# the prediction is 
# pred_y = int + beta*X
# beta*X = pred_y - int
# X = (pred_y - int)/beta
SD_ME_X <- (SD_ME_95 - SD_ME_skew_cor$coefficients["(Intercept)"]) / 
  SD_ME_skew_cor$coefficients["ME_skew"]

# now compute the standard error pf the estimate
# sy.x = sy*sqrt(1-r^2)
SD_ME_SE <- sd(dat[["SD"]]$esz_ME_dist[dat[["SD"]]$n != 313])*sqrt(1-0.04726)

SD_ME_pred_list <- list(full=SD_ME_full_model, cv=SD_ME_model, grp=SD_ME_grp_int,
                        skew=SD_ME_skew_cor, k=SD_ME_k_cor, mu_skew=SD_mu_skew,
                        X=SD_ME_X, SE=SD_ME_SE)

############################################################
# SD INT
############################################################
SD_int_df <- dat[["SD"]] %>% filter(n == tstN) %>%
                    select(esz_int_dist, int_skew,
                            int_k, sigma_mu, skew_mu, k_mu)
ggpairs(SD_int_df)
# esz_int_dist and int_k are both highly skewed
x <- 0.01 - min(SD_int_df$esz_int_dist)
apply(SD_int_df %>% select(int_k, k_mu), 2, min)
SD_int_mod_df <- SD_int_df %>%  mutate(esz_int_dist_t = log(esz_int_dist + x),
                                       int_k_t = log(int_k),
                                       k_mu_t = log(k_mu)) %>%
                                select(esz_int_dist_t, int_skew, int_k_t,
                                       sigma_mu, skew_mu, k_mu_t)

SD_int_model <- train(esz_int_dist_t ~ int_skew + int_k_t + sigma_mu + skew_mu + k_mu_t, 
                          data = SD_int_mod_df,
                          method = "lm",
                          trControl = train)
SD_int_model
SD_int_model$finalModel

SD_int_full_model <- lm(esz_int_dist_t ~ int_skew + int_k_t + sigma_mu + skew_mu + k_mu_t, 
                        data = SD_int_mod_df)
summary(SD_int_full_model)

# group comp:
SD_int_grp_df <- dat[["SD"]] %>% select(esz_int_dist, n, int_skew,
                                          int_k, sigma_mu, skew_mu, k_mu) %>%  
                                mutate(esz_int_dist_t = log(esz_int_dist + x),
                                          int_k_t = log(int_k),
                                          k_mu_t = log(k_mu)) %>%
                               select(esz_int_dist_t, n, int_skew, int_k_t,
                                          sigma_mu, skew_mu, k_mu_t)
SD_int_grp_int <- lm(esz_int_dist_t ~ int_skew + int_k_t + sigma_mu + skew_mu + 
                                    k_mu_t + int_skew*n + int_k_t*n + 
                                    sigma_mu*n + skew_mu*n +
                                    k_mu_t*n, data=SD_int_grp_df)
plot(SD_int_grp_int, which=1)
summary(SD_int_grp_int)

### now compute tthe correlation between effect size error
# and skew
SD_int_skew_cor <- lm(esz_int_dist_t ~ int_skew, data = SD_int_mod_df)
summary(SD_int_skew_cor)

# and kurtosis
SD_int_k_cor <- lm(esz_int_dist_t ~ int_k_t, data = SD_int_mod_df)
summary(SD_ME_k_cor)

SD_int_mu_skew <- lm(esz_int_dist_t ~ skew_mu, data = SD_int_mod_df)
summary(SD_int_mu_skew)

SD_int_pred_list <- list(full=SD_int_full_model, cv=SD_int_model, grp=SD_int_grp_int,
                          skew=SD_int_skew_cor, k=SD_int_k_cor, mu_skew=SD_int_mu_skew)

#############################################################
# CC
#############################################################
# first look at data for ME and int effects, transform variables
# to make them behave well
############################################################
# CC ME
############################################################
tstN <- 25
CC_ME_df <- dat[["CC"]] %>% filter(n == tstN) %>%
                select(esz_ME_dist, ME_skew,
                ME_k, sigma_mu, skew_mu, k_mu)
ggpairs(CC_ME_df) # ME_k, sigma_mu & k_mu are most
# skewed, gonna start by transforming those
apply(CC_ME_df, 2, min)
x = 0.01 - min(CC_ME_df$ME_k) 
CC_ME_df <- CC_ME_df %>% mutate(ME_k_t = log(ME_k + x),
                                sigma_mu_t = (sigma_mu),
                                k_mu_t = log(k_mu)) %>%
                         select(esz_ME_dist, ME_skew, ME_k_t, sigma_mu_t,
                                skew_mu, k_mu_t)

CC_ME_mod <- train(esz_ME_dist ~ ME_skew + ME_k_t + sigma_mu_t + skew_mu + 
                     k_mu_t, 
                    data = CC_ME_df,
                    method = "lm",
                    trControl = train)
CC_ME_mod
CC_ME_mod$finalModel

CC_ME_mod_full <- lm(esz_ME_dist ~ ME_skew + ME_k_t + sigma_mu_t + skew_mu + 
                       k_mu_t, 
                     data = CC_ME_df)
summary(CC_ME_mod_full)

# now do group comparison
CC_ME_grp_df <- dat[["CC"]] %>% select(esz_ME_dist, n, ME_skew, ME_k, 
                                         sigma_mu, skew_mu, k_mu) %>% 
                                mutate(ME_k_t = log(ME_k + x),
                                       sigma_mu_t = log(sigma_mu),
                                       k_mu_t = log(k_mu)) %>%
                               select(esz_ME_dist, n, ME_skew, ME_k_t, sigma_mu_t,
                                       skew_mu, k_mu_t)
CC_ME_grp_int <- lm(esz_ME_dist ~ ME_skew + ME_k_t + sigma_mu_t + skew_mu +
                                      k_mu_t + ME_skew*n + ME_k_t*n + 
                                      sigma_mu_t*n + skew_mu*n + k_mu_t*n,
                                data = CC_ME_grp_df) 
plot(CC_ME_grp_int, which=1)
summary(CC_ME_grp_int)

#### 16%
CC_ME_mu_skew <- lm(esz_ME_dist ~ sigma_mu_t, data = CC_ME_df)
summary(CC_ME_mu_skew)

CC_ME_pred_list <- list(full=CC_ME_mod_full, cv=CC_ME_mod, grp=CC_ME_grp_int, 
                        mu_skew=CC_ME_mu_skew)


###############################################################################
## CC Int
##############################################################################
CC_int_df <- dat[["CC"]] %>% filter(n == tstN) %>%
                select(esz_int_dist, int_skew,
                       int_k, sigma_mu, skew_mu, k_mu)
ggpairs(CC_int_df) # esz_int_dist, int_k, sigma_mu & k_mu are most
# skewed, gonna transform those
apply(CC_int_df, 2, min)
x = 0.01 - min(CC_int_df$esz_int_dist)
CC_int_df <- CC_int_df %>% mutate(esz_int_dist_t = log(esz_int_dist + x),
                                  int_k_t = log(int_k),
                                  sigma_mu_t = log(sigma_mu),
                                  k_mu_t = log(k_mu)) %>%
                  select(esz_int_dist_t, int_skew, int_k_t, sigma_mu_t,
                          skew_mu, k_mu_t)

CC_int_mod <- train(esz_int_dist_t ~ int_skew + int_k_t + sigma_mu_t + skew_mu + 
                     k_mu_t, 
                     data = CC_int_df,
                     method = "lm",
                     trControl = train)
CC_int_mod
CC_int_mod$finalModel

CC_int_mod_full <- lm(esz_int_dist_t ~ int_skew + int_k_t + sigma_mu_t + skew_mu + 
                        k_mu_t, 
                      data = CC_int_df)
summary(CC_int_mod_full)

# now do group comparison
CC_int_grp_df <- dat[["CC"]] %>% select(esz_int_dist, n, int_skew, int_k,
                                        sigma_mu, skew_mu, k_mu) %>%
                                  mutate(esz_int_dist_t = log(esz_int_dist + x),
                                         int_k_t = log(int_k),
                                         sigma_mu_t = log(sigma_mu),
                                         k_mu_t = log(k_mu))
CC_int_grp_int <- lm(esz_int_dist_t ~ int_skew + int_k_t + sigma_mu_t + skew_mu + 
                       k_mu_t + int_skew*n + int_k_t*n + sigma_mu_t*n + skew_mu*n +
                       k_mu_t*n,
                    data = CC_int_grp_df) 
plot(CC_int_grp_int, which=1)
summary(CC_int_grp_int) # no effect of N

# .004
CC_int_mu_skew <- lm(esz_int_dist_t ~ sigma_mu_t, data = CC_int_df)
summary(CC_int_mu_skew)

CC_int_pred_list <- list(full=CC_int_mod_full, cv=CC_int_mod, grp=CC_int_grp_int,
                         mu_skew=CC_int_mu_skew)
##############################################################################
# save a list for the manuscript
predictions_all_tasks <- list(AB_pred_list, SRT_pred_list, SD_ME_pred_list,
                              SD_int_pred_list, CC_ME_pred_list, 
                              CC_int_pred_list)
names(predictions_all_tasks) <- c("AB", "SRT", "SD_ME", "SD_int", "CC_ME", "CC_int")
save(predictions_all_tasks, file="../data/predicted_error_by_task.RData")
