#########################################################################
# plot fx size/p vs dist measures
# K. Garner 
# for the supfx project, plot (1) the probability of the effect size observation
# vs the distributional stats for that observation (e.g. mu, skew and 
# kurtosis of Xdiff) (2) the probability of the effect size observation and the 
# individual distributional stats (mean var, skew and kurtosis across participants
# in that sample)

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
source("efilids_functions.R") # custom functions written for this project
source("moments_functions.R")
source("R_rainclouds.R") # functions for plotting


################################################################################
# Run things
################################################################################
not_new <- T
# define sub Ns
datpath <- "../data"
tasks <- c("SD", "AB", "SRT")
zipnms <- c("SD_wv.zip", "EPSAB_wv.zip", "SRT_wv.zip")
sub.Ns <- round(exp(seq(log(13), log(313), length.out = 20)))
maxN = 313
tstN = 25
Ns <- c(tstN, maxN) # for initial testing

# first get a list, where each element is the data from one task
if (!not_new){ 
dat <- mapply(function(x,y, Ns) do.call(rbind, lapply(Ns, 
                                                  get_dat_4corrs, 
                                                  datpath=datpath, 
                                                  task=x, 
                                                  zipnm=y)),
                x=tasks,
                y = zipnms,
                MoreArgs = list(Ns = Ns))
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


#############################################################
# per task, perform 1)
# a stepwise linear regression to get vif and residuals
#############################################################
# AB
#############################################################
AB_model <- do_stp_n_prs(df = dat[["AB"]] %>% filter(n == tstN) %>%
                                     select(esz_dist, AB.skew,
                                            AB.k, AB.r, sigma_mu, k_mu),
                              dv = "esz_dist",
                              task = "AB",
                              prs_fnm = "AB",
                              rsd_fm = "AB_init")
AB_model[[1]] # all good
# gonna keep sigma mu (correlates highly with skew_mu and k_mu)
# initial model VIF
# AB.skew     AB.k     AB.r sigma_mu  skew_mu     k_mu
# 1.745574 1.602303 1.236187 5.184596 8.904735 3.439406
# AB.skew     AB.k     AB.r sigma_mu  skew_mu 
# 1.726115 1.601608 1.117864 4.569866 4.266436 
# AB.skew     AB.k     AB.r sigma_mu     k_mu # keeping k_mu as VIF lower
# 1.733866 1.601263 1.096272 1.785786 1.647888 
# at this point, residuals look good!
summary(AB_model[[2]]) 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept) -0.225060   0.018058 -12.463  < 2e-16 ***
#   AB.mu        1.689350   0.040806  41.400  < 2e-16 ***
#   AB.sigma    -6.038839   0.173207 -34.865  < 2e-16 ***
#   AB.r        -0.027334   0.009056  -3.018  0.00261 **
#   sigma_mu    -0.002633   0.001090  -2.415  0.01592 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.02954 on 995 degrees of freedom
# Multiple R-squared:  0.7725,	Adjusted R-squared:  0.7716
# F-statistic: 844.5 on 4 and 995 DF,  p-value: < 2.2e-16

#############################################################
# SRT
#############################################################
# SRT is a little more complex, got some work to do on the variables
# first
# 1) get EV of distribution
get_ev_non_norm(dat[["SRT"]]$esz_dist) #[1] 0.05454278
# 2) allocate data to a df and make variables behave
SRT_df <- dat[["SRT"]] %>% filter(n == tstN) %>%
                            select(esz_dist, SRT.mu, SRT.sigma, SRT.skew,
                                    SRT.k, SRT.r, sigma_mu, skew_mu, k_mu)
ggpairs(SRT_df)
hist(SRT_df$esz_dist, breaks = 30) # positively skewed
apply(SRT_df %>% select(SRT.k, sigma_mu, skew_mu, k_mu), 2, min) # all +ve 
# first make all the values positive and log transform
x = 0.01 - min(SRT_df$esz_dist)
SRT_df <- SRT_df %>% mutate(esz_dist_t = log(esz_dist + x),
                            SRT.r_t = scale(r2z(SRT.r), scale=F),
                            SRT.sigma_t = log(SRT.sigma), # assumes pwr function
                            SRT.k_t = log(SRT.k), # same
                            sigma_mu_t = log(sigma_mu), # same
                            skew_mu_t = log(skew_mu),
                            k_mu_t = log(k_mu)) # same

hist(SRT_df$esz_dist_t, breaks = 30) # this has squashed in the upper tail, some 
# large negative values
ggpairs(SRT_df %>% select(esz_dist_t, SRT.sigma_t, SRT.skew,
                          SRT.k_t, SRT.r_t, sigma_mu_t, k_mu_t))
# clearly skewy variables are SRT sigma, SRT.k, sigma_mu, skew_mu

#SRT_df <- SRT_df %>% select(esz_dist_t, SRT.sigma, SRT.k, sigma_mu, skew_mu)
SRT_df <- SRT_df %>% select(esz_dist_t, SRT.skew, SRT.k_t, SRT.r_t,
                            sigma_mu_t, k_mu_t)
SRT_model <- do_stp_n_prs(df = SRT_df,
                          dv = "esz_dist_t",
                          task = "SRT",
                          prs_fnm = "SRT",
                          rsd_fm = "SRT")
SRT_model[[1]]
# SRT.skew    SRT.k_t    SRT.r_t sigma_mu_t  skew_mu_t     k_mu_t 
# 2.650863   2.394736   1.824746   1.709689   4.147821   4.335670 
# with k_mu removed
# SRT.skew    SRT.k_t    SRT.r_t sigma_mu_t  skew_mu_t 
# 2.650330   2.393024   1.748995   1.609469   1.293854 
# with skew mu removed
# SRT.skew    SRT.k_t    SRT.r_t sigma_mu_t     k_mu_t # going with this one as 
# 2.605319   2.354499   1.823450   1.706594   1.352451 # did w AB and much of a much
summary(SRT_model[[2]])
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -4.68609    0.23587 -19.867  < 2e-16 ***
#   SRT.skew    -0.07786    0.02580  -3.017  0.00261 ** 
#   SRT.k_t     -0.08978    0.04310  -2.083  0.03751 *  
#   SRT.r_t      0.50777    0.03797  13.371  < 2e-16 ***
#   sigma_mu_t  -1.00907    0.04927 -20.481  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3468 on 995 degrees of freedom
# Multiple R-squared:  0.3415,	Adjusted R-squared:  0.3389 
# F-statistic:   129 on 4 and 995 DF,  p-value: < 2.2e-16

#############################################################
# SD
#############################################################
# first look at data for ME and int effects, transform variables
# to make them behave well
############################################################
# SD ME
############################################################
SD_ME_df <- dat[["SD"]] %>% filter(n == tstN) %>%
              select(esz_ME_dist, ME_skew,
                     ME_k, ME_r, sigma_mu, skew_mu, k_mu)
ggpairs(SD_ME_df)
apply(SD_ME_df, 2, min) 
# k and R are clearly non-normal
SD_ME_mod_df <- SD_ME_df %>% select(esz_ME_dist, ME_skew,
                                    ME_k, ME_r, sigma_mu, skew_mu, k_mu) %>%
                              mutate(ME_k_t = log(ME_k),
                                     ME_r_t = r2z(ME_r)) %>%
                              select(esz_ME_dist, ME_skew,
                                     ME_k_t, ME_r_t, sigma_mu, k_mu)
ggpairs(SD_ME_mod_df)
SD_ME_mod <- do_stp_n_prs(df = SD_ME_mod_df,
                          dv = "esz_ME_dist",
                          task = "SD",
                          prs_fnm = "SD_ME",
                          rsd_fm = "SD_ME")
SD_ME_mod[[1]]
# ME_skew   ME_k_t   ME_r_t sigma_mu  skew_mu     k_mu 
# 1.607146 1.566770 1.105848 1.247206 4.198354 4.421814 
# ME_skew   ME_k_t   ME_r_t sigma_mu  skew_mu #### without k
# 1.606887 1.566577 1.103301 1.160747 1.068549 
# ME_skew   ME_k_t   ME_r_t sigma_mu     k_mu 
# 1.596962 1.562803 1.096520 1.230073 1.125424 #### without skew
summary(SD_ME_mod[[2]])
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.098640 -0.018380  0.002502  0.019298  0.067126 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.063152   0.014221  -4.441 9.97e-06 ***
#   ME_skew     -0.022371   0.002217 -10.093  < 2e-16 ***
#   ME_k_t       0.020647   0.003615   5.712 1.48e-08 ***
#   ME_r_t       0.026333   0.004021   6.549 9.29e-11 ***
#   sigma_mu     0.207436   0.141691   1.464    0.144    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.02838 on 995 degrees of freedom
# Multiple R-squared:  0.1398,	Adjusted R-squared:  0.1363 
# F-statistic: 40.42 on 4 and 995 DF,  p-value: < 2.2e-16
############################################################
# SD INT
############################################################
SD_int_df <- dat[["SD"]] %>% filter(n == tstN) %>%
                    select(esz_int_dist, int_skew,
                            int_k, int_r, sigma_mu, skew_mu, k_mu)
ggpairs(SD_int_df)
# esz_int_dist and int_k are both highly skewed
x <- 0.01 - min(SD_int_df$esz_int_dist)
apply(SD_int_df %>% select(int_k, k_mu), 2, mean)
SD_int_mod_df <- SD_int_df %>%  mutate(esz_int_dist_t = log(esz_int_dist + x),
                                       int_k_t = log(int_k),
                                       k_mu_t = log(k_mu)) %>%
                                select(esz_int_dist_t, int_skew, int_k_t,
                                       int_r, sigma_mu, k_mu_t)
SD_int_mod <- do_stp_n_prs(df = SD_int_mod_df,
                           dv = "esz_int_dist_t",
                           task = "SD",
                           prs_fnm = "SD_int",
                           rsd_fm = "SD_int")
SD_int_mod[[1]]
# int_skew  int_k_t    int_r sigma_mu   k_mu_t 
# 1.205284 1.231812 1.028055 1.107811 1.112723 
summary(SD_int_mod[[2]])
# lm(formula = esz_int_dist_t ~ int_skew + int_k_t + int_r + sigma_mu, 
#    data = df)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.8350 -0.4520  0.1674  0.6112  1.8183 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.66239    0.36145  -1.833  0.06716 .  
# int_skew      0.16763    0.05934   2.825  0.00482 ** 
#   int_k_t       0.56474    0.10666   5.295 1.47e-07 ***
#   int_r         0.91010    0.11847   7.682 3.73e-14 ***
#   sigma_mu    -21.30041    4.04767  -5.262 1.74e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8502 on 995 degrees of freedom
# Multiple R-squared:  0.1158,	Adjusted R-squared:  0.1123 
# F-statistic: 32.58 on 4 and 995 DF,  p-value: < 2.2e-16


### result is that for all tasks, skew, kurtosis, correlation, and sigma_mu
### predict error in effect size estimate
### now I need to implement each process for each task and each N,
### collecting the models in a list