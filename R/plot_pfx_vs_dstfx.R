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
Ns <- list(c(42, 313), c(25, 313), c(36, 313)) 

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


#############################################################
# per task, perform 1)
# a stepwise linear regression to get vif and residuals
#############################################################
# AB

#############################################################
# Things I can say so far:
# for all tasks, skew, kurtosis, correlation between 
# variables and subject variability predicted error in effect 
# size estimates
#############################################################
tstN <- 25
AB_model <- do_stp_n_prs(df = dat[["AB"]] %>% filter(n == tstN) %>%
                                     select(esz_dist, AB.skew,
                                            AB.k, sigma_mu, k_mu),
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
# (Intercept) -0.298680   0.032354  -9.232  < 2e-16 ***
#   AB.skew     -0.051905   0.005644  -9.197  < 2e-16 ***
#   AB.k         0.013495   0.002912   4.634 4.06e-06 ***
#   sigma_mu     0.016322   0.001759   9.276  < 2e-16 ***
#   k_mu        -0.827387   0.122728  -6.742 2.65e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05018 on 995 degrees of freedom
# Multiple R-squared:  0.3436,	Adjusted R-squared:  0.3409 
# F-statistic: 130.2 on 4 and 995 DF,  p-value: < 2.2e-16
#### now get standardised residuals for winning model
tmp = dat[["AB"]] %>% filter(n == tstN) %>%
                        select(esz_dist, AB.skew,
                        AB.k, sigma_mu, k_mu)
summary(lm(scale(esz_dist) ~ scale(AB.skew) + scale(AB.k) + scale(sigma_mu) +
                         scale(k_mu), data=tmp))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      2.539e-17  2.567e-02   0.000        1    
#   scale(AB.skew)  -3.081e-01  3.350e-02  -9.197  < 2e-16 *** (2)
#   scale(AB.k)      1.504e-01  3.246e-02   4.634 4.06e-06 *** (4)
#   scale(sigma_mu)  3.145e-01  3.391e-02   9.276  < 2e-16 *** (1)
#   scale(k_mu)     -2.211e-01  3.280e-02  -6.742 2.65e-11 *** (3)
# order = 
# now take the winning formula and add group interaction terms
AB_dat <- dat[["AB"]]
AB_grp_int <- lm(esz_dist ~ AB.skew + AB.k + sigma_mu + k_mu +
                            AB.skew*n + AB.k*n + sigma_mu*n + k_mu*n ,
                        data = AB_dat) # no group interaction
plot(AB_grp_int, which=1)
summary(AB_grp_int)
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
                                    SRT.k, sigma_mu, k_mu)
ggpairs(SRT_df)
hist(SRT_df$esz_dist, breaks = 30) # positively skewed
apply(SRT_df %>% select(SRT.k, sigma_mu, k_mu), 2, min) # all +ve 
# first make all the values positive and log transform
x = 0.01 - min(SRT_df$esz_dist)
SRT_df <- SRT_df %>% mutate(esz_dist_t = log(esz_dist + x),
                            SRT.k_t = log(SRT.k), # same
                            sigma_mu_t = log(sigma_mu), # same
                            k_mu_t = log(k_mu)) # same

SRT_df <- SRT_df %>% select(esz_dist_t, SRT.skew, SRT.k_t, 
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
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   -4.68609    0.23587 -19.867  < 2e-16 ***
#   SRT.skew    -0.07786    0.02580  -3.017  0.00261 **
#   SRT.k_t     -0.08978    0.04310  -2.083  0.03751 *
#   SRT.r_t      0.50777    0.03797  13.371  < 2e-16 ***
#   sigma_mu_t  -1.00907    0.04927 -20.481  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# GET STANDARDISED
summary(lm(scale(esz_dist_t) ~ scale(SRT.skew) + scale(SRT.k_t) + scale(sigma_mu_t) + 
            scale(k_mu_t), data = SRT_df))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       -5.564e-16  2.800e-02   0.000  1.00000    
#   scale(SRT.skew)   -4.324e-01  4.221e-02 -10.244  < 2e-16 *** (1)
#   scale(SRT.k_t)     1.476e-01  4.116e-02   3.587  0.00035 *** (3)
#   scale(sigma_mu_t) -3.925e-01  3.184e-02 -12.329  < 2e-16 *** (2)
#   scale(k_mu_t)     -7.421e-02  3.082e-02  -2.408  0.01624 *   (4)


# now make group level dataframe
SRT_grp_df <- dat[["SRT"]] %>%
            select(esz_dist, n, SRT.skew,
              SRT.k, sigma_mu, k_mu)
apply(SRT_grp_df %>% select(SRT.k, sigma_mu, k_mu), 2, min) # all +ve 
# first make all the values positive and log transform
x = 0.01 - min(SRT_df$esz_dist)
SRT_grp_df <- SRT_grp_df %>% mutate(esz_dist_t = log(esz_dist + x),
                                    SRT.k_t = log(SRT.k), # same
                                    sigma_mu_t = log(sigma_mu), # same
                                    k_mu_t = log(k_mu)) # same
SRT_grp_int <- lm(esz_dist ~ SRT.skew + SRT.k_t + sigma_mu_t + k_mu_t +
                                        SRT.skew*n + SRT.k_t*n + sigma_mu_t*n + 
                                        k_mu_t*n ,
                         data = SRT_grp_df) 
# Call:
#   lm(formula = esz_dist ~ SRT.skew + SRT.k_t + sigma_mu_t + k_mu_t + 
#        SRT.skew * n + SRT.k_t * n + sigma_mu_t * n + k_mu_t * n, 
#      data = SRT_grp_df)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.79681 -0.12668 -0.01028  0.11069  1.30911 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -2.103e+00  2.114e-01  -9.948  < 2e-16 ***
#   SRT.skew     -2.421e-01  2.028e-02 -11.937  < 2e-16 ***
#   SRT.k_t       8.365e-02  3.231e-02   2.589  0.00970 ** 
#   sigma_mu_t   -5.106e-01  3.913e-02 -13.049  < 2e-16 ***
#   k_mu_t       -8.885e-02  2.746e-02  -3.235  0.00123 ** 
#   n            -2.748e-04  1.376e-03  -0.200  0.84174    
# SRT.skew:n   -5.786e-05  1.460e-04  -0.396  0.69190    
# SRT.k_t:n     6.497e-04  2.215e-04   2.933  0.00339 ** 
#   sigma_mu_t:n  2.193e-04  2.616e-04   0.838  0.40201    
# k_mu_t:n      1.075e-04  1.769e-04   0.608  0.54349    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2185 on 1990 degrees of freedom
# Multiple R-squared:  0.2415,	Adjusted R-squared:  0.2381 
# F-statistic:  70.4 on 9 and 1990 DF,  p-value: < 2.2e-16
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
                     ME_k, sigma_mu, k_mu)

apply(SD_ME_df, 2, min) 
# k and R are clearly non-normal
SD_ME_mod_df <- SD_ME_df %>% select(esz_ME_dist, ME_skew,
                                    ME_k, sigma_mu, k_mu) %>%
                              mutate(ME_k_t = log(ME_k)) %>%
                              select(esz_ME_dist, ME_skew,
                                     ME_k_t, sigma_mu, k_mu)
ggpairs(SD_ME_mod_df)
SD_ME_mod <- do_stp_n_prs(df = SD_ME_mod_df,
                          dv = "esz_ME_dist",
                          task = "SD",
                          prs_fnm = "SD_ME",
                          rsd_fm = "SD_ME")
SD_ME_mod[[1]]
# ME_skew   ME_k_t sigma_mu     k_mu 
# 1.585676 1.567407 1.179480 1.171554 
# summary(SD_ME_mod[[2]])

summary(SD_ME_mod[[2]])
# Call:
#   lm(formula = esz_ME_dist ~ ME_skew + ME_k_t, data = df)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.124190 -0.014643  0.001763  0.017530  0.060017 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.018665   0.003901  -4.785 1.97e-06 ***
#   ME_skew     -0.022481   0.002404  -9.350  < 2e-16 ***
#   ME_k_t       0.022554   0.003794   5.944 3.83e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.02505 on 997 degrees of freedom
# Multiple R-squared:  0.08079,	Adjusted R-squared:  0.07895 
# F-statistic: 43.81 on 2 and 997 DF,  p-value: < 2.2e-16
# GET STANDARDISED CO-EFFICIENTS
summary(lm(scale(esz_ME_dist) ~ scale(ME_skew) + scale(ME_k_t), data = SD_ME_mod_df))
# (Intercept)     9.072e-17  3.035e-02   0.000        1    
#   scale(ME_skew) -3.542e-01  3.788e-02  -9.350  < 2e-16 *** (1)
#   scale(ME_k_t)   2.252e-01  3.788e-02   5.944 3.83e-09 *** (2)

# now do group comparison
SD_ME_grp_df <- dat[["SD"]] %>% select(esz_ME_dist, n, ME_skew,
                                        ME_k) %>% 
                                 mutate(ME_k_t = log(ME_k)) %>%
                                 select(esz_ME_dist, n, ME_skew,
                                        ME_k_t)
SD_ME_grp_int <- lm(esz_ME_dist ~ ME_skew + ME_k_t + ME_skew*n + ME_k_t*n,
                    data = SD_ME_grp_df)
plot(SD_ME_grp_int, which=1)
# summary(SD_ME_grp_int)
# Call:
#   lm(formula = esz_ME_dist ~ ME_skew + ME_k_t + ME_skew * n + ME_k_t * 
#        n, data = SD_ME_grp_df)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.124190 -0.009785  0.001125  0.011977  0.060017 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.880e-02  3.679e-03  -5.110 3.53e-07 ***
#   ME_skew     -2.263e-02  2.269e-03  -9.973  < 2e-16 ***
#   ME_k_t       2.277e-02  3.574e-03   6.371 2.32e-10 ***
#   n            3.252e-06  2.494e-05   0.130    0.896    
# ME_skew:n    3.508e-06  1.546e-05   0.227    0.820    
# ME_k_t:n    -5.231e-06  2.388e-05  -0.219    0.827    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01979 on 1994 degrees of freedom
# Multiple R-squared:  0.08018,	Adjusted R-squared:  0.07787 
# F-statistic: 34.76 on 5 and 1994 DF,  p-value: < 2.2e-16
############################################################
# SD INT
############################################################
SD_int_df <- dat[["SD"]] %>% filter(n == tstN) %>%
                    select(esz_int_dist, int_skew,
                            int_k, sigma_mu, k_mu)
ggpairs(SD_int_df)
# esz_int_dist and int_k are both highly skewed
x <- 0.01 - min(SD_int_df$esz_int_dist)
apply(SD_int_df %>% select(int_k, k_mu), 2, min)
SD_int_mod_df <- SD_int_df %>%  mutate(esz_int_dist_t = log(esz_int_dist + x),
                                       int_k_t = log(int_k),
                                       k_mu_t = log(k_mu)) %>%
                                select(esz_int_dist_t, int_skew, int_k_t,
                                       sigma_mu, k_mu_t)
SD_int_mod <- do_stp_n_prs(df = SD_int_mod_df,
                           dv = "esz_int_dist_t",
                           task = "SD",
                           prs_fnm = "SD_int",
                           rsd_fm = "SD_int")
SD_int_mod[[1]]
# int_skew  int_k_t    int_r sigma_mu   k_mu_t 
# 1.205284 1.231812 1.028055 1.107811 1.112723 
# summary(SD_int_mod[[2]])
# Call:
#   lm(formula = esz_int_dist_t ~ int_skew + int_k_t + sigma_mu + 
#        k_mu_t, data = df)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.8281 -0.3510  0.1543  0.5120  1.5610 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -2.73986    0.97678  -2.805  0.00513 ** 
#   int_skew      0.20335    0.06625   3.069  0.00220 ** 
#   int_k_t       0.87406    0.11512   7.592 7.22e-14 ***
#   sigma_mu    -22.25101    4.52420  -4.918 1.02e-06 ***
#   k_mu_t        1.60878    0.67740   2.375  0.01774 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7563 on 995 degrees of freedom
# Multiple R-squared:  0.1089,	Adjusted R-squared:  0.1053 
# F-statistic:  30.4 on 4 and 995 DF,  p-value: < 2.2e-16
# GET STANDARDIZED
summary(lm(scale(esz_int_dist_t) ~ scale(int_skew) + scale(int_k_t) + 
                                     scale(sigma_mu) +
                                     scale(k_mu_t), data=SD_int_mod_df)) 
#   (Intercept)      4.930e-17  2.991e-02   0.000   1.0000    
#   scale(int_skew)  1.058e-01  3.446e-02   3.069   0.0022 **  (3)
#   scale(int_k_t)   2.650e-01  3.490e-02   7.592 7.22e-14 *** (1)
#   scale(sigma_mu) -1.591e-01  3.236e-02  -4.918 1.02e-06 *** (2)
#   scale(k_mu_t)    7.687e-02  3.237e-02   2.375   0.0177 *   (4)


# group comp:
SD_int_grp_df <- dat[["SD"]] %>% select(esz_int_dist, n, int_skew,
                                          int_k, sigma_mu, k_mu) %>%  
                                mutate(esz_int_dist_t = log(esz_int_dist + x),
                                          int_k_t = log(int_k),
                                          k_mu_t = log(k_mu)) %>%
                               select(esz_int_dist_t, n, int_skew, int_k_t,
                                          sigma_mu, k_mu_t)
SD_int_grp_int <- lm(esz_int_dist_t ~ int_skew + int_k_t + sigma_mu + k_mu_t +
                                    int_skew*n + int_k_t*n + sigma_mu*n +
                                    k_mu_t*n, data=SD_int_grp_df)
plot(SD_int_grp_int, which=1)
summary(SD_int_grp_int)

# Call:
#   lm(formula = esz_int_dist_t ~ int_skew + int_k_t + sigma_mu + 
#        k_mu_t + int_skew * n + int_k_t * n + sigma_mu * n + k_mu_t * 
#        n, data = SD_int_grp_df)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.82815 -0.24897  0.06122  0.33855  1.56096 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2.633e+00  8.984e-01  -2.931  0.00342 ** 
#   int_skew     1.791e-01  6.082e-02   2.945  0.00327 ** 
#   int_k_t      8.357e-01  1.058e-01   7.897 4.67e-15 ***
#   sigma_mu    -2.419e+01  4.155e+00  -5.822 6.76e-09 ***
#   k_mu_t       1.671e+00  6.232e-01   2.681  0.00740 ** 
#   n           -2.542e-03  6.145e-03  -0.414  0.67913    
# int_skew:n   5.778e-04  4.066e-04   1.421  0.15544    
# int_k_t:n    9.140e-04  7.182e-04   1.273  0.20332    
# sigma_mu:n   4.613e-02  2.791e-02   1.653  0.09847 .  
# k_mu_t:n    -1.477e-03  4.275e-03  -0.346  0.72970    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5822 on 1990 degrees of freedom
# Multiple R-squared:  0.1353,	Adjusted R-squared:  0.1314 
# F-statistic: 34.59 on 9 and 1990 DF,  p-value: < 2.2e-16
