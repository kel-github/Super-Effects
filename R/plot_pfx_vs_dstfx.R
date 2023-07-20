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
source("R_rainclouds.R") # functions for plotting

##############################################################################
## functions
############################################################################
get_dat_4corrs <- function(datpath, task, zipnm, N){
  # for a given N, compile fx size and distributional data
  # datpath [str] - where is the data - e.g. "../data"
  # task [str] - which task are you getting data for?
  # zipnm [str] - what is the name of the zipped folder?
  # N [int] - which sample size are we working on?
  
  # first, list the files in the zip rxv
  tmp = unzip(zipfile=paste(datpath, task, zipnm, sep="/"), list = TRUE)
  tmp = tmp$Name[!is.na(str_extract(tmp$Name, paste("N-", N, "_", sep = "")))]
  dst_idx <- do.call(cbind, lapply(tmp, str_detect, pattern = "diststats"))
  fnms <- tmp[!dst_idx]
  # now get the fx data
  fx_dat <- do.call(rbind, lapply(fnms, unzp_fx_n_spew, datpath=datpath, task=task, zipnm=zipnm))
  # and reshape/relabel it if it has 2 effects in it
  if(task == "SD"){
    fx_dat$mod[fx_dat$mod == "RM-AN"] = "ME"
    fx_dat$mod[fx_dat$mod == "LME"] = "int"
    fx_dat <- fx_dat %>% pivot_wider(id_expand = FALSE, names_from = mod, values_from = c(esz, p))
  }
  # and the mu diff and var stats
  dst_nms <- tmp[dst_idx]
  dst_dat <- do.call(rbind, lapply(dst_nms, unzp_mu_stats_n_spew, 
                                   datpath=datpath,
                                   task=task,
                                   zipnm=zipnm))
  
  # join up fx_dat and dst_dat and return
  inner_join(fx_dat, dst_dat, by = "parent")
}


unzp_fx_n_spew <- function(datpath, task, zipnm, fnm){
  # given a list of filenames, unextract
  # each one, open, and add to dataframe
  # output the key variables from the 
  # dataframe
  f <- unzip(paste(datpath, task, zipnm, sep="/"), 
                file=fnm,
                exdir=paste(datpath, task, sep="/"))
  load(f)
  row_idx <- !is.na(out$p)
  out <- out[row_idx, ]
  # add parent #
  pN <- strsplit(fnm, split = "-")
  pNidx <- !is.na(sapply(pN, str_extract, pattern = ".RData"))
  nustr <- pN[[1]][pNidx]
  parent <- as.numeric(strsplit(nustr, split="R")[[1]][1])
  out$parent <- parent
  
  # will need to pivot wider here for tasks with 2 fx
  
  # output
  out %>% select(n, p, esz, mod, parent)
}

unzp_mu_stats_n_spew <- function(datpath, task, zipnm, fnm){
  # given a list of filenames, unextract
  # each one, open, and add to dataframe
  # output the key variables from the 
  # dataframe
  f <- unzip(paste(datpath, task, zipnm, sep="/"), 
             file=fnm,
             exdir=paste(datpath, task, sep="/"))
  load(f)
  
  if (task == "AB" | task == "SRT") {
    dst_dat <- cbind(dist_out[[1]][1], dist_out[[2]])
  } else if (task == "SD"){
    dst_dat <- cbind(dist_out[[1]]$ME, dist_out[[1]]$int, dist_out[[2]])
    names(dst_dat) <- c("ME_mu", "ME_sigma", "ME_skew", "ME_k", "ME_r",
                        "int_mu", "int_sigma", "int_skew", "int_k", "int_r",
                        "sigma_mu", "skew_mu", "k_mu")
  } 
  # now get parent number
  
  # add parent #
  pN <- strsplit(fnm, split = "-")
  pNidx <- !is.na(sapply(pN, str_extract, pattern = "_diststats"))
  nustr <- pN[[1]][pNidx]
  parent <- as.numeric(strsplit(nustr, split="_")[[1]][1])
  dst_dat$parent <- parent
  
  # output
  dst_dat
}

get_ev_non_norm <- function(x){
  # given a random variable x, compute the expected value
  # note, check with collaborators they're happy with this
  d <- density(x)
  mu <- sum(d$x*d$y)*diff(d$x)[1]
  mu
}

do_stp_n_prs <- function(df, dv, task, datpath = "../data/", prs_fnm, rsd_fm){
  # given the dataframe df, perform an AIC stepwise
  # regression process, with dv as the dv
  # dv [str], corresponding to a named variable in df
  # the remaining columns in df will serve as the predictor 
  # variables
  # task [str] - which task are you performing analysis on?
  # datpath [str] - where is the parent data folder?
  # prs_fnm [str] - what base do you want the plot fnms to have?
  # rsd_fnm [str] - same as above but for residuals
  # first do pairs plot and save
  prs <- ggpairs(df)
  ggsave(paste(prs_fnm, "pairs.pdf", sep="_"), plot=prs, 
         path=paste(datpath, task, sep="/"),
         width=30,
         height=30,
         units = "cm")
  
  # perform linear regression and test the VIF factor
  full_model <- lm(as.formula(paste(dv, "~ .")), data=df)
  vif_out <- vif(full_model)
  step_model <- stepAIC(full_model, direction = "both", 
                        trace = FALSE)
  pdf(paste(datpath, task, paste(rsd_fm, "resid.pdf", sep="_"), sep="/"))
  plot(step_model, which = c(1, 2))
  dev.off()
  
  pdf(paste(datpath, task, paste(rsd_fm, "aV.pdf", sep="_"), sep="/"))
  avPlots(step_model)
  dev.off()

  list(vif_out, step_model)
}

r2z <- function(r){
  # Fisher R to Z transform https://en.wikipedia.org/wiki/Fisher_transformation
  0.5*log((1+r)/(1-r))
}
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
