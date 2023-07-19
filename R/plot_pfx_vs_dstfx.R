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
  list(vif_out, step_model)
}

################################################################################
# Run things
################################################################################

# define sub Ns
datpath <- "../data"
tasks <- c("SD", "AB", "SRT")
zipnms <- c("SD_wv.zip", "EPSAB_wv.zip", "SRT_wv.zip")
sub.Ns <- round(exp(seq(log(13), log(313), length.out = 20)))
maxN = 313
tstN = 25
Ns <- c(tstN, maxN) # for initial testing

# first get a list, where each element is the data from one task
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
AB_init_model <- do_stp_n_prs(df = dat[["AB"]] %>% filter(n == tstN) %>%
                                     select(esz_dist, AB.mu, AB.sigma, AB.skew,
                                            AB.k, AB.r, sigma_mu),
                              dv = "esz_dist",
                              task = "AB",
                              prs_fnm = "AB",
                              rsd_fm = "AB_init")
AB_init_model[[1]] # sigma and skew_mu are an issue, but residuals looks ok
# gonna keep sigma mu (correlates highly with skew_mu and k_mu)
# initial model VIF
# AB.mu  AB.sigma   AB.skew      AB.k      AB.r  sigma_mu   skew_mu      k_mu 
# 2.360360  2.252898  2.173843  1.976478  1.669528  7.360513 11.926566  5.276604
# AB.mu AB.sigma  AB.skew     AB.k     AB.r sigma_mu 
# 2.140984 1.465127 1.929311 1.866657 1.185319 1.939769
# at this point, residuals look good!
summary(AB_init_model[[2]]) 
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

# for mu and sigma of Xdiff
# apply this over tasks
#############
#############################################################
# SRT
#############################################################
SRT_init_model <- do_stp_n_prs(df = dat[["SRT"]] %>% filter(n == tstN) %>%
                                      select(esz_dist, SRT.mu, SRT.sigma, SRT.skew,
                                             SRT.k, SRT.r, sigma_mu, skew_mu),
                                    dv = "esz_dist",
                                    task = "SRT",
                                    prs_fnm = "SRT",
                                    rsd_fm = "SRT")
SRT_init_model[[1]]
# skew mu and k_mu a problem. 
# SRT.mu SRT.sigma  SRT.skew     SRT.k     SRT.r  sigma_mu   skew_mu      k_mu 
# 1.766901  2.639953  2.839966  2.581512  2.039296  2.182733  7.050590  6.510980 
# correlation skew_mu & k_mu = 0.914, ditching k_mu
# second model
#    SRT.mu SRT.sigma  SRT.skew     SRT.k     SRT.r  sigma_mu   skew_mu 
# 1.753563  2.616690  2.834695  2.568372  2.039295  2.168606  1.296429 
# residuals show we're mis-fitting high and low values, will take a 
# look at all histograms
dat[["SRT"]] %>% filter(n == tstN) %>% 
                  select(esz_dist, SRT.mu, SRT.sigma, SRT.skew,
                          SRT.k, SRT.r, sigma_mu, skew_mu) %>%
                  pivot_longer(cols = everything(), names_to = "meas", values_to = "values") %>%
                    ggplot(aes(x=values)) + geom_histogram() + facet_wrap(~meas, scales = "free")
# sigma_mu, skew_mu, SRT.k and SRT.sigma are positively skewed
# SRT.r is negatively skewed
# decide re: transforms - see if that will help
# https://rcompanion.org/handbook/I_12.html

# perform the linear regression at each level of N, and take the
# beta co-efficients for the key features for each task
# do a correlation between N and beta, as well as plotting (each
# feature is a panel, each line is a task)
# what general/conclusive statements can be made?

