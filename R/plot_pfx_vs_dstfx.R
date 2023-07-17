#########################################################################
# plot fx size/p vs dist measures
# K. Garner 
# for the supfx project, plot (1) the probability of the effect size observation
# vs the distributional stats for that observation (e.g. mu, skew and 
# kurtosis of Xdiff) (2) the probability of the effect size observation and the 
# individual distributional stats (mean var, skew and kurtosis across participants
# in that sample)

rm(list=ls())

library(tidyverse) # for data wrangling
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(lme4) # for mixed effects modelling
library(ggridges)
library(car)
library(parallel)
library(rstatix)
library(moments)
library(stringr)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting

# define sub Ns
datpath <- "../data"
task <- "SD"
zipnm <- "EPSSD_wv.zip"
sub.Ns <- round(exp(seq(log(13), log(313), length.out = 20)))

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
  
  if (task == "AB") {
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

Ns <- c(21, 313)

dat <- do.call(rbind, lapply(Ns, 
                             get_dat_4corrs, 
                             datpath=datpath, 
                             task=task, 
                             zipnm=zipnm))

# plan per task
# convert each effect size to a probability (integrate)
# for the median N, get the data for each task and perform 
# a linear regression to see what accounts for p, after controlling
# for mu and sigma of Xdiff
# perform the linear regression at each level of N, and take the
# beta co-efficients for the key features for each task
# do a correlation between N and beta, as well as plotting (each
# feature is a panel, each line is a task)
# what general/conclusive statements can be made?


### plotting
me_sml <- dat %>% filter(n == 21) %>% select("esz_ME", "ME_mu", "ME_sigma",
                                             "ME_skew", "ME_k", "ME_r", "sigma_mu",
                                             "skew_mu", "k_mu")
me_lrg <- dat %>% filter(n == 313) %>% select("esz_ME", "ME_mu", "ME_sigma",
                                             "ME_skew", "ME_k", "ME_r", "sigma_mu",
                                             "skew_mu", "k_mu")
int_sml <- dat %>% filter(n == 21) %>% select("esz_int", "int_mu", "int_sigma",
                                              "int_skew", "int_k", "int_r", "sigma_mu",
                                              "skew_mu", "k_mu")
int_lrg <- dat %>% filter(n == 313) %>% select("esz_int", "int_mu", "int_sigma",
                                              "int_skew", "int_k", "int_r", "sigma_mu",
                                              "skew_mu", "k_mu")
## compute the density function over each set of effect sizes, and get the 
# probability of each value
pairs(me_sml)
pairs(me_lrg)

pairs(int_sml)
pairs(int_lrg)


with(me_sml, summary(lm(esz_ME ~ ME_mu + ME_sigma + ME_skew + ME_k + ME_r + sigma_mu + skew_mu + k_mu)))
with(me_lrg, summary(lm(esz_ME ~  ME_mu + ME_sigma + ME_skew + ME_k + ME_r + sigma_mu + skew_mu + k_mu)))
with(int_sml, summary(lm(esz_int ~ int_mu + int_sigma + int_skew + int_k + int_r + sigma_mu + skew_mu + k_mu)))
with(int_lrg, summary(lm(esz_int ~ int_mu + int_sigma + int_skew + int_k + int_r + sigma_mu + skew_mu + k_mu)))


# 
# with(sml, summary(lm(AB.sigma ~ AB.skew + AB.k + AB.r + sigma_mu + skew_mu + k_mu)))
# with(lrg, summary(lm(AB.mu ~ AB.skew + AB.k + AB.r + sigma_mu + skew_mu + k_mu)))
# with(lrg, summary(lm(AB.sigma ~ AB.skew + AB.k + AB.r + sigma_mu + skew_mu + k_mu)))
