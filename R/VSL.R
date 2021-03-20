### written by K. Garner, April 2020
### edited by Z. Nott, August 2020
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

rm(list=ls())
### run analysis of sample size x effect size variability on the SRT data
# ----------------------------------------------------------------------------------------------------
# load packages and source function files

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the location of this file
# uncomment the below and run if you need to install the packages
# install.packages("tidyverse")
# install.packages("wesanderson")
# install.packages("cowplot")
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(lme4) # for mixed effects modelling
library(ggridges)
library(parallel)
library(tidyverse) # for data wrangling
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting

set.seed(42) # testing diff seeds on output
# ----------------------------------------------------------------------------------------------------
# load data and wrangle into tidy form (see https://r4ds.had.co.nz/tidy-data.html), plus relabel to make
# labels a little simpler
# ----------------------------------------------------------------------------------------------------
dat = read.csv("../data/total_of_313_subs_VSL_task_trial_level_data.csv", header=TRUE)

# ----------------------------------------------------------------------------------------------------
# Create dataframe for analysis
# ----------------------------------------------------------------------------------------------------

# data frame contains TRUE ordering
prev.dat <- dat %>% select(Subj.No, Trial.No, Response, Target.Order, Accuracy) 
prev.dat$Response <- as.factor(prev.dat$Response)
prev.dat$Target.Order <- as.factor(prev.dat$Target.Order)

prev.dat <- prev.dat %>% mutate(Response = recode(Response,
                                      "122" = "Novel",
                                      "109" = "Repeat"),
                                Target.Order = recode(Target.Order,
                                      "1" = "Novel",
                                      "2" = "Repeat"))

# ----------------------------------------------------------------------------------------------------
# define levels for simulations
# ----------------------------------------------------------------------------------------------------

sub.Ns = round(exp(seq(log(13), log(313), length.out = 20)))
sub.Ns = 313
n.perms = 1000# for each sample size, we will repeat our experiment n.perms times
k = 2 #outer loop
Np = 1000
cores = 2
# 13, 10, 10, 10 = 4.384 mins
# 13, 15, 15, 15 = 10.075 
# 13, 20, 20, 20 = 14.606 
# 313, 10, 10, 10 = 1.227 hrs
# 313, 1000, 1, 1000 = 
# to run 
# ----------------------------------------------------------------------------------------------------
# run simulations for t-test model, getting p values from t.tests, and cohen's d values, and save results to a list
# ----------------------------------------------------------------------------------------------------

subs  <- unique(prev.dat$Subj.No)
start  <-  Sys.time()

lapply(sub.Ns, function(x) run.outer(in.data=prev.dat, subs=subs, N=x, k=n.perms, j=k, cores=cores, ffx.f=run.os.t.test.sim, rfx.f=run.prev.test, fstem="VSL_N-%d_parent-%d.RData"))
end <-  Sys.time()
end - start

# ----------------------------------------------------------------------------------------------------
# save the data of import to an RData file
# ----------------------------------------------------------------------------------------------------
save(sims.dat, prev.res, file="VS_sim_data.RData")