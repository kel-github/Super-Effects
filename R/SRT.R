### written by K. Garner, April 2020
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
library(tidyverse) # for data wrangling
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(lme4) # for mixed effects modelling
library(ggridges)
library(car)
library(parallel)
library(rstatix)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting

args <- commandArgs(trailingOnly=TRUE)
Nind <- NA
if (length(args) == 0) {
  fname <- "../data/total_of_313_subs_SRT_task_trial_level_data.csv"
  outpath <- "$HOME/tmp"
} else if (length(args) == 1) {
  fname <- args[1]
  outpath <- "$HOME/tmp"
} else if (length(args) >= 2) {
  fname <- args[1]
  outpath <- args[2]
}
if (length(args) == 3) {
  Nind <- as.integer(args[3])
}

set.seed(42) 
# -----------------------------------------------------------------------------
# load data and wrangle into tidy form
# (see https://r4ds.had.co.nz/tidy-data.html), plus relabel to make
# labels a little simpler
# -----------------------------------------------------------------------------
dat <- read.csv(fname,
                header=TRUE)

# ----------------------------------------------------------------------------------------------------
# Create dataframes 
# ----------------------------------------------------------------------------------------------------

# Create a summary of the data for fixed fx modelling
min.RT <- 200 # in msec
sd.crit <- 2.5

ffx.dat <- dat %>% filter(Block.No > 2) %>%
              group_by(Subj.No, Block.No.Names) %>%
              filter(Accuracy == 1) %>%
              filter(RT.ms > min.RT) %>%
              filter(RT.ms < (mean(RT.ms) + sd.crit*sd(RT.ms))) %>%
              summarise(RT=mean(RT.ms))


# ----------------------------------------------------------------------------------------------------
# define levels for simulations
# ----------------------------------------------------------------------------------------------------

sub.Ns <- round(exp(seq(log(13), log(313), length.out = 20)))
if (!is.na(Nind)) {
  sub.Ns <- sub.Ns[Nind]
}
n.perms =1000# for each sample size, we will repeat our experiment n.perms times
cores = 20
subs  <- unique(ffx.dat$Subj.No)

# ----------------------------------------------------------------------------------------------------
# run simulations, getting p values from linear models, and cohen's d values, and save results to a list, using intermediate sampling
# ----------------------------------------------------------------------------------------------------
fstem <- paste(outpath, "/imm_SRT_N-%d_parent-%d.RData", sep="")
lapply(sub.Ns, function(x) run.outer(in.data=ffx.dat, subs=subs, N=x, 
                                     k=1, j=n.perms, cores=cores, 
                                     f=get.ps.srt, fstem=fstem, samp="imm"))

# ----------------------------------------------------------------------------------------------------
# run simulations, getting p values from linear models, and cohen's d values, and save results to a list, using intermediate sampling
# ----------------------------------------------------------------------------------------------------
fstem <- paste(outpath, "/SRT_N-%d_parent-%d.RData", sep="")
lapply(sub.Ns, function(x) run.outer(in.data=ffx.dat, subs=subs, N=x, 
                                     k=n.perms, j=n.perms, cores=cores, 
                                     f=get.ps.srt, fstem=fstem, samp="int"))

# ----------------------------------------------------------------------------------------------------
# get outta here
# ----------------------------------------------------------------------------------------------------
quit()
