#########################################################################
# run moments models on effect size errors for all N
# K. Garner 
# see plot_pfx_vs_dstfx.R for working out that lead to this
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

print("running AB models")
run_models(f=run_AB_models, datpath=datpath, task="AB", zipnm="EPSAB_wv.zip",
           maxN=313, Ns=sub.Ns)

print("running SRT models")
run_models(f=run_SRT_models, datpath=datpath, task="SRT", zipnm="SRT_wv.zip",
           maxN=313, Ns=sub.Ns)

print("running SD ME models")
run_models(f=run_SD_ME, datpath=datpath, task="SD", zipnm="SD_wv.zip",
           maxN=313, Ns=sub.Ns)

print("running SD int models")
run_models(f=run_SD_int, datpath=datpath, task="SD", zipnm="SD_wv.zip",
           maxN=313, Ns=sub.Ns)

quit()
