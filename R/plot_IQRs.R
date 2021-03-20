### written by K. Garner, March 2021
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

# ----------------------------------------------------------------------------------------------------
rm(list=ls())
# ----------------------------------------------------------------------------------------------------
### ATTAIN DENSITIES FOR ANY TASK GIVEN RANGE INPUTS
# ----------------------------------------------------------------------------------------------------
# load packages and source function files
# ----------------------------------------------------------------------------------------------------

rm(list=ls())

# ----------------------------------------------------------------------------------------------------
# load packages and source function files
# ----------------------------------------------------------------------------------------------------

library(tidyverse) # for data wrangling
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(lme4) # for mixed effects modelling
library(ggridges)
library(car)
library(parallel)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R")
# ----------------------------------------------------------------------------------------------------
# define session variables
# ----------------------------------------------------------------------------------------------------

task = "AB"
d_scale_ffx = 2
d_scale_rfx = 4
p_scale_ffx = 2
p_scale_rfx = 2
px_rng_d = c(0,3)
px_rng_p_ffx = c(-800,0)
px_rng_p_rfx = c(-800,0)
width = 8
height = 8

# ----------------------------------------------------------------------------------------------------
# define datas and load d's
# ----------------------------------------------------------------------------------------------------

fnames = list(paste("../data/", task, "_d", "_d.RData", sep=""), paste("../data/", task, "_int_1_d", "_d.RData", sep=""))

get.dat <- function(fname){
  load(fname)
  out <- d
  out
}
dat <- lapply(fnames, get.dat)

# ----------------------------------------------------------------------------------------------------
# define subs for plotting and get data into one dataframe
# ----------------------------------------------------------------------------------------------------
subs <- unique(dat[[1]]$Nsz)
subs.2.plt <- c(13, 25, 50, 115)

relabel <- function(dat, samps, x){
  data <- dat[[x]]
  data$samp = samps[x]
  data
}

samps <- c("A", "B")
dat <- lapply(c(1:length(samps)), relabel, dat=dat, samps=samps)
dat <- do.call(rbind, dat)
dat <- dat %>% filter(Nsz %in% subs.2.plt)

# ----------------------------------------------------------------------------------------------------
# now plot
# ----------------------------------------------------------------------------------------------------


