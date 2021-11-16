rm(list = ls())
# -----------------------------------------------------------------
# load packages and source function files
# -----------------------------------------------------------------

library(tidyverse) # for data wrangling
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(lme4) # for mixed effects modelling
library(ggridges)
library(car)
library(parallel)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting
# -----------------------------------------------------------------
# define session variables
# -----------------------------------------------------------------
task <- "VSL"
subfol <- "VSL"
sub_Ns <- round(exp(seq(log(13), log(313), length.out = 20)))
convert <- NA
rxvnme <- "VSL"
svnme <- 'tmp'

fstem <- "_N-%d_parent-%d.RData"
N <- sub_Ns
j <- 1000

datpath <- "../data/"
fname_add <- NULL # string or NULL

# -----------------------------------------------------------------
# define functions
# -----------------------------------------------------------------
get_data <- function(f, datpath, rxvnme, svnme){
  unzip(paste(datpath, rxvnme, "/", rxvnme, ".zip", sep=""),
        files = f, exdir = paste(datpath, svnme, sep = ""))
  load(paste(datpath, svnme, "/", f, sep=""))
  out
}

sv_data <- function(f, sv_path, rxvnme, dat){
  out = dat
  save(out, file=paste(datpath, sv_path, "/", f, sep=""))
}

data_proc <- function(f, datpath, rxvnme, svnme, k=1000, alpha = 0.05) {
  # see get.data for input arg info
  # : -- convert = model name for which conversion is required
  
  # do data preprocessing, largely for passing 
  # dat into the plotting or the stats function
  dat <- get_data(f, datpath, rxvnme, svnme)
  
  idx <- dat$mod == "LME"
  n <- dat$n[1]
  
  # recalc ps
  p <- dat[idx, "p"]
  p <- ((p*k)+1)/(k+1)
  
  # recalc gamma
  esz <- dat[idx, "esz"]
  esz = (alpha^(1/n) - p^(1/n)) / (1 - p^(1/n))
  # if gamma_zero = - inf, it means that all the perms were greater than the min accuracy
  # thus, the prevalence should be 0 in this case
  esz[p > alpha] = 0
  
  # now put back
  dat$esz[idx] <- esz
  dat$p[idx] <- p
  
  sv_data(f, svnme, rxvnme, dat)
}

# -----------------------------------------------------------------
# get list of files to extract
# -----------------------------------------------------------------
fs <- unzip(paste(datpath, "VSL/VSL.zip", sep=""), list = TRUE)
fs <- fs$Name[2:length(fs$Name)] # this will become 2:end
lapply(fs, data_proc, datpath = datpath, rxvnme = rxvnme, svnme = svnme)

