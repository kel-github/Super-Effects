# binding together VSL data
# this script was written to put data back together after a bug was found in the application
# of the t-test data (the bug was that I used a tidyverse function to group by subjects from 
# a sample without renumbering the subject numbers first, meaning that only unique numbers
# were being retained).
# this function assumes that the .tar.gz data archive has been converted to a zip and placed
# in out/
# out/ should also contain out.zip - the data that has the t-test output only

# there should be a folder in data/VSL called VSL
# this code takes the existing data from the VSL folder (data/VSL) and matches the filenames
# to those from the new datasets (out/), loads and selects the correct data from both,
# before binding to replace the files in data/VSL - ready for zipping and downloading.
# this was run manually on Mon 8th Nov, 2021

rm(list=ls())

# -----------------------------------------------------------------
# load libraries and functions
# -----------------------------------------------------------------
library(tidyverse)
source('efilids_functions.R')

# -----------------------------------------------------------------
# define session variables
# -----------------------------------------------------------------
save_dir <- '../data/VSL/'
tmp_dir <- '../data/VSL/out/'
task <- "VSL"
subfol <- "VSL"
sub_Ns <- round(exp(seq(log(13), log(313), length.out = 20)))
rxvnme <- "VSL"

# -----------------------------------------------------------------
# relatively constant settings
# -----------------------------------------------------------------
fstem <- "_N-%d_parent-%d.RData"
N <- sub_Ns
j <- 1:1000
j <- j[-c(66,119,152)]

datpath <- "../data/"
fname_add <- NULL # string or NULL

# -----------------------------------------------------------------
# get list of files from each folder
# -----------------------------------------------------------------
fprev <- unzip(paste(tmp_dir, task, ".zip", sep=""), list=TRUE)
ft <- unzip(paste(tmp_dir, "out", ".zip", sep=""), list=TRUE)
fprev <- str_extract(fprev$Name, "([^/]+$)")
ft <- str_extract(ft$Name, "([^/]+$)")
fs <- print(intersect(fprev, ft)) # this is our list of files we'll load and unzip from both archives
fs <- fs[!is.na(fs)]

# -----------------------------------------------------------------
# for each string in a, is it in b? if so load both, combine and 
# save, if not print out name of string as a warning message
# -----------------------------------------------------------------

data_binds <- function(fs, rxv_a, rxv_b, tmp_dir, save_dir){ # this will be applied over list a

  unzip(zipfile=rxv_a, files = paste(task, fs, sep="/"), exdir = paste(tmp_dir, "tmp", sep=""))
  unzip(zipfile=rxv_b, files = paste("out", task, fs, sep = "/"), exdir = paste(tmp_dir, "tmp", sep=""))
  load(paste(tmp_dir, "tmp/", "VSL/", fs, sep=""))
  o <- out %>% filter(mod == "LME")
  rm(out)
  load(paste(tmp_dir, "tmp/", "out/", "VSL/", fs, sep = ""))
  n <- out %>% filter(mod == "RM-AN")
  rm(out)
  out <- rbind(o, n)
  save(out, file=paste(save_dir, task, "/", fs, sep=""))
}

rxv_a <- paste(tmp_dir, task, ".zip", sep = "")
rxv_b <- paste(tmp_dir, "out", ".zip", sep = "")

lapply(fs, data_binds, rxv_a = rxv_a, rxv_b = rxv_b, tmp_dir = tmp_dir, save_dir = save_dir)

# -----------------------------------------------------------------
# manually zip data rxivs and move to appropriate places
# -----------------------------------------------------------------
