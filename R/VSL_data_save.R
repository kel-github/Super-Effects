# binding together VSL data
# this script was written to put data back together after a bug was found in the application
# of the t-test data (the bug was that I used a tidyverse function to group by subjects from 
# a sample without renumbering the subject numbers first, meaning that only unique numbers
# were being retained).
# this code takes the existing data from the VSL folder (data/VSL) and matches the filenames
# to those from the new datasets (out/), loads and selects the correct data from both,
# before binding to replace the files in data/VSL.
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
data_dir <- '../data/VSL/'
tmp_dir <- '../data/out/'
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
# unzip files of interest
# -----------------------------------------------------------------
lapply(sub_Ns, unzp, datpath = data_dir, rxvnme = paste(rxvnme, ".zip", sep=""), rxvsub = rxvnme, task = task, j = j)

# -----------------------------------------------------------------
# grab a list of files from both directories 
# -----------------------------------------------------------------
original <- list.files(paste(data_dir, task, sep = ""))
nu <- list.files(paste(tmp_dir, task, sep=""))

# -----------------------------------------------------------------
# for each string in a, is it in b? if so load both, combine and 
# save, if not print out name of string as a warning message
# -----------------------------------------------------------------

data_binds <- function(a, list_b){ # this will be applied over list a
  tryCatch({
    if (as.logical(sum(str_detect(a, list_b)))){
      load(paste(data_dir, task, "/", a, sep=""))
      o <- out %>% filter(mod == "LME")
      rm(out)
      load(paste(tmp_dir, task, "/", a, sep = ""))
      n <- out %>% filter(mod == "RM-AN")
      rm(out)
      out <- rbind(o, n)
      save(out, file=paste(data_dir, task, "/", a, sep=""))
    }
  },
  warning = function(cond){
    message(paste(a))
  })
}

lapply(original, data_binds, list_b = nu)

# -----------------------------------------------------------------
# check which from the new list are not in the old list, and
# move those datafiles to out
# -----------------------------------------------------------------
missing <- setdiff(nu, original)

deal_w_missing_files <- function(m, orig_dir, nu_dir){
  
  file.copy(from = paste(orig_dir, m, sep = ""), to = paste(nu_dir, m, sep = ""))
  file.remove(paste(orig_dir, m, sep = ""))  
}
lapply(missing, deal_w_missing_files, orig_dir = paste(data_dir, "VSL/", sep=""), nu_dir = paste(tmp_dir, "VSL/", sep=""))

# -----------------------------------------------------------------
# manually zip data rxivs and move to appropriate places
# -----------------------------------------------------------------
