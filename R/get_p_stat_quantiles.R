## written by K. Garner (2022)
## get quantiles for p-distribution from each N
######################################################################
rm(list=ls())

############ LOAD PACKAGES ##################
library(tidyverse)
library(gridExtra)

############ EXTRACT DIRECTORY ##############
fnms_fx <- "../data/CC/EPSCC.zip"
extfnms_fx <- unzip(fnms_fx)

########### DEFINE VARIABLES ################
Ns <- round(exp(seq(log(13), log(313), length.out = 20)))
task <- "CC"
ftmplt = "../data/%s/%s_pquants.RData"

##########  FUNCTIONS #####################
get_dat_fx <- function(extfnms_fx, x){
  # get one dataframe given list of names and single idx
  # extfnms [char, list] list of extracted files from directory
  # x [int] idx
  load(extfnms_fx[x])
  out
}

get_fx_dat4_n <- function(task, N, extfnms_fx){
  # for one task, and N, get the data and make it tidy
  # task [str] - which task are you extracting data from?
  # N [int] - which sample size are you extracting?
  # extfnms_fx [list] list of extracted filenames
  
  tmp <- sprintf("imm_%s_N-%d_", task, N)
  fnms2opn <- grep(tmp, extfnms_fx)
  dat <- lapply(fnms2opn, get_dat_fx, extfnms=extfnms_fx)
  dat <- do.call(rbind, dat)
  nrow <- nrow(dat)
  nmod <- length(unique(dat$mod))
  dat$perm <- rep(1:(nrow/nmod), each = nmod)
  dat <- dat[!is.na(dat$p),]
  dat
} 

############################################################################## 
###### get all the data from all the Ns
dat <- lapply(Ns, get_fx_dat4_n, task=task, extfnms_fx=extfnms_fx)
dat <- do.call(rbind, dat)


# convert p values
dat$p <- qnorm(dat$p)
dat <- dat %>% filter(is.finite(p))
############################################################################## 
###### get the stats on the max sample size

sum <- dat %>% 
    group_by(n, mod) %>%
    summarise(m=mean(p),
              s=sd(p),
              lower=quantile(p, .025),
              upper=quantile(p, .975))


############################################################################## 
###### save the output
save(sum, file= sprintf(ftmplt, task, task))
