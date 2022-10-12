## written by K. Garner (2022)
## get percentage of values that fall within 1 and 1.96 z-scores of the 
## best estimate, for each task x N 
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
ftmplt = "../data/%s/%s_phitbest.RData"

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
  if (task == "CC") {
    dat$mod[dat$mod == "RM-AN"] <- 'b*c'
    dat$mod[dat$mod == "LME"] <- 'c'
 #   dat <- dat %>% pivot_wider(id_cols = perm, names_from = "mod", values_from = "esz")
  } else if (task == "SD"){
    dat$mod[dat$mod == "RM-AN"] <- 'tt'
    dat$mod[dat$mod == "LME"] <- 'ta*tt'
#    dat <- dat %>% pivot_wider(id_cols = perm, names_from = "mod", values_from = "esz")
  }
  dat
} 

############################################################################## 
###### get all the data from all the Ns
dat <- lapply(Ns, get_fx_dat4_n, task=task, extfnms_fx=extfnms_fx)
dat <- do.call(rbind, dat)

############################################################################## 
###### get the stats on the max sample size
if (task == "SD" | task == "CC"){
  sum <- dat %>% filter(n==313) %>%
                 group_by(mod) %>%
                 summarise(m=mean(esz),
                           s=sd(esz),
                           lower1=quantile(esz, .25),
                           upper1=quantile(esz, .75),
                           lower2=quantile(esz, .025),
                           upper2=quantile(esz, .975))
} else {
  sum <- dat %>% filter(n==313) %>%
                 summarise(m=mean(esz),
                           s=sd(esz),
                           lower1=quantile(esz, .25),
                           upper1=quantile(esz, .75),
                           lower2=quantile(esz, .025),
                           upper2=quantile(esz, .975))
}

if (task == "SD" | task == "CC"){
  dat <- inner_join(dat, sum, by="mod")
} else {
  dat <- cbind(dat, sum)
}

if (task == "SD" | task == "CC"){
  out <- dat %>% mutate(sdev1 = ifelse(esz > lower1 & esz < upper1, 1, 0),
                        sdev2 = ifelse(esz > lower2 & esz < upper2, 1, 0)) %>%
                  group_by(n, mod) %>%
                  summarise(p_lower = sum(sdev1)/length(sdev1),
                            p_upper = sum(sdev2)/length(sdev2))
} else {
  out <- dat %>% mutate(sdev1 = ifelse(esz > lower1 & esz < upper1, 1, 0),
                        sdev2 = ifelse(esz > lower2 & esz < upper2, 1, 0)) %>%
                        group_by(n) %>%
                        summarise(p_lower = sum(sdev1)/length(sdev1),
                                  p_upper = sum(sdev2)/length(sdev2))
}

############################################################################## 
###### save the output
save(out, file= sprintf(ftmplt, task, task))

