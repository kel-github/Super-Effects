## K. Garner 2022
## getting bits of information attained at each N for each task
###############################################################
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
fx <- "b*c" # this is only necessary when selecting conditions
bns <- c(15, 20, 25) # bin values to test

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
    dat <- dat %>% pivot_wider(id_cols = perm, names_from = "mod", values_from = "esz")
  } else if (task == "SD"){
    dat$mod[dat$mod == "RM-AN"] <- 'tt'
    dat$mod[dat$mod == "LME"] <- 'ta*tt'
    dat <- dat %>% pivot_wider(id_cols = perm, names_from = "mod", values_from = "esz")
  }
  dat
} 

compute_entropy <- function(esz, bn, min_esz, max_esz){
  # compute the entropy for a given vector esz
  # KWARGS
  # -- esz [vector] a vector of values, over which entropy will be computed
  # -- bn [int]: how many bins for the data?
  # -- min_esz [int]: lower bound value from which to begin computing bins
  # -- max_esz [int]: as above, but max
  # RETURNS
  # -- H = -sum(p*log(p))
  
  # use the hist function to do the binning and counring
  h <- hist(esz, breaks = seq(min_esz, max_esz, length.out=bn+1)) 
  counts <- h$counts
  counts[counts == 0] <- .00001/length(counts[counts==0]) # hacky way to get rid of the zeros
  p <- counts/sum(counts)
  -sum(p*log(p))
}


compute_entropy_N <- function(bn, min_esz, max_esz, task, extfnms_fx, N, fx = NA){
  # use the get_data and compute_entropy functions to load one N and get 
  # the entropy for that distribution
  # KWARGS
  # -- esz [vector] a vector of values, over which entropy will be computed
  # -- bn [int]: how many bins for the data?
  # -- min_esz [int]: lower bound value from which to begin computing bins
  # -- max_esz [int]: as above, but max
  # -- fx [str] - select a condition for models that have multiple fx
  #                 e.g. "tt" for SD
  #                
  # RETURNS
  # -- H = -sum(p*log(p))
  
  dat <- get_fx_dat4_n(task, N, extfnms_fx)
  if (task == "SD"|task == "CC"){
    esz <- dat[,fx][[1]]
  } else {
    esz <- dat$esz
  }
  compute_entropy(esz, bn, min_esz, max_esz)
}

######## RUNNING! #####################
######## first there is a manual section to check binning and range works ######## 
######## for min and max N ######################################################## 
tmp_min <- get_fx_dat4_n(task, Ns[1], extfnms_fx)

tmp_min$esz <- tmp_min[,fx][[1]] #### this line is for dealing with models with ME and interactions

min_esz <- min(tmp_min$esz) #- .9
max_esz <- max(tmp_min$esz) #+ .9
# get some reasonable bin values for this one

hist(tmp_min$esz, breaks = seq(min_esz, max_esz, length.out=bns[1]+1))
hist(tmp_min$esz, breaks = seq(min_esz, max_esz, length.out=bns[2]+1))
hist(tmp_min$esz, breaks = seq(min_esz, max_esz, length.out=bns[3]+1))

rm(tmp_min)
tmp_max <- get_fx_dat4_n(task, Ns[length(Ns)], extfnms_fx)

tmp_max$esz <- tmp_max[,fx][[1]] # as above

hist(tmp_max$esz, breaks = seq(min_esz, max_esz, length.out=bns[1]+1))
hist(tmp_max$esz, breaks = seq(min_esz, max_esz, length.out=bns[2]+1))
hist(tmp_max$esz, breaks = seq(min_esz, max_esz, length.out=bns[3]+1))

# AB bns[3] looks best
AB_bn <- 3
SRT_bn <- 3 
SD_ME <- 3
SD_int <- 3
CC_ME <- 3
CC_int <- 3

###### RUN AB ###### 
AB_H <- unlist(lapply(Ns, compute_entropy_N, bn=bns[AB_bn], min_esz=min_esz, max_esz=max_esz,
                      task=task, extfnms_fx=extfnms_fx))
save(AB_H, file="../data/AB/AB_H.RData")

###### RUN SRT ###### 
SRT_H <- unlist(lapply(Ns, compute_entropy_N, bn=bns[SRT_bn], min_esz=min_esz, max_esz=max_esz,
                      task=task, extfnms_fx=extfnms_fx))
save(SRT_H, file="../data/SRT/SRT_H.RData")

###### RUN SD ME ###### 
SD_ME_H <- unlist(lapply(Ns, compute_entropy_N, bn=bns[SD_ME], min_esz=min_esz, max_esz=max_esz,
                       task=task, extfnms_fx=extfnms_fx, fx="tt"))
save(SD_ME_H, file="../data/SD/SD_ME_H.RData")

###### RUN SD INT ###### 
SD_INT_H <- unlist(lapply(Ns, compute_entropy_N, bn=bns[SD_int], min_esz=min_esz, max_esz=max_esz,
                          task=task, extfnms_fx=extfnms_fx, fx="ta*tt"))
save(SD_INT_H, file="../data/SD/SD_INT_H.RData")

###### RUN CC ME ########
CC_ME_H <- unlist(lapply(Ns, compute_entropy_N, bn=bns[CC_ME], min_esz=min_esz, max_esz=max_esz,
                          task=task, extfnms_fx=extfnms_fx, fx="c"))
save(CC_ME_H, file="../data/CC/CC_ME_H.RData")

###### RUN CC INT ########
CC_INT_H <- unlist(lapply(Ns, compute_entropy_N, bn=bns[CC_int], min_esz=min_esz, max_esz=max_esz,
                          task=task, extfnms_fx=extfnms_fx, fx="b*c"))
save(CC_INT_H, file="../data/CC/CC_INT_H.RData")
