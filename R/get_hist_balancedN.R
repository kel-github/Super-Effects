## K. Garner 2022
## plot histograms of balanced N (more studies with small N, less with big N)
## for each paradigm
###############################################################
rm(list=ls())

# balancing
# 18 x 50, 25 x 36, 30 x 30 , 36 x 25, 50 x 18.

############ LOAD PACKAGES ##################
library(tidyverse)
library(gridExtra)

########### DEFINE VARIABLES ################
Ns <- c(18, 25, 30, 36, 50)
samps <- c(50, 36, 30, 25, 18)
tasks <- c("AB", "SD", "SRT", "CC")


##########  FUNCTIONS #####################
get_dat_fx <- function(extfnms_fx, x){
  # get one dataframe given list of names and single idx
  # extfnms [char, list] list of extracted files from directory
  # x [int] idx
  load(extfnms_fx[x])
  out
}

get_fx_dat4_n <- function(task, N, samp, extfnms_fx){
  # for one task, and N, get the data and make it tidy
  # task [str] - which task are you extracting data from?
  # N [int] - which sample size are you extracting?
  # samp [int] - how many samples of the N do you want?
  # extfnms_fx [list] list of extracted filenames
  
  tmp <- sprintf("imm_%s_N-%d_", task, N)
  fnms2opn <- grep(tmp, extfnms_fx)
  fnms2opn <- sample(fnms2opn, samp) # replace = FALSE is the default
  dat <- lapply(fnms2opn, get_dat_fx, extfnms=extfnms_fx)
  dat <- do.call(rbind, dat)
  nrow <- nrow(dat)
  nmod <- length(unique(dat$mod))
  dat$perm <- rep(1:(nrow/nmod), each = nmod)
  dat <- dat[!is.na(dat$p),]
  if (task == "CC") {
    dat$mod[dat$mod == "RM-AN"] <- 'b*c'
    dat$mod[dat$mod == "LME"] <- 'c'
  } else if (task == "SD"){
    dat$mod[dat$mod == "RM-AN"] <- 'tt'
    dat$mod[dat$mod == "LME"] <- 'ta*tt'
  }
  dat$n <- N
  dat$samp <- samp
  dat$task <- task
  dat %>% select(n, esz, mod, samp, task)
} 

get_fx_dat_4N_4task <- function(task, Ns, samps){
  
  ############ EXTRACT DIRECTORY ##############
  fnms_fx <- sprintf("../data/%s/EPS%s.zip", task, task)
  extfnms_fx <- unzip(fnms_fx)
  ############ GET DATA ##############
  do.call(rbind, mapply(get_fx_dat4_n, N=Ns, samp=samps, 
                        MoreArgs=list(task=task, 
                                      extfnms_fx=extfnms_fx),
                        SIMPLIFY = FALSE))
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
  
  # use the hist function to do the binning and counting
  h <- hist(esz, breaks = seq(min_esz, max_esz, length.out=bn+1)) 
  counts <- h$counts
  counts[counts == 0] <- .00001/length(counts[counts==0]) # hacky way to get rid of the zeros
  p <- counts/sum(counts)
  -sum(p*log(p))
}

######## RUNNING! #####################
baldat <- do.call(rbind, lapply(tasks, function(x) get_fx_dat_4N_4task(task=x, 
                                                             Ns=Ns, 
                                                             samps=samps)))

###### just tidy up baldat ############
baldat$mod[baldat$task == "AB"] <- "AB"
baldat$mod[baldat$task == "SRT"] <- "SRT"

###### now plot ############
baldat %>% ggplot(aes(x=esz, group=as.factor(n), colour=as.factor(n))) +
           geom_histogram(alpha=0.2) +
           facet_wrap(~mod, scales="free")


get_entropy_per_mod <- function(Ns, task, min_esz, max_esz){
      e <- unlist(lapply(Ns, function(x) compute_entropy(baldat$esz[baldat$n == x & baldat$mod == task],
                                           bn=30,
                                           min_esz=min_esz,
                                           max_esz=max_esz)))
      tibble(e=e, task=task, n=Ns)
}

ent <- do.call(rbind, mapply(get_entropy_per_mod, task=unique(baldat$mod),
                                           min_esz=c(0.3, 0.7, -0.1,  0, -0.1, -0.1),
                                           max_esz=c(0.9, 1,  1, 5,  0.3, 0.7),
                                           MoreArgs=list(Ns=Ns),
                                           SIMPLIFY=FALSE))

ent %>% ggplot(aes(x=n, y=e, group=task, colour=task)) + geom_line()
