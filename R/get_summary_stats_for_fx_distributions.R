## K. Garner 2022
##  quick sample of viability of comparing variability ratios
###############################################################
### assumes library(tidyverse) has been loaded in the .Rmd

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
    #dat <- dat %>% pivot_wider(id_cols = perm, names_from = "mod", values_from = "esz")
  } else if (task == "SD"){
    dat$mod[dat$mod == "RM-AN"] <- 'tt'
    dat$mod[dat$mod == "LME"] <- 'ta*tt'
    #dat <- dat %>% pivot_wider(id_cols = perm, names_from = "mod", values_from = "esz")
  }
  dat$task = task
  dat
} 

############ EXTRACT DIRECTORY ##############
get_summary_data_fx_p <- function(delete_files_after){
  fnms_fx <- "../data/%s/EPS%s.zip"
  #extfnms_fx <- unzip(fnms_fx)
  Ns <- round(exp(seq(log(13), log(313), length.out = 20)))
  tasks <- c("AB", "SD", "SRT", "CC")
  
  get_dat_acrs_tsks <- function(Ns, task, fnms_fx){
    f <- sprintf(fnms_fx, task, task)
    if(file.exists(sprintf("EPS%s", task))){
      extfnms_fx <- list.files(sprintf("EPS%s", task), full.names = TRUE)
    } else {
      extfnms_fx <- unzip(f)
    }
    do.call(rbind, lapply(Ns, get_fx_dat4_n, task=task, extfnms_fx = extfnms_fx))
  }
  dat <- do.call(rbind, lapply(tasks, get_dat_acrs_tsks, Ns=Ns, fnms_fx=fnms_fx))
  
  
  if (delete_files_after) {
    unlink(c("EPSAB", "EPSCC",
             "EPSSD", "EPSSRT"), 
           recursive=TRUE)
  }
  
  ########## GET SUMMARY STATS ################
  dat %>% select(n, p, esz, mod, task) %>%
    pivot_longer(p:esz, names_to = "meas", values_to = "value") %>%
    group_by(n, mod, task, meas) %>% 
    summarise(mu=mean(value),
              med=median(value),
              sd=sd(value),
              LB=quantile(value, .025),
              UB=quantile(value, .975))
}