## written by K. Garner (2022)
## plot percentage of values that fall within 1 and 1.96 z-scores of the 
## best estimate, for each task x N 
######################################################################

########### DEFINE VARIABLES ################
Ns <- round(exp(seq(log(13), log(313), length.out = 20)))
tasks <- c("AB", "CC", "SRT", "SD")
ftmplt = "../data/%s/%s_phitbest.RData"

########### LOAD DATASETS ################
get_prct_dats <- function(ftmplt, task){
  load(sprintf(ftmplt, task, task))
  if(task == "AB" | task == "SRT"){ 
    out$mod <- NA
    out$mod <- task
  }
  out
}

plt_dat <- do.call(rbind, lapply(tasks, get_prct_dats, ftmplt=ftmplt))

