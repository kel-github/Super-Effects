### written by K. Garner, April 2020
### edited by Z. Nott, July 2020
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nolan, C., Nydam, A, Nott, Z., Bowman, H., & Dux, PE 

rm(list=ls())
### run analysis of sample size x effect size variability on the SRT data
# ----------------------------------------------------------------------------------------------------
# load packages and source function files
library(tidyverse) # for data wrangling
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(lme4) # for mixed effects modelling
library(ggridges)
library(car)
library(parallel)
library(rstatix)
library(moments)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting

args <- commandArgs(trailingOnly=TRUE)

n.outer <- 1000
n.inner <- 1000
i.outer <- NA
cores <- 10
sub.Ns <- round(exp(seq(log(13), log(313), length.out = 20)))

if (length(args) == 0) {
  fname <- "../data/total_of_313_subs_SingDual_task_trial_level_data.csv"
  outpath <- "../data/SD"
} else if (length(args) == 1) {
  fname <- args[1]
  outpath <- "$HOME/tmp"
} else if (length(args) >= 2) {
  fname <- args[1]
  outpath <- args[2]
}
if (length(args) == 3) {
  Nind <- as.integer(args[3])
}

set.seed(42) 
seeds <- sample(1:n.outer, n.outer, replace=FALSE)

# ----------------------------------------------------------------------------------------------------
# load data and wrangle into tidy form (see https://r4ds.had.co.nz/tidy-data.html), plus relabel to make
# labels a little simpler
# ----------------------------------------------------------------------------------------------------
dat <- read.csv(fname,
                header=TRUE)

# ----------------------------------------------------------------------------------------------------
# Create dataframe
# ----------------------------------------------------------------------------------------------------
# Create a summary of the data for ffx and rfx modelling
min.RT <- .200 # in sec
sd.crit <- 2.5

rfx.dat <- dat %>% filter(Overall.Accuracy == 1) %>%
                   select(c('Subj.No', 'Trial.Type.Name', 'Task.1.RT.Sec', 'Task.2.RT.Sec', 'Task.1.Response', 'Task.2.Response'))  %>%
                   pivot_longer(c('Task.1.RT.Sec', 'Task.2.RT.Sec'), names_to = "task", values_to="RT") %>%
                   drop_na()
rfx.dat$task[rfx.dat$Trial.Type.Name == 'single_auditory'] = 'sound'
rfx.dat$task[rfx.dat$Trial.Type.Name == 'single_visual'] = 'vis'
rfx.dat$task[rfx.dat$Trial.Type.Name == 'dual_task' & rfx.dat$task == 'Task.1.RT.Sec'] = 'vis'
rfx.dat$task[rfx.dat$Trial.Type.Name == 'dual_task' & rfx.dat$task == 'Task.2.RT.Sec'] = 'sound'

rfx.dat <- rfx.dat %>% mutate(trialtype = fct_recode(Trial.Type.Name,
                                                                    'single' = 'single_auditory',
                                                                    'single' = 'single_visual',
                                                                    'dual' = 'dual_task')) 
rfx.dat$task.stim <- NA
rfx.dat$task.stim[rfx.dat$trialtype == "single"] = rfx.dat$Task.1.Response[rfx.dat$trialtype == "single"]
rfx.dat$task.stim[rfx.dat$trialtype == "dual" & rfx.dat$task == "vis"] = rfx.dat$Task.1.Response[rfx.dat$trialtype == "dual" & rfx.dat$task == "vis"]
rfx.dat$task.stim[rfx.dat$trialtype == "dual" & rfx.dat$task == "sound"] = rfx.dat$Task.2.Response[rfx.dat$trialtype == "dual" & rfx.dat$task == "sound"]

ffx.dat <- rfx.dat %>% select(-c("Trial.Type.Name", "Task.1.Response", "Task.2.Response")) %>%
                        group_by(Subj.No, task, trialtype) %>%
                        filter(RT > min.RT) %>%
                        filter(RT < (mean(RT)+sd.crit*sd(RT))) %>%
                        summarise(RT = mean(RT))

sub_var <- rfx.dat %>% select(-c("Trial.Type.Name", "Task.1.Response", "Task.2.Response")) %>%
                        group_by(Subj.No, task, trialtype) %>%
                        filter(RT > min.RT) %>%
                        filter(RT < (mean(RT)+sd.crit*sd(RT))) %>%
                        group_by(Subj.No) %>%
                        summarise(sigma_sq = var(RT),
                                  skew = skewness(RT),
                                  k = kurtosis(RT))
write_csv(sub_var, file = "../data/SD/SD_sub_var_stats.csv")
subs  <- unique(ffx.dat$Subj.No)

ffx.dat <- inner_join(ffx.dat, sub_var, by = "Subj.No")

# ----------------------------------------------------------------------------------------------------
# run simulations, getting p values from linear models, and cohen's d values, and save results to a list, using immediate sampling
# ----------------------------------------------------------------------------------------------------
fstem <- paste(outpath, "/imm_SD_N-%d_parent-%d.RData", sep="")
lapply(sub.Ns, function(x) run.outer(in.data=ffx.dat, subs=subs, N=x, k=1,
                                     j=n.outer, outer_index=i.outer,
                                     cores=cores, 
                                     f=get.ps.SD, 
                                     fstem=fstem, 
                                     samp="imm",
                                     seeds=seeds))

# # ----------------------------------------------------------------------------------------------------
# get outta here
# ----------------------------------------------------------------------------------------------------
quit()
