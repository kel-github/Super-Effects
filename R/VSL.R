### written by K. Garner, April 2020
### edited by Z. Nott, August 2020
### for the project 'On the detectability of effects in executive
### function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

rm(list=ls())
### run analysis of sample size x effect size variability on the SRT data
# -----------------------------------------------------------------------------
# load packages and source function files
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(lme4) # for mixed effects modelling
library(ggridges)
library(parallel)
library(tidyverse) # for data wrangling
source("efilids_functions.R") # custom functions written for this project
#source("R_rainclouds.R") # functions for plotting

args <- commandArgs(trailingOnly=TRUE)

n.inner <- 1000
n.outer <- 1000
i.outer <- NA
cores <- 1
sub.Ns <- round(exp(seq(log(13), log(313), length.out = 20)))

if (length(args) == 0) {
  fname <- "../data/total_of_313_subs_VSL_task_trial_level_data.csv"
  outpath <- "../out/VSL"
} else if (length(args) == 1) {
  fname <- args[1]
  outpath <- "../out/VSL"
} else if (length(args) >= 2) {
  fname <- args[1]
  outpath <- args[2]
}
if (length(args) == 3) {
  i.outer <- as.integer(args[3])
}
set.seed(42)
seeds <- sample(1:n.outer, n.outer, replace=FALSE)

# -----------------------------------------------------------------------------
# load data and wrangle into tidy form
# (see https://r4ds.had.co.nz/tidy-data.html), plus relabel to make
# labels a little simpler
# -----------------------------------------------------------------------------
dat <- read.csv(fname, header=TRUE)

# -----------------------------------------------------------------------------
# Create dataframe for analysis
# -----------------------------------------------------------------------------

# data frame contains TRUE ordering
prev.dat <- dat %>% select(Subj.No, Trial.No, Response, Target.Order, Accuracy)
prev.dat$Response <- as.factor(prev.dat$Response)
prev.dat$Target.Order <- as.factor(prev.dat$Target.Order)

prev.dat <- prev.dat %>% mutate(Response = recode(Response,
                                      "122" = "Novel",
                                      "109" = "Repeat"),
                                Target.Order = recode(Target.Order,
                                      "1" = "Novel",
                                      "2" = "Repeat"))

# -----------------------------------------------------------------------------
# run simulations for t-test model, getting p values from t.tests, and
# cohen's d values, and save results to a list
# -----------------------------------------------------------------------------

fstem <- paste(outpath, "/VSL_N-%d_parent-%d.RData", sep="")
subs  <- unique(prev.dat$Subj.No)
start  <-  Sys.time()
lapply(sub.Ns, function(x) run.outer(in.data=prev.dat, subs=subs, N=x,
                                     k=n.inner, j=n.outer, outer_index=i.outer,
                                     cores=cores,
                                     f=get.ps.vsl,
                                     fstem=fstem,
                                     samp="int",
                                     seeds=seeds))
end <-  Sys.time()
end - start

