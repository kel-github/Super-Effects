# this code runs the plots for the effect size and p-value distributions
# this code also calculates the stats and saves to a list
# both run info and functions are contained in the code
# written by K. Garner, 2021

### Call plotting functions for observed effect sizes and p values
rm(list = ls())
# -----------------------------------------------------------------
# load packages and source function files
# -----------------------------------------------------------------

library(tidyverse) # for data wrangling
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(lme4) # for mixed effects modelling
library(ggridges)
library(car)
library(parallel)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting

# -----------------------------------------------------------------
# task settings
# -----------------------------------------------------------------
# AB
# convert = LME

# CC
# convert = LME

# SRT & VSL
# have both coded in

# SD
# convert = LME

# VSL
# convert RM-AN (prev = % of pop)
# -----------------------------------------------------------------
# define session variables
# -----------------------------------------------------------------
task <- "SRT"
subfol <- "SRT"
sub_Ns <- round(exp(seq(log(13), log(313), length.out = 20)))
convert <- NA
rxvnme <- "SRT"

# -----------------------------------------------------------------
# relatively constant settings
# -----------------------------------------------------------------
fstem <- "_N-%d_parent-%d.RData"
N <- sub_Ns
j <- 1000
if (task == "VSL"){
  j <- 1:1000
# #  j <- j[-66]
  j <- j[-c(66,119,152)]
}

datpath <- "../data/"
fname_add <- NULL #"IMM" #NULL # string or NULL

# -----------------------------------------------------------------
# define functions
# -----------------------------------------------------------------
get_data <- function(fstem, n, j, datpath, rxvnme) {
  # concatenate data for one task
  # -- fstem: fstem = filestem to be sprintf'd with N and j
  # -- n: vector of sub sample sizes
  # -- j: outer loop size
  # -- dv: which dv do you want to know about (d or p)?
  # -- datpath: where is the data? (relative path)
  # -- rxvnme: name of zipped rxv folder e.g. "CC"
  # -- model: "rfx" (LME) or "ffx" (e.g. ANOVA)
  if (nchar(rxvnme) <= nchar("VSL")) {
       dn <- lapply(n, function(x) unzp(paste(datpath, rxvnme, "/", sep = ""),
                                   paste(rxvnme, ".zip", sep = ""),
                                   rxvnme,
                                   rxvnme,
                                   j,
                                   x))
  } else {
      dn <- lapply(n, function(x)
                       unzp(paste(datpath, sub("IMM", "", rxvnme), "/", sep = ""),
                            paste(rxvnme, ".zip", sep = ""),
                            rxvnme,
                            paste("imm_", sub("IMM", "", rxvnme), sep = ""),
                            j,
                            x))
  }
  dn <- do.call(rbind, dn)

  get.dat <- function(f) {
    load(f)
    out
  }
  ds <- do.call(rbind, lapply(dn, get.dat))
  ds
}

d2r <- function(dat, m) {
  # apply d2r transform for data in dat corresponding to model m
  dat %>% filter(mod == m) %>%
    group_by(n) %>%
    mutate(esz = (esz / (sqrt((esz^2) + 4)))^2) %>% # Cohen 1988 equation 2.2.6
    ungroup()
}

data_proc <- function(fstem, n, j, datpath, rxvnme, convert) {
  # see get.data for input arg info
  # : -- convert = model name for which conversion is required

  # do data preprocessing, largely for passing 
  # dat into the plotting or the stats function
  dat <- get_data(fstem, n, j, datpath, rxvnme)
  
  if (!is.na(convert)){
  dat <- rbind(dat %>% filter(mod != convert),
               d2r(dat, convert)) 
  }
  
  dat <- dat %>% mutate(p = qnorm(p))
  dat %>% ungroup()
}

compute_stats <- function(dat) {
   # given the output of the data preprocessing (data.proc), 
   # compute the required stats
# -----------------------------------------------------------------
   # for each model produce a density for
   # the effect sizes and for the p-values
# -----------------------------------------------------------------
   do_dens <- function(dat, x, dv) {
     tryCatch(
      {
        density(dat[dat$mod == x, dv])
      },
        error = function(cond) {
        NULL
      }
     )
   }
   mods <- unique(dat$mod)
   dens_fx <- lapply(mods, do_dens, dat = dat, dv = "esz")
   names(dens_fx) <- mods
   dens_p <- lapply(mods, do_dens, dat = dat, dv = "p")
   names(dens_p) <- mods

# -----------------------------------------------------------------
   # for each model, get the central tendency, sd,
   # and .025 & .975 quantiles
# -----------------------------------------------------------------
   get_stats <- function(y) {
     mu <- mean(y)
     med <- median(y)
     sd <- sd(y)
     qs <- quantile(y, probs = c(.025, .975))
     list(mu=mu, med=med, sd=sd, qs=qs)
   }
   stats_fx <- sapply(mods,
              function(x) get_stats(dat$esz[dat$mod == x & is.finite(dat$esz)]))
   stats_p <- sapply(mods,
              function(x) get_stats(pnorm(dat$p[dat$mod == x & is.finite(dat$esz)])))

  # -------------------------------------------------------------
   # for each model, get the percent of 'significant results'
   # -----------------------------------------------------------
   psig <- function(data,y) {
      sum(data$p[data$mod == y] < qnorm(.05)) / length(data$p[data$mod == y])
   }
   sig <- sapply(mods, psig, data = dat)

  # -----------------------------------------------------
   # get the mean effect size from 'sig results'
   # ----------------------------------------------------
   stats_sig <- sapply(mods,
                       function(x) get_stats(dat$esz[dat$mod == x & is.finite(dat$esz) & dat$p < qnorm(.05)]) )

   #--------------------------------------------------
   # return it all!
   # -------------------------------------------------
   res <- list(dens_fx, dens_p, stats_fx, stats_p, sig, stats_sig)
   names(res) <- c("dens_fx", "dens_p", 
                   "stats_fx", "stats_p", 
                   "sig", "stats_sig")
   res
}

stats_4_subs <- function(fstem, n, j, datpath, rxvnme, convert) {
  # for a subject group, apply the functions to preprocess the data
  # compute stats
  # return the list of results
  # for use in application over each level of subject
  dat <- data_proc(fstem, n, j, datpath, rxvnme, convert)
  compute_stats(dat)
}

# ------------------------------------------------------------
# run the code across each subject group
# ------------------------------------------------------------
res <- lapply(sub_Ns, stats_4_subs,
                      fstem = fstem,
                      j = j,
                      datpath = datpath,
                      rxvnme = rxvnme,
                      convert = convert)

names(res) <- paste(sub_Ns)
res <- do.call(rbind, res) # makes it neater for reffing

# ------------------------------------------------------------
# save the output
# ------------------------------------------------------------
if (is.null(fname_add)){
  save(res, file = paste(datpath, task, "/", task, "stats.RData", sep = ""))
} else {
  save(res, file = paste(datpath, task, "/", task, fname_add, "stats.RData", sep = ""))
}
# ------------------------------------------------------------
# delete the unzipped files
# ------------------------------------------------------------
if (nchar(rxvnme) <= nchar("VSL")) {
  unlink(paste(datpath, task, "/", task, sep = ""), recursive = TRUE)
} else {
  unlink(paste(datpath, task, "/", rxvnme, sep = ""), recursive = TRUE)
}
