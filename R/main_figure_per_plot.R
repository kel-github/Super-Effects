# # use this code to compile all the subplots for
# the behaviour and fx size into one 6 panel plot

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
source("plotting.R") # functions for plotting density and p value info
source("behaviour_plots.R") # functions for plotting task specific behaviours
source("fig_label.R")
# ----------------------------------------------------
# LISTS OF SETTINGS
# ----------------------------------------------------
# AB med = 24
# CC med = 23
# SRT med = 39
# SD med = 24
# ----------------------------------------------------
# CONSTANT SETTINGS
# ----------------------------------------------------
sub_Ns <- paste(round(exp(seq(log(13), log(313), length.out = 20))))
w = 1.96 * 3 # width of the plot, in inches
h = 2.36 * 2 # height
task <- "CC"

# ----------------------------------------------------
# behavioural data
# ----------------------------------------------------
fname <- "../data/total_of_313_subs_CC_task_trial_level_data.csv"

# ----------------------------------------------------
# density variables
# ----------------------------------------------------
fx <- list(datpath = "../data/",
           task = "CC",
           jmax = 2,
           dv = "dens_fx",
           sel_n = paste(c(25, 59, 136, 313)),
           w = 1.96,
           h = 2.36 * 2,
           xlabs = c(expression(eta[p]^2), expression("r"^2)),
           xl = c(0, 0.25),
           max_idx = c(20, 20),
           leg_id = 2,
           leg_locs = c(0.1, 300),
           figlabel = "B",
           figlabelon = TRUE)

p <- list(datpath = "../data/",
          task = "CC",
          jmax = 2,
          dv = "dens_p",
          sel_n = paste(c(25, 59, 136, 313)),
          w = 1.96,
          h = 2.36 * 2,
          xlabs = c("p", "p"),
          xl = c(-10, 10),
          max_idx = c(5, 20),
          leg_id = 1,
          leg_locs = c(0, 0.45))

# ----------------------------------------------------
# KL divergence
# ----------------------------------------------------

kl <- list(datpath = "../data/",
           task = "CC",
           dv = "dens_fx",
           ratio_type = "KL",
           origin = "313",
           sub_Ns = sub_Ns,
           w = 1.96,
           h = 2.36,
           leg_id = FALSE,
           leg_locs = c(5, 20),
           leg_txt = c("RM-AN", "LME"),
           ylabel = expression(italic("KL p||q")))

# ----------------------------------------------------
# fx sz ratio between models
# ----------------------------------------------------
           
model_rats <- list(datpath = "../data/",
              task = "CC",
              dv = "stats_fx",
              ratio_type = "model",
              origin = "",
              sub_Ns = sub_Ns,
              w = 1.96,
              h = 2.36,
              leg_id = TRUE,
              leg_locs = c(5, 20),
              leg_txt = "",
              ylabel = expression(italic("RM-AN / LME")),
              mods = c("RM-AN", "LME"))

# ----------------------------------------------------
# meta-analytic vs observed fx sz ratio
# ----------------------------------------------------
sig <- list(datpath = "../data/",
            task = "CC",
            dv = "stats_sig",
            ratio_type = "stats_sig",
            origin = "",
            sub_Ns = sub_Ns,
            w = 1.96,
            h = 2.36,
            leg_id = TRUE,
            leg_locs = c(5, 0.5),
            leg_txt = c("RM-AN", "LME"),
            ylabel = expression(italic("meta / sim")),
            mods = c("RM-AN", "LME"))

# ----------------------------------------------------
# plotting
# ----------------------------------------------------
pdf(paste("../images/", task, "_", "fx_main", ".pdf", sep = ""),
          width = w, height = h)
par(mfrow = c(2, 3), mar = c(3, 3, 0, 0),
    oma = c(2, 1, 2, 1),
    mgp = c(3, 1, 0), las = 0)
plot_CC_results(fname)
fig_label("A", cex = 2)
plot_dens(fx)
plot_ratios(kl)
fig_label("C", cex = 2)
plot_ratios(sig)
fig_label("D", cex = 2)
plot_ratios(model_rats)
fig_label("E", cex = 2)
dev.off()

#plot_dens(p)
