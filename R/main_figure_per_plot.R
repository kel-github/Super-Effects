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

# ----------------------------------------------------
# TASK SETTINGS
# ----------------------------------------------------
task = "AB"
load(paste("../data/", task, "/",
            task, "_plot_settings.RData", sep = ""))

# ----------------------------------------------------
# plotting
# ----------------------------------------------------
pdf(paste("../images/", task, "_", "fx_main", ".pdf", sep = ""),
          width = w, height = h)
par(mfrow = c(2, 3), mar = c(3, 3, 1, 1),
    oma = c(1, 2, 1, 1),
    mgp = c(2, 1, 0), las = 0)
plot_AB_results(fname)
fig_label("A", cex = 2)
plot_dens(fx)
plot_ratios(kl)
fig_label("C", cex = 2)
plot_ratios(sig)
fig_label("D", cex = 2)
plot_ratios(model_rats)
fig_label("E", cex = 2)
dev.off()

pdf(paste("../images/", task, "_", "ps", ".pdf", sep = ""),
          width = w*.5, height = h)
par(mfrow = c(2, 1), mar = c(4, 3, 0, 0),
oma = c(1, 1, 1, 1),
mgp = c(2, 1, 0), las = 0)
plot_dens(p)
fig_label("B", cex = 2)
dev.off()