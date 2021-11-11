# # use this code to compile all the subplots for
# the behaviour and fx size into one 6 panel plot

rm(list = ls())
# -----------------------------------------------------------------
# load packages and source function files
# -----------------------------------------------------------------
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(lme4) # for mixed effects modelling
library(ggridges)
library(car)
library(parallel)
library(tidyverse) # for data wrangling
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
# VSL med = 
# ----------------------------------------------------
# CONSTANT SETTINGS
# ----------------------------------------------------
sub_Ns <- paste(round(exp(seq(log(13), log(313), length.out = 20))))
w = 1.96 * 3 # width of the plot, in inches
h = 2.36 * 2 # height

# ----------------------------------------------------
# TASK SETTINGS
# ----------------------------------------------------
task = "VSL"
load(paste("../data/", task, "/",
            task, "_plot_settings.RData", sep = ""))

# ----------------------------------------------------
# plotting
# ----------------------------------------------------
# use this for plots containing only one test
# i.e. AB, SD & SRT
pdf(paste("../images/", task, "_", "fx_main", ".pdf", sep = ""),
    width = w, height = h)
plot.mat = matrix(c(1, 1, 1, 2, 2, 2,
                    3, 3, 4, 4, 5, 5),
                  nrow = 2, byrow = T)
layout(plot.mat)
layout.show(n = 5)

plot_AB_results(fname)
fig_label("A", cex = 2)



pdf(paste("../images/", task, "_", "fx_main", ".pdf", sep = ""),
          width = w, height = h)
par(mfrow = c(2, 3), mar = c(3, 3, 1, 1),
    oma = c(1, 2, 1, 1),
    mgp = c(2, 1, 0), las = 0)
#plot_AB_results(fname)
#plot_CC_results(fname)
#plot_MT_results(fname)
#plot_SRT_results(fname)
plot_VSL_results(fname)
fig_label("A", cex = 2)
plot_dens(fx)
plot_qq_med_vs_best(qq_inputs)
fig_label("C", cex = 2)
plot_ratios(kl)
fig_label("D", cex = 2)
#plot_mean_diff_between_mods(model_mu_diff)
plot_mean_vs_meta(meta_mu)
fig_label("E", cex = 2)
dev.off()

pdf(paste("../images/", task, "_", "ps", ".pdf", sep = ""),
          width = w*.5, height = h)
par(mfrow = c(3, 1), mar = c(4, 3, 0, 0),
oma = c(1, 1, 1, 1),
mgp = c(2, 1, 0), las = 0)
plot_dens(p)
fig_label("B", cex = 2)
plot_ratios(p_rat)
abline(h=1, lty = 2, col = "grey48")
fig_label("C", cex = 2)
dev.off()
