# # use this code to compile all the subplots for
# the behaviour and fx size into one 6 panel plot

rm(list = ls())
# -----------------------------------------------------------------
# load packages and source function files
# -----------------------------------------------------------------
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(lme4) # for mixed effects modelling
#library(ggridges)
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
tasks = c("AB", "SD", "SRT", "CC")


# ----------------------------------------------------
# BEHAVIOURAL RESULTS PLOT
# ----------------------------------------------------
fnm_tmplt <- "../data/total_of_313_subs_%s_task_trial_level_data.csv"
pdf(paste("../images/", "EPS", "_", "all_tasks", "_", "behav", ".pdf", sep = ""),
    width = w, height = h)  
plot.mat = matrix(c(1, 1, 1, 2, 2, 2,
                    3, 3, 3, 4, 4, 4),
                  nrow = 2, byrow = T)
layout(plot.mat)
par(las=1)
plot_AB_results(sprintf(fnm_tmplt, "AB"))
fig_label("A", cex = 2) 
plot_MT_results(sprintf(fnm_tmplt, "SingDual"))
fig_label("B", cex = 2) 
plot_SRT_results(sprintf(fnm_tmplt, "SRT"))
fig_label("C", cex = 2) 
plot_CC_results(sprintf(fnm_tmplt, "CC"), type="mean")
fig_label("D", cex = 2) 
dev.off()

# ----------------------------------------------------
# EFFECT SIZES PLOT
# ----------------------------------------------------
st_tmplt <- "../data/%s/EPS%s_plot_settings.RData"

pdf(paste("../images/", "EPS", "_", "all_tasks", "_", "fx_sz", ".pdf", sep = ""),
    width = w, height = h)  
plot.mat = matrix(c(1, 1, 2, 2, 3, 3,
                    4, 4, 5, 5, 6, 6),
                  nrow = 3, byrow = T)
layout(plot.mat)
par(las=1)

load(sprintf(st_tmplt, "AB", "AB"))
plot_dens(fx)
fig_label("A", cex = 2)
load(sprintf(st_tmplt, "SD", "SD"))
plot_dens(fx)
fig_label("B", cex = 2)
load(sprintf(st_tmplt, "SRT", "SRT"))
plot_dens(fx)
fig_label("C", cex = 2)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)

palette_choice <- wes_palette("IsleofDogs1")[c(1, 2, 3, 5)]
leg_cols <-adjustcolor(palette_choice, alpha = 0.75)
legend(0.65, 1.2, legend = fx$sel_n,
       col = leg_cols, lty = 1, lwd = 3, bty = "n", cex = 1.5)



dev.off()
