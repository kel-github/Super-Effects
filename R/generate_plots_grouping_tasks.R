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
AB_med = 24
CC_med = 23
SRT_med = 39
SD_med = 24
med_line_col <- "grey81"
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

# ----------------------------------------------------
# EXECUTIVE FUNCTION
# ----------------------------------------------------
w = 12/2.54
h = 17/2.54
pdf(paste("../images/", "EPS", "_", "EF_tasks", "_", "fx_sz", ".pdf", sep = ""),
    width = w, height = h)  
plot.mat = matrix(c(1, 1, 2, 2, 
                    3, 3, 4, 4, 
                    5, 5, 6, 6),
                  nrow = 3, byrow = T)
layout(plot.mat)
par(las=1)

# Ab effect sizes and over N
load(sprintf(st_tmplt, "AB", "AB"))
plot_dens(fx)
fig_label("A", cex = 2)
plot_mus_n_quants(fxovrN)
# add legend
pal <- wesanderson::wes_palette("IsleofDogs2")[c(4,3)]
names(pal) <- c("mu", "q")
legend(189, 0.4, legend=c(expression(mu), "q"), col=pal, lty=1, bty = "n")
# add the median N line
abline(v=AB_med, lty=2, lwd=2, col=med_line_col)
fig_label("B", cex = 2)
# MT effect for the main effect of task
load(sprintf(st_tmplt, "SD", "SD"))
plot_dens(fx)
fig_label("C", cex = 2)
plot_mus_n_quants(fxovrN)
abline(v=SD_med, lty=2, lwd=2, col=med_line_col)
fig_label("D", cex = 2)
# MT effect for task x trial type interaction
plot_dens(fx, mods="LME")
fig_label("E", cex = 2)
fxovrN$mod <- "LME"
plot_mus_n_quants(fxovrN)
abline(h=0, lty=2, col="grey")
abline(v=SD_med, lty=2, lwd=2, col=med_line_col)
fig_label("F", cex = 2)

dev.off()

# ----------------------------------------------------
# IMPLICIT LEARNING
# ----------------------------------------------------

pdf(paste("../images/", "EPS", "_", "IL_tasks", "_", "fx_sz", ".pdf", sep = ""),
    width = w, height = h)  
plot.mat = matrix(c(1, 1, 2, 2, 
                    3, 3, 4, 4, 
                    5, 5, 6, 6),
                  nrow = 3, byrow = T)
layout(plot.mat)
par(las=1)

# first SRT
load(sprintf(st_tmplt, "SRT", "SRT"))
plot_dens(fx)
fig_label("A", cex = 2)
plot_mus_n_quants(fxovrN)
legend(189, 5, legend=c(expression(mu), "q"), col=pal, lty=1, bty = "n")
abline(v=SRT_med, lty=2, lwd=2, col=med_line_col)
fig_label("B", cex = 2)

# now contextual cueing interaction
load(sprintf(st_tmplt, "CC", "CC"))
plot_dens(fx)
fig_label("C", cex = 2)
plot_mus_n_quants(fxovrN)
abline(h=0, lty=2, col="grey")
abline(v=CC_med, lty=2, lwd=2, col=med_line_col)
fig_label("D", cex = 2)

# and contextual cueing ME
plot_dens(fx, mods="LME")
fig_label("E", cex = 2)
fxovrN$mod <- "LME"
plot_mus_n_quants(fxovrN)
abline(h=0, lty=2, col="grey")
abline(v=CC_med, lty=2, lwd=2, col=med_line_col)
fig_label("F", cex = 2)

dev.off()

# ----------------------------------------------------
# P-VALUES EXECUTIVE FUNCTION TASKS
# ----------------------------------------------------
pdf(paste("../images/", "EPS", "_", "EF_tasks", "_", "ps", ".pdf", sep = ""),
    width = w, height = h)  
plot.mat = matrix(c(1, 1, 2, 2, 
                    3, 3, 4, 4, 
                    5, 5, 6, 6),
                  nrow = 3, byrow = T)
layout(plot.mat)
par(las=1)

# Ab effect sizes and over N
load(sprintf(st_tmplt, "AB", "AB"))
plot_dens(p)
fig_label("A", cex = 2)
plot_mus_n_quants(povrN, dv="stats_p")
# add legend
pal <- wesanderson::wes_palette("IsleofDogs2")[c(4,3)]
names(pal) <- c("mu", "q")
legend(224, 0, legend=c(expression(mu), "q"), col=pal, lty=1, bty = "n")
# add the median N line
abline(v=AB_med, lty=2, lwd=2, col=med_line_col)
abline(h=qnorm(.05), lty=2, lwd=1, col='grey')
fig_label("B", cex = 2)

# now the main effect for multitasking
load(sprintf(st_tmplt, "SD", "SD"))
plot_dens(p)
fig_label("C", cex=2)
plot_mus_n_quants(povrN, dv="stats_p")
abline(v=SD_med, lty=2, lwd=2, col=med_line_col)
abline(h=qnorm(.05), lty=2, lwd=1, col='grey')
fig_label("D", cex=2)

p$leg_id <- 1
p$leg_locs <- c(-60, 0.3)
plot_dens(p, mod="LME")
fig_label("E", cex=2)
povrN$mod <- "LME"
plot_mus_n_quants(povrN, dv="stats_p")
abline(v=SD_med, lty=2, lwd=2, col=med_line_col)
abline(h=qnorm(.05), lty=2, lwd=1, col='grey')
fig_label("F", cex=2)
dev.off()

# ----------------------------------------------------
# P-VALUES IMPLICIT LEARNING TASKS
# ----------------------------------------------------
pdf(paste("../images/", "EPS", "_", "IL_tasks", "_", "ps", ".pdf", sep = ""),
    width = w, height = h)  
plot.mat = matrix(c(1, 1, 2, 2, 
                    3, 3, 4, 4, 
                    5, 5, 6, 6),
                  nrow = 3, byrow = T)
layout(plot.mat)
par(las=1)

# SRT effect sizes and over N
load(sprintf(st_tmplt, "SRT", "SRT"))
plot_dens(p)
fig_label("A", cex = 2)
plot_mus_n_quants(povrN, dv="stats_p")
#legend(13, -20, legend=c(expression(mu), "q"), col=pal, lty=1, bty = "n")
# add the median N line
abline(v=SRT_med, lty=2, lwd=2, col=med_line_col)
abline(h=qnorm(.05), lty=2, lwd=1, col='grey')
fig_label("B", cex = 2)

# now CC interaction
load(sprintf(st_tmplt, "CC", "CC"))
plot_dens(p)
fig_label("C", cex=2)
plot_mus_n_quants(povrN, dv="stats_p")
abline(v=CC_med, lty=2, lwd=2, col=med_line_col)
abline(h=qnorm(.05), lty=2, lwd=1, col='grey')
fig_label("D", cex=2)

p$leg_id <- FALSE
plot_dens(p, mod="LME")
fig_label("E", cex=2)
povrN$mod <- "LME"
plot_mus_n_quants(povrN, dv="stats_p")
abline(v=CC_med, lty=2, lwd=2, col=med_line_col)
abline(h=qnorm(.05), lty=2, lwd=1, col='grey')
fig_label("F", cex=2)
legend(59, -10, legend=c(expression(mu), "q"), col=pal, lty=1, bty = "n")
dev.off()

# ----------------------------------------------------
# FILE DRAWER TASK CONTEXTUAL CUEING
# ----------------------------------------------------
w <- 8/2.54
h <- 8/2.54
pdf(paste("../images/", "EPS", "_", "IL_tasks", "_", "ps", ".pdf", sep = ""),
    width = w, height = h)  
load(sprintf(st_tmplt, "CC", "CC"))
meta_mu$eps <- TRUE
plot_mean_vs_meta(meta_mu)
