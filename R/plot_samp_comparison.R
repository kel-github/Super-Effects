### written by K. Garner, Jan 2021
### Call plotting functions for observed effect sizes and p values
rm(list = ls())

# ----------------------------------------------------
# load packages and source function files
# ----------------------------------------------------

library(tidyverse) # for data wrangling
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(lme4) # for mixed effects modelling
library(ggridges)
library(car)
library(parallel)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting

# ----------------------------------------------------
# LIST OF SETTINGS
# ----------------------------------------------------
# AB med = 24
# CC med = 23
# SRT med = 39
# SD med = 24
sub_Ns <- paste(round(exp(seq(log(13), log(313), length.out = 20))))
# ----------------------------------------------------
# density variables
# ----------------------------------------------------
fx <- list(datpath = "../data/",
           task = "AB",
           jmax = 2,
           dv = "dens_fx",
           sel_n = paste(c(25, 59, 136, 313)),
           w = 1.96,
           h = 2.36 * 2,
           xlabs = c(expression(eta[p]^2), expression("r"^2)),
           xl = c(0, 0.25),
           max_idx = c(20, 20),
           leg_id = 2,
           leg_locs = c(0.1, 300))


# ----------------------------------------------------
# define functions
# ----------------------------------------------------
plot_dens <- function(inputs4plot) {
    # plot the density functions for select n
    # saves to a png file in the project image folder
    # kwargs
    # inputs4plot is a list containing the following
    # fields
    # -- datpath: e.g. "../data/"
    # -- task: "CC" or "AB" etc
    # -- jmax: max number of plots for loop
    # -- dv: "dens_fx" or "dens_p"
    # -- sel_n: plot densities from which subjects? c(n1, n2, ..., nn)
    # -- w: width in inches
    # -- h: height in inches
    # -- xlabs: vector of labels for x axis, 1 per plot
    # -- xl: xlims
    # -- max_idx: from which element should the max y value be taken from
    #             number corresponds to idx for sub.Ns
    # -- leg_id: on which plot would you like the legend? (1, or, 2)
    # -- leg_locs: x and y coordinates for the legend
    datpath <- inputs4plot$datpath
    task <- inputs4plot$task
    jmax <- inputs4plot$jmax
    dv <- inputs4plot$dv
    sel_n <- inputs4plot$sel_n
    w <- inputs4plot$w
    h <- inputs4plot$h
    xlabs <- inputs4plot$xlabs
    xl <- inputs4plot$xl
    max_idx <- inputs4plot$max_idx
    leg_id <- inputs4plot$leg_id
    leg_locs <- inputs4plot$leg_locs

    # ----------------------------------------------------
    # load data
    # ----------------------------------------------------
    fs <- list.files(paste(datpath, task, sep=""), "*.RData")
    get_datas <- function(datpath, task){    
        load(paste(datpath, task, "/", fs[1], sep=""))
        d <- 
    }
    load(paste(datpath, task, "/", task, "stats.RData", sep = ""))

    pdf(paste("../images/", task, "_", dv, ".pdf", sep = ""),
        width = w, height = h)
    par(mfrow = c(jmax, 1), mar = c(3, 3, 3, 1), mgp = c(2, 1, 0), las = 1)
    for (j in 1:jmax) {
        plot(res[sel_n[1], ][[dv]][[j]],
             col = wes_palette("IsleofDogs1")[1],
             lwd = 3,
             ylim = c(0,
                      max(res[max_idx[j], ][[dv]][[j]]["y"]$y)),
             xlim = xl,
             main = " ", ylab = "density",
             xlab = xlabs[j],
             bty = "n",
             cex.lab = 0.75,
             cex.axis = 0.5)
        polygon(res[sel_n[1], ][[dv]][[j]],
                col = adjustcolor(wes_palette("IsleofDogs1")[1],
                alpha.f = 0.5))
    for (i in c(2:length(sel_n))) {
        lines(res[sel_n[i], ][[dv]][[j]],
              col = wes_palette("IsleofDogs1")[i], lwd = 2)
        polygon(res[sel_n[i], ][[dv]][[j]],
              col = adjustcolor(wes_palette("IsleofDogs1")[i],
              alpha.f = 0.5))
    }
    # if a p statistic plot, add the criteria for significance
    if (dv == "dens_p") {
        abline(v = qnorm(.05), col = "#161616",
               lty = 2, lwd = 1)
    }
    # do you want a legend, and if so where to put it?
    if (j == leg_id) {
        leg_cols <- adjustcolor(wes_palette("IsleofDogs1")[c(1:4)],
                                alpha.f = 0.5)
        legend(leg_locs[1], leg_locs[2], legend = sel_n,
               col = leg_cols, lty = 1, bty = "n", cex = 0.75)
    }
 }
 dev.off()
}

# ----------------------------------------------------
# plotting
# ----------------------------------------------------
plot_dens(fx)
plot_dens(p)
plot_ratios(kl)
plot_ratios(model_rats)
plot_ratios(sig)
