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
           w = 1.96 * 2,
           h = 2.36 * 2,
           xlabs = expression(eta[p]^2),
           xl = c(0, 1),
           max_idx = c(20, 20),
           leg_id = 2,
           leg_locs = c(0.01, 10))


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
    fs <- list.files(paste(datpath, task, sep=""), "*stats.RData")

    get_a_d <- function(datpath, f, n){
        load(paste(datpath, task, "/", f, sep=""))
        res[,"dens_fx"][[n]][["RM-AN"]]
    }

    ds <- lapply(fs, function(x)
                 lapply(sel_n, function(y)
                 get_a_d(datpath, x, paste(y))))
    ds <- do.call(rbind, ds)
    colnames(ds) <- paste(sel_n)
    rownames(ds) <- fs

    pdf(paste("../images/", task, "_", "sampling", ".pdf", sep = ""),
        width = w, height = h)
    par(mfrow = c(2, 2), mar = c(3, 3, 3, 1), mgp = c(2, 1, 0), las = 1)

    # plot by subject
    make_ds <- function(ds, fs, ns) {
        plot(ds[fs[1], ns][[1]],
             col = wes_palette("IsleofDogs1")[1],
             lwd = 3,
             ylim = c(0,
                     max(c(ds[,ns][[1]]$y, ds[,ns][[2]]$y, ds[,ns][[3]]$y))),
             xlim = xl,
             ylab = "density",
             xlab = xlabs,
             bty = "n",
             cex.lab = 0.75,
             cex.axis = 0.5,
             main = paste(ns))
        polygon(ds[fs[1], ns][[1]],
                 col = adjustcolor(wes_palette("IsleofDogs1")[1],
                 alpha.f = 0.5))
        for (i in 2:length(fs)){
            lines(ds[fs[i], ns][[1]],
                 col = wes_palette("IsleofDogs1")[i], lwd = 2)
            polygon(ds[fs[i], ns][[1]],
                 col = adjustcolor(wes_palette("IsleofDogs1")[i],
                 alpha.f = 0.5))
        }
        if (ns == "25") {
            leg_cols <- adjustcolor(wes_palette("IsleofDogs1")[c(1:3)],
                                alpha.f = 0.5)
            legend(leg_locs[1], leg_locs[2], legend = c("j=1,k=1000", "j=k=1000", "j=0,k=1000"),
               col = leg_cols, lty = 1, bty = "n", cex = 0.75)
        }
    }

    lapply(paste(sel_n), make_ds, ds = ds, fs = fs)
    # do you want a legend, and if so where to put it?

    dev.off()
}

# ----------------------------------------------------
# plotting
# ----------------------------------------------------

plot_dens(fx)