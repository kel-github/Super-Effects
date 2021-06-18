### written by K. Garner, Jan 2021
### for the project 'On the detectability of effects
# in executive function and implicit learning tasks
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE

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
# define session variables
# ----------------------------------------------------
datpath <- "../data/"
task <- "CC"
med <- 23

# ----------------------------------------------------
# LIST OF SETTINGS
# ----------------------------------------------------
# AB med = 24
# CC med = 23
# SRT med = 39
# SD med = 24

# ----------------------------------------------------
# load data
# ----------------------------------------------------
load(paste(datpath, task, "/", task, "stats.RData", sep = ""))

# ----------------------------------------------------
# session variables
# ----------------------------------------------------
sel_n <- paste(c(25, 59, 136, 313))
xl <- c(-10, 10)
jmax <- 2
w <- 1.96
h <- 2.36
dv <- "dens_fx"
xlabs <- c("a", "b") # to be fixed
max_idx <- c(5, 20)

# ----------------------------------------------------
# define functions
# ----------------------------------------------------
plot_dens <- function(task, jmax, dv, w, h, xlabs, xl, max_idx) {
    # plot the density functions for select n
    # saves to a png file in the project image folder
    # kwargs
    # -- task: "CC" or "AB" etc
    # -- jmax: max number of plots for loop
    # -- dv: "dens_fx" or "dens_p"
    # -- w: width in inches
    # -- h: height in inches
    # -- xlabs: vector of labels for x axis, 1 per plot
    # -- xl: xlims
    # -- max_idx: from which element should the max y value be taken from
    #             number corresponds to idx for sub.Ns

    pdf(paste("../images/", task, "_", dv, ".pdf", sep = ""),
        width = w, height = h)
    par(mfrow = c(jmax, 1), mar = c(2, 3, 3, 1), mgp = c(2, 1, 0), las = 1)
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
 }
 dev.off()
}

# ----------------------------------------------------
# plotting
# ----------------------------------------------------

plot_dens(task, jmax, dv = "dens_p", w, h, xlabs, xl, max_idx)
