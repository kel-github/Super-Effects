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
           leg_locs = c(0.1, 300))

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
    # -- sel_n: overlay densities from which subjects?
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

calc_ratios_sing_origin <- function(rat_inputs, res){
    # calculate ratio between x_1 and y_1_to_j
    # NOTE: this is for the ratio between 95% of 
    # distribution ONLY!
    # kwargs
    # -- rat_inputs: a list comprising of the following fields
    # origin <- e.g. "313"
    # dv <- "stats_fx", "stats_p", "stats_sig"
    # sub_Ns <- names of subject groups (Ns)
    
    # -- res: stored results from get_statistics.R
    origin <- rat_inputs$origin
    sub_Ns <- rat_inputs$sub_Ns
    dv <- rat_inputs$dv
    # get models from which to return results
    mods <- unique(colnames(res[, dv][[origin]]))
    base <- do.call(cbind, lapply(mods, 
                                  function (x) abs(diff(res[, dv][[origin]][, x][[3]]))))
    colnames(base) <- mods # just labeling for sanity checks
    ratios4plotting <- do.call(rbind, lapply(sub_Ns, function(x) abs(diff(res[, dv][[x]][[3]])) / base))
    rownames(ratios4plotting) <- sub_Ns
    ratios4plotting
}

calc_KL_sing_origin <- function(rat_inputs, res){
    #### this project on hold
    # calculate ratio between x_1 and y_1_to_j
    # NOTE: this is for the ratio between 95% of 
    # distribution ONLY!
    # kwargs
    # -- rat_inputs: a list comprising of the following fields
    # origin <- e.g. "313"
    # dv <- "stats_fx", "stats_p", "stats_sig"
    # sub_Ns <- names of subject groups (Ns)
    
    # -- res: stored results from get_statistics.R
    origin <- rat_inputs$origin
    sub_Ns <- rat_inputs$sub_Ns
    dv <- rat_inputs$dv
    # get models from which to return results
    mods <- unique(colnames(res[, dv][[origin]]))
    base <- do.call(cbind, lapply(mods, 
                                  function (x) list(mu = res[ , dv][[origin]][ , x][1][[1]],
                                                    sd = res[ , dv][[origin]][ , x][2][[1]]) ))
    colnames(base) <- mods # just labeling for sanity checks
    compare <- lapply(sub_Ns, 
                      function (x) lapply(mods,
                                          function(y) list(mu = res[ , dv][[x]][ , y][1][[1]],
                                                           sd = res[ , dv][[x]][ , y][2][[1]])))
    compare <- do.call(rbind, compare)
    colnames(compare) <- mods
    rownames(compare) <- sub_Ns
    getKL <- function( mu_a, sd_a, mu_b, sd_b ) {
        # https://stats.stackexchange.com/questions/7440/kl-divergence-between-two-univariate-gaussians
        log((sd_b/sd_a) + ((sd_a + (mu_a - mu_b)^2) / (2*sd_b)) - 0.5)
    }
    
    KL4plotting <- lapply(sub_Ns, function(x) 
                                  lapply(mods, function(y) 
                                               getKL(compare[x , y][[1]]$mu,
                                                     compare[x , y][[1]]$sd,
                                                     base[, y]$mu,
                                                     base[, y]$sd)))
    KL4plotting <- do.call(rbind, KL4plotting)
    colnames(KL4plotting) <- mods
    rownames(KL4plotting) <- sub_Ns
    KL4plotting
}

calc_ratios_by_model <- function(rat_inputs, res){
    # calculate ratio between x_1_to_j and y_1_to_j
    # NOTE: this is for the ratio between 95% of 
    # distribution ONLY!
    # Gives mod_a - mod_b
    # kwargs
    # -- rat_inputs: a list comprising of the following fields
    # -- dv: "stats_fx", "stats_p", "stats_sig"
    # -- sub_Ns: names of subject groups (Ns)
    # -- mods: names of the two models to compare e.g c("RM-AN", "LME")
    ratios4plotting <- do.call(rbind, lapply(sub_Ns, 
                                      function(x) abs(diff(res[,dv][[x]][3, mods[1]][[1]])) /
                                                  abs(diff(res[,dv][[x]][3, mods[2]][[1]])) ))
    colnames(ratios4plotting) <- "ratio"
    rownames(ratios4plotting) <- sub_Ns
    ratios4plotting   
}

plot_ratios <- function(rat_inputs) {
    # plot ratio between x and y
    # kwargs
    # -- rat_inputs: a list comprising of:
    # -- datpath: e.g. "../data/
    # -- task: e.g. "CC"
    # -- dv: "stats_fx", "stats_p", "stats_sig"
    # -- ratio_type: "origin" or "model"
    # -- origin <- e.g. "313" i.e. calc ratio to this one
    # -- dv <- "stats_fx", "stats_p", "stats_sig"
    # -- sub_Ns <- names of subject groups (Ns)
    datpath <- rat_inputs$datpath
    task <- rat_inputs$task
    dv <- rat_inputs$dv
    ratio_type <- rat_inputs$ratio_type
    
    # ----------------------------------------------------
    # load data
    # ----------------------------------------------------
    load(paste(datpath, task, "/", task, "stats.RData", sep = ""))
    
    # ----------------------------------------------------
    # compute ratios
    # ----------------------------------------------------
    if (ratio_type == "origin") {
        ratios <- calc_ratios_sing_origin(rat_inputs, res)
    } else if (ratio_type == "model") {
        ratios <- calc_ratios_by_model(rat_inputs, res)
    }
    
    # ----------------------------------------------------
    # plot ratios
    # ----------------------------------------------------
}

# ----------------------------------------------------
# plotting
# ----------------------------------------------------
plot_dens(fx)
plot_dens(p)
