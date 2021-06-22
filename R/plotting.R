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

KL <- list(datpath = "../data/",
           task = "CC",
           dv = "dens_fx",
           ratio_type = "KL",
           origin = "313",
           sub_Ns = sub_Ns,
           w = 1.96,
           h = 2.36,
           leg_id = TRUE,
           leg_locs = c(5, 20),
           leg_txt = c("RM-AN", "LME"),
           ylabel = expression(italic("KL p||q")))
           
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
sig <- list(datpath = "../data/",
            task = "CC",
            dv = "sig",
            ratio_type = "stats_sig",
            origin = "",
            sub_Ns = sub_Ns,
            w = 1.96,
            h = 2.36,
            leg_id = TRUE,
            leg_locs = c(5, 20),
            leg_txt = "",
            ylabel = expression(italic("meta / sim")),
            mods = c("RM-AN", "LME"))

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
    mods <- unique(unique(colnames(res[, "stats_fx"][[origin]])))
    # first calculate the approximating function for the all the densities
    d.funcs <- lapply(sub_Ns, function(j) 
                              lapply(mods, function(z) 
                                     approxfun(x=res[ , dv][[j]][[z]]$x, 
                                     y=res[ , dv][[j]][[z]]$y/length(res[ , dv][[j]][[z]]$y), 
                                     rule=1)))
    names(d.funcs) <- sub_Ns
    
    # then select a series of bin widths, or xs
    xvals <- seq(0,1, by=.001) # cos our effect sizes are positive correlations
    # do the same for the second distribution
    # calculate the KL divergence between the two
    KLfunc <- function(pf, qf, x){
        # -- pf: function for p - P IS THE P DIST
        # -- pq: function for q - Q IS THE APPROXIMATING DIST
        # -- x: values of x
        p <- pf(x)
        q <- qf(x)
        sum(p*log2(p/q), na.rm=T)
    }
    # next, match binwidths between the 

    KL4plotting <- lapply(sub_Ns, function(x) 
                                  lapply(1:length(mods), 
                                         function(y) 
                                         KLfunc(pf=d.funcs[[origin]][[y]],
                                                qf=d.funcs[[x]][[y]],
                                                x=xvals)))
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
    dv <- rat_inputs$dv
    sub_Ns <- rat_inputs$sub_Ns
    mods <- rat_inputs$mods
    ratios4plotting <- do.call(rbind, lapply(sub_Ns, 
                                      function(x) abs(diff(res[,dv][[x]][3, mods[1]][[1]])) /
                                                  abs(diff(res[,dv][[x]][3, mods[2]][[1]])) ))
    colnames(ratios4plotting) <- "ratio"
    rownames(ratios4plotting) <- sub_Ns
    ratios4plotting   
}

calc_meta_vs_model <- function(rat_inputs, res){
    # calculate ratio between the mu's for the sig and
    # the sim distributions
    sub_Ns <- rat_inputs$sub_Ns
    mods <- rat_inputs$mods
    dv <- rat_inputs$dv
    ratios4plotting <- do.call(rbind, lapply(sub_Ns, function(y)
                                lapply( mods, function(x)
                                res[ , dv][[y]][ , x][[1]] /
                                res[ ,"stats_fx"][[y]][ , x][[1]])))
    colnames(ratios4plotting) <- mods
    rownames(ratios4plotting) <- sub_Ns
    ratios4plotting
}

plot_ratios <- function(rat_inputs) {
    # plot ratio between x and y
    # kwargs
    # -- rat_inputs: a list comprising of:
    # -- datpath: e.g. "../data/
    # -- task: e.g. "CC"
    # -- dv: e.g. "dens_fx", "stats_fx", "stats_p", "stats_sig"
    # -- ratio_type: "origin" or "model" or "KL"
    # -- origin <- e.g. "313" i.e. calc ratio to this one
    # -- sub_Ns <- names of subject groups (Ns)
    # -- w width of plot in inches
    # -- h: height of plot in inches
    # -- leg_id: legend? TRUE or FALSE
    # -- leg_locs: e.g. c(5, 20)
    # -- leg_txt: e.g. c("RM-AN", "LME")
    datpath <- rat_inputs$datpath
    task <- rat_inputs$task
    dv <- rat_inputs$dv
    ratio_type <- rat_inputs$ratio_type
    ylabel <- rat_inputs$y_label
    sub_Ns <- rat_inputs$sub_Ns
    leg_id <- rat_inputs$leg_id
    leg_locs <- rat_inputs$leg_locs
    leg_txt <- rat_inputs$leg_txt
    w <- rat_inputs$w
    h <- rat_inputs$h
    ylabel <- rat_inputs$ylabel
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
    } else if (ratio_type == "KL"){
        ratios <- calc_KL_sing_origin(rat_inputs, res) 
    } else if (ratio_type == "stats_sig") {
        ratios <- calc_meta_vs_model(rat_inputs, res)
    }

    # ----------------------------------------------------
    # plot ratios
    # ----------------------------------------------------
    pdf(paste("../images/", task, "_", ratio_type, ".pdf", sep = ""),
        width = w, height = h)
    par(mfrow = c(1, 1), mar = c(3, 3, 3, 1), mgp = c(2, 1, 0), las = 1)
    if (ncols(ratios) > 1) {
        plot(x=1:length(rownames(ratios)), y = ratios[,"RM-AN"], 
             xaxt = "n",
             bty = "n",
             type = "l", lty = 1,
             lwd = 2,
             col = wes_palette("IsleofDogs1")[6],
             ylim = c(0, 
                      max(cbind(do.call(cbind, ratios[ , "RM-AN"]), 
                                do.call(cbind, ratios[ , "LME"])))),
             ylab = ylabel,
             xlab = expression(italic("N")),
             cex.lab = 0.75,
             cex.axis = 0.5)
        axis(1, at = seq(1, 20, 2), 
             labels = sub_Ns[seq(1, 20, 2)], 
             cex.lab = 0.75,
             cex.axis = 0.5)
        lines(x=1:length(rownames(ratios)), y = ratios[, "LME"],
              lwd = 2,
              col = wes_palette("IsleofDogs1")[5])
        # do you want a legend, and if so where to put it?
        if (leg_id) {
            leg_cols <- wes_palette("IsleofDogs1")[c(6, 5)]
            legend(leg_locs[1], leg_locs[2], legend = leg_txt,
                   col = leg_cols, lty = 1, bty = "n", cex = 0.75)
        }
        if (ratio_type == "stats_sig"){
            abline(h=1, lty=2, col="grey48")
        }
    } else if (ncols(ratios) == 1){
        plot(x=1:length(rownames(ratios)), y = ratios[,"ratio"], 
             xaxt = "n",
             bty = "n",
             type = "l", lty = 1,
             lwd = 2,
             col = wes_palette("IsleofDogs1")[4],
             ylim = c(0, 
                      max(ratios[,"ratio"])+0.5),
             ylab = ylabel,
             xlab = expression(italic("N")),
             cex.lab = 0.75,
             cex.axis = 0.5)
        axis(1, at = seq(1, 20, 2), 
             labels = sub_Ns[seq(1, 20, 2)], 
             cex.lab = 0.75,
             cex.axis = 0.5)
    } else if (ratio_type == "stats_sig"){
        
        
    }
    dev.off()
}

# ----------------------------------------------------
# plotting
# ----------------------------------------------------
plot_dens(fx)
plot_dens(p)
plot_ratios(KL)
plot_ratios(model_rats)
