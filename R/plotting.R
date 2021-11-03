### written by K. Garner, Jan 2021
### Call plotting functions for observed effect sizes and p values
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
    figlabel <- inputs4plot$figlabel
    figlabelon <- inputs4plot$figlabelon

    # ----------------------------------------------------
    # load data
    # ----------------------------------------------------
    load(paste(datpath, task, "/", task, "stats.RData", sep = ""))

    for (j in 1:jmax) {
        plot(res[sel_n[length(sel_n)], ][[dv]][[j]],
             col = wes_palette("IsleofDogs1")[1],
             lwd = 3,
             ylim = c(0,
                      max(res[max_idx[j], ][[dv]][[j]]["y"]$y)),
             xlim = xl,
             main = " ", ylab = "density",
             xlab = xlabs[j],
             bty = "n",
             cex.lab = 1,
             cex.axis = 1)
        polygon(res[sel_n[length(sel_n)], ][[dv]][[j]],
                col = adjustcolor(wes_palette("IsleofDogs1")[length(sel_n)],
                alpha.f = 0.5))
    for (i in c((length(sel_n)-1):1)) {
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
               col = leg_cols, lty = 1, bty = "n", cex = 1)
    }
    # add a label if desired
    if (j == 1 & figlabelon) {
        fig_label(figlabel, cex = 2)
    }
 }
}

calc_ratios_sing_origin <- function(rat_inputs, res) {
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
                           function(x)
                           abs(diff(res[, dv][[origin]][, x][[3]]))))
    colnames(base) <- mods # just labeling for sanity checks
    ratios4plotting <- do.call(rbind,
                                lapply(sub_Ns, 
                                function(x) lapply(mods, function(y)
                                abs(diff(res[, dv][[x]][,y][[3]]) / base[,y]))))
    rownames(ratios4plotting) <- sub_Ns
    ratios4plotting
}

calc_kl_sing_origin <- function(rat_inputs, res) {
    # compute kl divergence between density at given N relative to origin density
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
    d_funcs <- lapply(sub_Ns, function(j)
                              lapply(mods, function(z)
                                approxfun(x = res[, dv][[j]][[z]]$x,
                                y=res[ , dv][[j]][[z]]$y/length(res[ , dv][[j]][[z]]$y),
                                rule = 1)))
    names(d_funcs) <- sub_Ns
    # then select a series of bin widths, or xs
    xvals <- seq(0, 1, by = .001)
    # cos our effect sizes are positive correlations
    # do the same for the second distribution
    # calculate the KL divergence between the two
    kl_func <- function(pf, qf, x) {
        # -- pf: function for p - P IS THE P DIST
        # -- pq: function for q - Q IS THE APPROXIMATING DIST
        # -- x: values of x
        p <- pf(x)
        q <- qf(x)
        kl <- p * log2(p / q)
        sum(kl, na.rm = T)
    }
    kl4plotting <- lapply(sub_Ns, function(x)
                                  lapply(1:length(mods),
                                         function(y)
                                         kl_func(pf = d_funcs[[origin]][[y]],
                                                 qf = d_funcs[[x]][[y]],
                                                 x = xvals)))
    kl4plotting <- do.call(rbind, kl4plotting)
    colnames(kl4plotting) <- mods
    rownames(kl4plotting) <- sub_Ns
    kl4plotting
}

# calc_ratios_by_model <- function(rat_inputs, res){
#     # calculate ratio between x_1_to_j and y_1_to_j
#     # NOTE: this is for the ratio between 95% of
#     # distribution ONLY!
#     # Gives mod_a - mod_b
#     # kwargs
#     # -- rat_inputs: a list comprising of the following fields
#     # -- dv: "stats_fx", "stats_p", "stats_sig"
#     # -- sub_Ns: names of subject groups (Ns)
#     # -- mods: names of the two models to compare e.g c("RM-AN", "LME")
#     dv <- rat_inputs$dv
#     sub_Ns <- rat_inputs$sub_Ns
#     mods <- rat_inputs$mods
#     r4p <- do.call(rbind,
#                    lapply(sub_Ns,
#                      function(x) abs(diff(res[, dv][[x]]["qs", mods[1]][[1]])) /
#                                  abs(diff(res[, dv][[x]]["qs", mods[2]][[1]]))))
#     colnames(r4p) <- "ratio"
#     rownames(r4p) <- sub_Ns
#     r4p
# }

# calc_meta_vs_model <- function(rat_inputs, res){
#     # calculate ratio between the mu's for the sig and
#     # the sim distributions
#     sub_Ns <- rat_inputs$sub_Ns
#     mods <- rat_inputs$mods
#     dv <- rat_inputs$dv
#     ratios4plotting <- do.call(rbind, lapply(sub_Ns, function(y)
#                                 lapply( mods, function(x)
#                                 res[, dv][[y]][, x][[1]] /
#                                 res[, "stats_fx"][[y]][, x][[1]])))
#     colnames(ratios4plotting) <- mods
#     rownames(ratios4plotting) <- sub_Ns
#     ratios4plotting
# }

# plot_diagonal <- function(rat_inputs) {
#     # plot x vs y, draw a dotted line on
#     # the diagonal
#     # kwargs
#     # -- rat_inputs: a list comprising of
#     # -- datpath e.g. "../data/
#     # -- task: e.g. "CC"
#     # -- dvs: e.g. c("stats_fx", "stats_sig")
#     # -- sub_Ns <- names of subject groups (Ns)
#     # -- w width of plot in inches
#     # -- h: height of plot in inches
#     # -- axl: axis labels - e.g. c("sim", "M-A")
#     datpath <- rat_inputs$datpath
#     task <- rat_inputs$task
#     dvs <- rat_inputs$dvs
#     sub_Ns <- rat_inputs$sub_Ns
#     mods <- rat_inputs$mods
#     w <- rat_inputs$w
#     h <- rat_inputs$h
#     axl <- rat_inputs$axl
# 
#     # ----------------------------------------------------
#     # load data
#     # ----------------------------------------------------
#     load(paste(datpath, task, "/", task, "stats.RData", sep = ""))
#     # ----------------------------------------------------
#     # get xs and ys
#     # ----------------------------------------------------
#     xys <- do.call(rbind, lapply(sub_Ns, function(j)
#                           lapply(mods, function(z)
#                              data.frame(x = res[, dvs[1]][[j]][, z][[1]],
#                                         y = res[, dvs[2]][[j]][, z][[1]]))))
#     colnames(xys) <- mods
#     rownames(xys) <- sub_Ns
#     xyls <- list (do.call(rbind, xys[, mods[1]]),
#                   do.call(rbind, xys[, mods[2]]))
#     names(xyls) <- mods
#     # ----------------------------------------------------
#     # now plot the things
#     # ----------------------------------------------------
#     plot(x = xyls[[mods[1]]][, "x"],
#          y = xyls[[mods[1]]][, "y"],
#          bty = "n",
#          type = "l", lty = 1,
#          lwd = 2,
#          col = wes_palette("IsleofDogs1")[6],
#          ylim = c(0,
#                   max(do.call(rbind, xyls))),
#          xlim = c(0,
#                   max(do.call(rbind, xyls))),
#          xlab = axl[1],
#          ylab = axl[2],
#          cex.lab = 1,
#          cex.axis = 1)
#          lines(x = xyls[[mods[2]]][, "x"],
#                y = xyls[[mods[2]]][, "y"],
#                lwd = 2,
#                col = wes_palette("IsleofDogs1")[5])
#          abline(a = 0, b = 1, lty = 2, col = "grey48")
# }

plot_ratios <- function(rat_inputs) {
    # plot ratio between x and y
    # kwargs
    # -- rat_inputs: a list comprising of:
    # -- datpath: e.g. "../data/
    # -- task: e.g. "CC"
    # -- dv: e.g. "dens_fx", "stats_fx", "stats_p", "stats_sig"
    # -- ratio_type: "origin" or "model" or "KL" or "stats_sig"
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
    yl <- rat_inputs$yl
    if (ratio_type == "origin" | ratio_type == "KL") origin <- rat_inputs$origin
    
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
    } else if (ratio_type == "KL") {
        ratios <- calc_kl_sing_origin(rat_inputs, res)
    } else if (ratio_type == "stats_sig") {
        ratios <- calc_meta_vs_model(rat_inputs, res)
    }

    # ----------------------------------------------------
    # plot ratios
    # ----------------------------------------------------

    if (ncol(ratios) > 1) {
      if (is.null(yl)){
        nuyl = c(0,
                 max(cbind(do.call(cbind, ratios[, 1]),
                           do.call(cbind, ratios[, 2]))))
      } else {
        nuyl = yl
      }
        plot(x = 1:length(rownames(ratios)), y = ratios[, 1],
             xaxt = "n",
             bty = "n",
             type = "l", lty = 1,
             lwd = 2,
             col = wes_palette("IsleofDogs1")[6],
             ylim = nuyl,
             ylab = ylabel,
             xlab = expression(italic("N")),
             cex.lab = 1,
             cex.axis = 1)
        axis(1, at = seq(1, 20, 2),
             labels = sub_Ns[seq(1, 20, 2)],
             cex.lab = 1,
             cex.axis = 1)
        if (task == "SRT" & ratio_type == "KL"){
          points(x = seq_len(length(rownames(ratios))), y = ratios[, 2],
                 pch = 19, col = wes_palette("IsleofDogs1")[5])
        } else {
        lines(x = 1:length(rownames(ratios)), y = ratios[, 2],
              lwd = 2,
              col = wes_palette("IsleofDogs1")[5])
        }
        # do you want a legend, and if so where to put it?
        if (leg_id) {
            leg_cols <- wes_palette("IsleofDogs1")[c(6, 5)]
            legend(leg_locs[1], leg_locs[2], legend = leg_txt,
                   col = leg_cols, lty = 1, bty = "n", cex = 1)
        }
        if (ratio_type == "stats_sig"){
            abline(h=1, lty=2, col="grey48")
        }
    } else if (ncol(ratios) == 1) {
        if (is.null(yl)){
          nuyl = c(0, 
                   max(ratios[,"ratio"])+0.5)
        } else {
          nuyl = yl
        }
        plot(x = 1:length(rownames(ratios)), y = ratios[,"ratio"], 
             xaxt = "n",
             bty = "n",
             type = "l", lty = 1,
             lwd = 2,
             col = wes_palette("IsleofDogs1")[4],
             ylim = nuyl,
             ylab = ylabel,
             xlab = expression(italic("N")),
             cex.lab = 1,
             cex.axis = 1)
        axis(1, at = seq(1, 20, 2), 
             labels = sub_Ns[seq(1, 20, 2)], 
             cex.lab = 1,
             cex.axis = 1)
    }
}

# ----------------------------------------------------
# function to get pooled standard error
# ----------------------------------------------------
get_pooled <- function(vars_a, vars_b, n=1000^2){
  # return 2.5 times pooled standard error of the mean difference
  (sqrt(((n-1)*vars_a + (n-1)*vars_b) / (n + n - 2)) / sqrt(n))*1.96
}

plot_mean_vs_meta <- function(mu_vs_meta_inputs){
  # this function will plot the difference between
  # the observed and meta-analytic means
  # with means plotted as dots, and standard error
  # of the difference plotted as
  # error bars
  # puts a line at the bottom to denote sig
  # differences
  # plots a separate line for RM-AN & LME
  
  # -- mu_vs_meta_inputs: a list comprising of:
  # ----- datpath: e.g. "../data/
  # ----- task: e.g. "CC"
  # ----- sub_Ns <- names of subject groups paste(Ns)
  # ----- mods <- names of models to extract info for
  #               e.g. c("RM-AN", "LME")
  # ----- w width of plot in inches
  # ----- h: height of plot in inches
  # ----- leg_id: legend? TRUE or FALSE
  # ----- leg_locs: e.g. c(5, 20)
  # ----- leg_txt: e.g. c("RM-AN", "LME")
  # ----- yl: ylims
  # ----- sig_lines a list (for each model) of a list of 
  #        min and max of sig lines (pass in NULL if none)
  # ----- sig_y - where on the y-axis to put sig values line
  #        should be a value for each model
  # ----------------------------------------------------
  # assign variables
  # ----------------------------------------------------
  datpath <- mu_vs_meta_inputs$datpath
  task <- mu_vs_meta_inputs$task
  mods <- mu_vs_meta_inputs$mods
  sub_Ns <- mu_vs_meta_inputs$sub_Ns
  yl <- mu_vs_meta_inputs$yl
  leg_locs <- mu_vs_meta_inputs$leg_locs
  leg_id <- mu_vs_meta_inputs$leg_id
  sig_lines <- mu_vs_meta_inputs$sig_lines
  sig_y <- mu_vs_meta_inputs$sig_y
  # ----------------------------------------------------
  # load data
  # ----------------------------------------------------
  load(paste(datpath, task, "/", task, "stats.RData", sep = ""))
  
  # ----------------------------------------------------
  # get dvs
  # ----------------------------------------------------
  if (task == "SRT" | task == "VSL"){
    # just rename the models to get the dvs, will reset prior to plotting
    ol_mods <- mods
    mods <- c("RM-AN", "LME")
  } 
  ys <- lapply(mods, function(y)
    data.frame(a = do.call(rbind, lapply(sub_Ns, function(x) res[,"stats_fx"][[x]][["mu", y]]))-
                   do.call(rbind, lapply(sub_Ns, function(x) res[,"stats_sig"][[x]][["mu", y]])),
               se = get_pooled(do.call(rbind, lapply(sub_Ns, function(x) res[,"stats_fx"][[x]][["sd", y]]^2)),
                               do.call(rbind, lapply(sub_Ns, function(x) res[,"stats_sig"][[x]][["sd", y]]^2)))))
  
  if (task == "SRT" | task == "VSL") mods <- ol_mods
  
  names(ys) <- mods

  # ----------------------------------------------------
  # do plot
  # ----------------------------------------------------
  plot(x = 1:length(sub_Ns),
       y = t(ys[[mods[[1]]]]["a"]),
       xaxt = "n",
       bty = "n",
       pch = 20,
       cex = 1,
       col = wes_palette("IsleofDogs1")[6],
       xlim = c(0, 21),
       ylim = yl,
       ylab = expression(italic(paste(mu, "diff", sep = " "))),
       xlab = expression(italic("N")),
       cex.lab = 1,
       cex.axis = 1)
  axis(side=1, at=1:length(sub_Ns), labels = sub_Ns)
  
  points(x = jitter(1:length(sub_Ns)),
         y = jitter(t(ys[[mods[[2]]]]["a"])),
         pch = 20,
         cex = 1,
         col = wes_palette("IsleofDogs1")[5])
  
  arrows(x0 = 1:length(sub_Ns),
         y0 = t(ys[[mods[[1]]]]["a"]) - 1.96*(t(ys[[mods[[1]]]]["se"])),
         x1 = 1:length(sub_Ns),
         y1 = t(ys[[mods[[1]]]]["a"]) + 1.96*(t(ys[[mods[[1]]]]["se"])),
         code = 3,
         col = wes_palette("IsleofDogs1")[6],
         angle = 90,
         length = .05)
  
  arrows(x0 = 1:length(sub_Ns),
         y0 = t(ys[[mods[[2]]]]["a"]) - (1.96*t(ys[[mods[[2]]]]["se"])),
         x1 = 1:length(sub_Ns),
         y1 = t(ys[[mods[[2]]]]["a"]) + (1.96*t(ys[[mods[[2]]]]["se"])),
         code = 3,
         col = wes_palette("IsleofDogs1")[5],
         angle = 90,
         length = .05)
  
  if (leg_id){
    leg_cols <- wes_palette("IsleofDogs1")[c(6, 5)]
    legend(x=leg_locs[1], y=leg_locs[2], legend = mods,
           col = leg_cols, pch = 19, bty = "n", cex = 1)
  }
  abline(h=0, lty=2, col="grey48")
  # add signifance lines
  if (!is.null(sig_y)){
    sig_line_cols <- wes_palette("IsleofDogs1")[c(6, 5)]
    names(sig_line_cols) <- mods
    names(sig_y) <- mods
    for (m in mods) {
      lapply(sig_lines[[m]], function(x) segments(x0 = x[1], y0 = sig_y[m],
                                                  x1 = x[2], 
                                                  col = sig_line_cols[m]))
    }
  }
}

# ----------------------------------------------------
# functions to plot mean differences
# ----------------------------------------------------

plot_mean_diff_between_mods <- function(mu_z_inputs) {
  # this function will plot the means and 95% CIs
  # for the diff between mean effect size observations
  # between the two model choices
  # -- mu_z_inputs: a list comprising of:
  # ----- datpath: e.g. "../data/
  # ----- task: e.g. "CC"
  # ----- sub_Ns <- names of subject groups paste(Ns)
  # ----- w width of plot in inches
  # ----- h: height of plot in inches
  # ----- leg_id: legend? TRUE or FALSE
  # ----- leg_locs: e.g. c(5, 20)
  # ----- leg_txt: e.g. c("RM-AN", "LME")
  # ----- yl: ylims

  # ----------------------------------------------------
  # assign variables
  # ----------------------------------------------------
  datpath <- mu_z_inputs$datpath
  task <- mu_z_inputs$task
  mods <- mu_z_inputs$mods
  sub_Ns <- mu_z_inputs$sub_Ns
  yl <- mu_z_inputs$yl
  leg_locs <- mu_z_inputs$leg_locs
  leg_id <- mu_z_inputs$leg_id
  sig_lines <- mu_z_inputs$sig_lines
  sig_y <- mu_z_inputs$sig_y
  # ----------------------------------------------------
  # load data
  # ----------------------------------------------------
  load(paste(datpath, task, "/", task, "stats.RData", sep = ""))

  # ----------------------------------------------------
  # compute ys and std error
  # ----------------------------------------------------
  if (task == "SRT"){
    ol_mods <- mods
    mods <- c("RM-AN", "LME")
  }
  ys <- data.frame(a = do.call(rbind, lapply(sub_Ns, function(x) res[,"stats_fx"][[x]][["mu", mods[1]]]))-
                       do.call(rbind, lapply(sub_Ns, function(x) res[,"stats_fx"][[x]][["mu", mods[2]]])),
                   se = get_pooled(do.call(rbind, lapply(sub_Ns, function(x) res[,"stats_fx"][[x]][["sd", mods[1]]]^2)),
                                   do.call(rbind, lapply(sub_Ns, function(x) res[,"stats_fx"][[x]][["sd", mods[2]]]^2))))
  
  if(task == "SRT") mods <- ol_mods
  # 
  plot(x = 1:length(sub_Ns),
       y = t(ys[["a"]]),
       xaxt = "n",
       bty = "n",
       pch = 20,
       cex = 1,
       col = wes_palette("IsleofDogs1")[2],
       xlim = c(0, 21),
       ylim = yl,
       ylab = expression(italic(paste(mu, "diff", sep = " "))),
       xlab = expression(italic("N")),
       cex.lab = 1,
       cex.axis = 1)
  axis(side = 1, at = 1:length(sub_Ns), labels = sub_Ns)

  arrows(x0 = 1:length(sub_Ns),
         y0 = t(ys[["a"]]) - (1.96*t(ys[["se"]])),
         x1 = 1:length(sub_Ns),
         y1 = t(ys[["a"]]) + (1.96*t(ys[["se"]])),
         code = 3,
         col = wes_palette("IsleofDogs1")[2],
         angle = 90,
         length = .05)
  abline(h=0, lty=2, col="grey48")

  if (!is.null(sig_y)){
  segments(x0 = sig_lines[1], y0 = sig_y,
           x1 = sig_lines[2],
           col = wes_palette("IsleofDogs1")[2])
  }
}

plot_qq_med_vs_best <- function(qq_inputs){
  # this function will plot the qq plot between the
  # median for the field and for the best estimate
  # -- qq_inputs: a list comprising of:
  # ----- datpath: e.g. "../data/
  # ----- task: e.g. "CC"
  # ----- sub_Ns <- names of subject groups paste(Ns)
  # ----- median_N <- median N for the field e.g. 25
  # ----- leg_id: legend? TRUE or FALSE
  # ----- leg_txt: e.g. c("RM-AN", "LME")
  # ----- xl, yl = xlim, ylim
  
  datpath <- qq_inputs$datpath
  task <- qq_inputs$task
  mods <- qq_inputs$mods
  sub_Ns <- qq_inputs$sub_Ns
  median_N <- paste(qq_inputs$median_N)
  yl <- qq_inputs$yl
  xl <- qq_inputs$xl
  leg_locs <- qq_inputs$leg_locs
  leg_id <- qq_inputs$leg_id
  leg_txt <- qq_inputs$leg_txt

  # sig_y <- mu_z_inputs$sig_y
  # ----------------------------------------------------
  # load data and assign variables
  # ----------------------------------------------------
  load(paste(datpath, task, "/", task, "stats.RData", sep = ""))
  medNdens <- res[,"dens_fx"][[median_N]]
  maxNdens <- res[,"dens_fx"][['313']]
  
  # ----------------------------------------------------
  # define function that will extract quantiles from 
  # density object
  # ----------------------------------------------------
  get_quantiles_from_dens <- function(d){
    # d = an object generated by density()

    ps <- d$y/sum(d$y) # convert from density to probability
    # get culmulative probabilities
    culm_ps <- ps[1]
    nps <- length(ps)
    for (i in 2:nps){
      culm_ps[i] <- sum(culm_ps[i-1], ps[i])
    }
    # now get indexes for desired probs
    probs <- seq(0, 1, 1/nps)
    idx <- c()
    for (iP in 1:length(probs)){
      idx[iP] <- which.min(abs(culm_ps - probs[iP]))
    }
    # return y-values for desired probs
    d$x[idx]
  }
  
  medYs <- lapply(names(medNdens), function(x) get_quantiles_from_dens(medNdens[[x]]))
  names(medYs) <- names(medNdens)
  maxYs <- lapply(names(maxNdens), function(x) get_quantiles_from_dens(maxNdens[[x]]))
  names(maxYs) <- names(maxNdens)
  
  # ----------------------------------------------------
  # # plot qq's
  # # ----------------------------------------------------

  # plot max vs median for RM-AN 
  plot(x = medYs[['RM-AN']],
       y = maxYs[['RM-AN']],
       xaxt = "n",
       yaxt = "n",
       ylim = yl,
       xlim = xl,
       bty = "n",
       pch = 20,
       cex = 1,
       col = wes_palette("IsleofDogs1")[6],
       ylab = expression(italic("N"[313])),
       xlab = expression(italic("N"[25])),
       cex.lab = 1,
       cex.axis = 1)
  abline(0, 1, lty=2, col="grey48")
  axis(1, at=xl, labels=xl)
  axis(2, at=yl, labels=yl)
  
  # add LME approach
  points(x = medYs[['LME']],
         y = maxYs[['LME']],
         col = wes_palette("IsleofDogs1")[5])
  
  if (leg_id){
    leg_cols <- wes_palette("IsleofDogs1")[c(6, 5)]
    leg_locs = leg_locs
    legend(x=leg_locs[1], y=leg_locs[2], legend = leg_txt,
           col = leg_cols, pch = 19, bty = "n", cex = 1)
  }
} 
