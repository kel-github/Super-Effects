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
  if (length(xl) == 2){
    xl <- matrix(c(xl, xl), ncol = 2)
  } else {
    xl <- matrix(xl, ncol = 2)
  }
  max_idx <- inputs4plot$max_idx
  leg_id <- inputs4plot$leg_id
  leg_locs <- inputs4plot$leg_locs
  figlabel <- inputs4plot$figlabel
  figlabelon <- inputs4plot$figlabelon
  mods <- c("RM-AN", "LME")
  imm <- inputs4plot$imm
  
  palette_choice <- wes_palette("IsleofDogs1")[c(1, 2, 3, 5)]
  # ----------------------------------------------------
  # load data
  # ----------------------------------------------------
  if (imm) {
    load(paste(datpath, task, "/", "IMM", task, "stats.RData", sep = ""))
  } else {
    load(paste(datpath, task, "/", task, "stats.RData", sep = ""))
  }
  
  for (j in 1:jmax) {
    plot(res[sel_n[length(sel_n)], ][[dv]][[mods[j]]],
         col = adjustcolor(palette_choice[length(sel_n)], alpha.f = 0.75),
         lwd = 2,
         ylim = c(0,
                  max(res[max_idx[j], ][[dv]][[mods[j]]]["y"]$y)),
         xlim = xl[,j],
         main = " ", ylab = "density",
         xlab = xlabs[j],
         bty = "n",
         cex.lab = 1,
         cex.axis = 1)
    # polygon(res[sel_n[length(sel_n)], ][[dv]][[mods[j]]],
    #         col = adjustcolor(wes_palette("IsleofDogs1")[length(sel_n)],
    #         alpha.f = 0.5),
    #         border = NA)
    for (i in c((length(sel_n)-1):1)) {
      lines(res[sel_n[i], ][[dv]][[mods[j]]],
            col = adjustcolor(palette_choice[i], alpha.f = 0.75),
            lwd = 2)
      # polygon(res[sel_n[i], ][[dv]][[mods[j]]],
      #       col = adjustcolor(wes_palette("IsleofDogs1")[i],
      #       alpha.f = 0.5),
      #       border = NA)
    }
    # if a p statistic plot, add the criteria for significance
    if (dv == "dens_p") {
      abline(v = qnorm(.05), col = "#161616",
             lty = 2, lwd = 1)
    }
    # do you want a legend, and if so where to put it?
    if (j == leg_id) {
      leg_cols <-adjustcolor(palette_choice, alpha = 0.75)
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
  mods <- rat_inputs$mods
  # get models from which to return results
  #mods <- unique(colnames(res[, dv][[origin]]))
  #mods <- unique(names(res[, dv][[origin]]))
  # rem.infin <- function(res, dv){
  #   # a quick function to replace infinite values with a really 
  #   # low value
  #   for (i in 1:length(res[,dv])){
  #     res[,dv][[i]][,"RM-AN"]$qs[is.infinite(res[,dv][[i]][,"RM-AN"]$qs)] <- qnorm(1e-300)
  #   }
  #   res
  # }
  # res <- rem.infin(res,dv)
  base <- do.call(cbind, lapply(mods,
                                function(x)
                                  abs(diff(res[, dv][[origin]][, x][["qs"]]))))
  colnames(base) <- mods # just labeling for sanity checks
  ratios4plotting <- do.call(rbind,
                             lapply(sub_Ns, 
                                    function(x) lapply(mods, function(y)
                                      abs(diff(res[, dv][[x]][,y][["qs"]]) / base[,y]))))
  rownames(ratios4plotting) <- sub_Ns
  colnames(ratios4plotting) <- mods
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
  xvals <- rat_inputs$xvals
  mods <- unique(unique(colnames(res[, "stats_fx"][[origin]])))
  # first calculate the approximating function for the all the densities
  d_funcs <- lapply(sub_Ns, function(j)
    lapply(mods, function(z)
      approxfun(x = res[, dv][[j]][[z]]$x,
                y=res[ , dv][[j]][[z]]$y/length(res[ , dv][[j]][[z]]$y),
                rule = 1)))
  names(d_funcs) <- sub_Ns
  for (i in sub_Ns) names(d_funcs[[i]]) <- mods
  # then select a series of bin widths, or xs
  if (is.null(xvals)) {
    xvals <- list(seq(0, 1, by = .001), seq(0,1, by = .001))
    
  } else { 
    xvals <- matrix(xvals, ncol = 2)
    sa <- seq(xvals[1,1], xvals[2,1], by = diff(xvals[,1])/1000)
    sb <- seq(xvals[1,2], xvals[2,2], by = diff(xvals[,2])/1000)
    xvals <- list(sa, sb)
  }
  names(xvals) <- c("RM-AN", "LME")
  # cos our effect sizes are positive correlations
  # do the same for the second distribution
  # calculate the KL divergence between the two
  kl_func <- function(pf, qf, x) {
    # -- pf: function for p - P IS THE P DIST
    # -- pq: function for q - Q IS THE APPROXIMATING DIST
    # -- x: values of x
    p <- pf(x) + (rnorm(length(x), 0, 1e-06)^2)
    q <- qf(x) + (rnorm(length(x), 0, 1e-06)^2)
    kl <- p * log2(p / q)
    sum(kl, na.rm = T)
  }
  kl4plotting <- lapply(sub_Ns, function(x)
    lapply(mods,
           function(y)
             kl_func(pf = d_funcs[[origin]][[y]],
                     qf = d_funcs[[x]][[y]],
                     x = xvals[[y]])))
  kl4plotting <- do.call(rbind, kl4plotting)
  colnames(kl4plotting) <- mods
  rownames(kl4plotting) <- sub_Ns
  kl4plotting
}


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
  mods <- rat_inputs$mods # model ratios to plot
  imm <- rat_inputs$imm # TRUE or FALSE for sampling type
  if (ratio_type == "origin" | ratio_type == "KL") origin <- rat_inputs$origin
  if (dv == "stats_p") yl <- NULL # a hack setting because I don't want to go back and resave all the 
  # plot settings and I want yl to be set by the data range rather than by a pre-defined value
  # ----------------------------------------------------
  # load data
  # ----------------------------------------------------
  if (imm) {
    load(paste(datpath, task, "/", "IMM", task, "stats.RData", sep = ""))
  } else {
    load(paste(datpath, task, "/", task, "stats.RData", sep = ""))
  }
  
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
  mnmes <- colnames(ratios)
  ratios <- matrix(unlist(ratios), ncol=length(ratios)/length(sub_Ns)) 
  colnames(ratios) <- mnmes
  
  if (length(mods) > 1) {
    if (is.null(yl)){
      if (ncol(ratios) != 2){
        if (is.list(ratios)){
        nuyl = c(0,
                 max(cbind(do.call(cbind, ratios[, "RM-AN"]),
                           do.call(cbind, ratios[, "LME"]))))
        } else {
          nuyl = c(0, max(ratios))
        } 
      } else {
        nuyl = c(0, max(ratios))
        rownames(ratios) <- sub_Ns
      }
    } else {
      nuyl = yl
    }
    if (ratio_type == "KL" | dv == "stats_p" ){
      xs = as.numeric(sub_Ns)
    } else {
      xs = 1:length(rownames(ratios))
    }
    plot(x = xs, y = ratios[, "RM-AN"],
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
    if (ratio_type == "KL" | dv == "stats_p"){
      axis(1, at = as.numeric(sub_Ns),
           labels = as.numeric(sub_Ns),
           cex.lab = 1,
           cex.axis = 1)
    } else {
      axis(1, at = seq(1, 20, 2),
           labels = sub_Ns[seq(1, 20, 2)],
           cex.lab = 1,
           cex.axis = 1)
    }
    if (task == "SRT" & ratio_type == "KL"){
      points(x = xs, y = ratios[, "LME"],
             pch = 19, col = wes_palette("IsleofDogs1")[5])
    } else {
      if (ratio_type == "KL" | dv == "stats_p"){
        xs = as.numeric(sub_Ns)
      } else {
        xs = 1:length(rownames(ratios))
      }
      lines(x = xs, y = ratios[, "LME"],
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
  } else if (ncol(ratios) == 1 | length(mods == 1)) {
    if (is.null(yl)){
      nuyl = c(0, 
               max(ratios[,1]+0.5))
    } else {
      nuyl = yl
    }
    if (ratio_type == "KL" | dv == "stats_p"){
      xs = as.numeric(sub_Ns)
    } else {
      xs = 1:length(rownames(ratios))
    }
    plot(x = xs, y = ratios[, 1], 
         xaxt = "n",
         yaxt = "n",
         bty = "n",
         type = "l", lty = 1,
         lwd = 2,
         col = wes_palette("IsleofDogs1")[4],
         ylim = nuyl,
         ylab = ylabel,
         xlab = expression(italic("N")),
         cex.lab = 1,
         cex.axis = 1)
    if (ratio_type == "KL" | dv == "stats_p"){
      axis(1, at = as.numeric(sub_Ns), 
           labels = sub_Ns, 
           cex.lab = 1,
           cex.axis = 1)
      if (dv == "stats_p"){
        ylabels <- sprintf("%.1e", nuyl)
      } else {
        ylabels <- sprintf("%.2f", nuyl)
      }
      axis(2, at = nuyl,
           labels = ylabels,
           cex.lab = 1,
           cex.axis = 1)
    } else {
      axis(1, at = seq(1, 20, 2), 
           labels = sub_Ns[seq(1, 20, 2)], 
           cex.lab = 1,
           cex.axis = 1)
      axis(2, at = nuyl,
           labels = nuyl,
           cex.lab = 1,
           cex.axis = 1)
    }
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
  leg_txt <- mu_vs_meta_inputs$leg_txt
  leg_id <- mu_vs_meta_inputs$leg_id
  sig_lines <- mu_vs_meta_inputs$sig_lines
  sig_y <- mu_vs_meta_inputs$sig_y
  imm <- mu_vs_meta_inputs$imm
  # ----------------------------------------------------
  # load data
  # ----------------------------------------------------
  if (imm) {
    load(paste(datpath, task, "/", "IMM", task, "stats.RData", sep = ""))
  } else {
    load(paste(datpath, task, "/", task, "stats.RData", sep = ""))
  }
  
  # ----------------------------------------------------
  # get dvs
  # ----------------------------------------------------
  if (task == "VSL"){
    # just rename the models to get the dvs, will reset prior to plotting
    ol_mods <- mods
    mods <- c("RM-AN", "LME")
  } else if (task == "SRT"){
    ol_mods <- mods
    mods <- c("RM-AN")
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
  
  arrows(x0 = 1:length(sub_Ns),
         y0 = t(ys[[mods[[1]]]]["a"]) - 1.96*(t(ys[[mods[[1]]]]["se"])),
         x1 = 1:length(sub_Ns),
         y1 = t(ys[[mods[[1]]]]["a"]) + 1.96*(t(ys[[mods[[1]]]]["se"])),
         code = 3,
         col = wes_palette("IsleofDogs1")[6],
         angle = 90,
         length = .05)
  
  if (length(mods) > 1){
  points(x = jitter(1:length(sub_Ns)),
         y = jitter(t(ys[[mods[[2]]]]["a"])),
         pch = 20,
         cex = 1,
         col = wes_palette("IsleofDogs1")[5])

  arrows(x0 = 1:length(sub_Ns),
         y0 = t(ys[[mods[[2]]]]["a"]) - (1.96*t(ys[[mods[[2]]]]["se"])),
         x1 = 1:length(sub_Ns),
         y1 = t(ys[[mods[[2]]]]["a"]) + (1.96*t(ys[[mods[[2]]]]["se"])),
         code = 3,
         col = wes_palette("IsleofDogs1")[5],
         angle = 90,
         length = .05)
  }
  
  if (leg_id){
    leg_cols <- wes_palette("IsleofDogs1")[c(6, 5)]
    legend(x=leg_locs[1], y=leg_locs[2], legend = leg_txt,
           col = leg_cols, pch = 19, bty = "n", cex = 1)
  }
  abline(h=0, lty=2, col= adjustcolor("grey48", alpha = 0.3))
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
  modn <- qq_inputs$modn
  imm <- qq_inputs$imm # TRUE or FALSE

  # sig_y <- mu_z_inputs$sig_y
  # ----------------------------------------------------
  # load data and assign variables
  # ----------------------------------------------------
  
  if (imm) {
    load(paste(datpath, task, "/", "IMM", task, "stats.RData", sep = ""))
  } else {
    load(paste(datpath, task, "/", task, "stats.RData", sep = ""))
  }
  
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
    if (task != "SRT") probs <- seq(0, 1, 1/nps)
    if (task == "SRT") probs <- seq(0, 3, 1/nps)
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
       xlab = expression(italic("N"[36])),
       cex.lab = 1,
       cex.axis = 1)
  abline(0, 1, lty=2, col="grey48")
  axis(1, at=xl, labels=xl)
  axis(2, at=yl, labels=yl)
  
  if (modn > 1){
 # add LME approach
  points(x = medYs[['LME']],
         y = maxYs[['LME']],
         col = wes_palette("IsleofDogs1")[5])

  if (leg_id) {
    leg_cols <- wes_palette("IsleofDogs1")[c(6, 5)]
    leg_locs = leg_locs
    legend(x = leg_locs[1], y = leg_locs[2], legend = leg_txt,
           col = leg_cols, pch = 19, bty = "n", cex = 1)
    }
  }
} 

mean_med_mode <- function(plot_set){
  # hacky function to plot mean, median and mode from selected n, using fx settngs
  # ----------------------------------------------------
  # settings - see fx for definitions
  # ----------------------------------------------------
  # imm = TRUE if immediate samples
  # datpath
  # task
  # sel_n
  # xlabs

  datpath = plot_set$datpath
  task = plot_set$task
  imm = plot_set$imm
  sel_n = plot_set$sel_n
  xlabs = plot_set$xlabs
  # ----------------------------------------------------
  # load data and assign variables
  # ----------------------------------------------------
  palette_choice <- wes_palette("IsleofDogs1")[c(1, 2, 3, 5)] # for colours
  if (imm) {
    load(paste(datpath, task, "/", "IMM", task, "stats.RData", sep = ""))
  } else {
    load(paste(datpath, task, "/", task, "stats.RData", sep = ""))
  }
# first get mode values    
  get_mode_idx <- function(res, N, mod){
    which(abs(max(res[N, "dens_fx"][[1]][[mod]]$y)-res[N, "dens_fx"][[1]][[mod]]$y) 
          == min(abs(max(res[N, "dens_fx"][[1]][[mod]]$y)-res[N, "dens_fx"][[1]][[mod]]$y)))
  }
  mode_idxs <- list(RM = sapply(sel_n, get_mode_idx, res=res, mod = "RM-AN"),
                    LME = sapply(sel_n, get_mode_idx, res=res, mod = "LME"))
  get_mode <- function(res, N, mod, idx){
    res[N, "dens_fx"][[1]][[mod]]$x[idx[N]]
  }
  modes <- list(RM = sapply(sel_n, get_mode, res=res, mod="RM-AN", idx=mode_idxs$RM),
                LME = sapply(sel_n, get_mode, res=res, mod="LME", idx=mode_idxs$LME))
  
  # now get medians
  meds <- list(RM = do.call(cbind, sapply(sel_n, function(x) res[x, "stats_fx"][[1]]["med", "RM-AN"])),
               LME = do.call(cbind, sapply(sel_n, function(x) res[x, "stats_fx"][[1]]["med", "LME"])))
  
  mus <- list(RM = do.call(cbind, sapply(sel_n, function(x) res[x, "stats_fx"][[1]]["mu", "RM-AN"])),
              LME = do.call(cbind, sapply(sel_n, function(x) res[x, "stats_fx"][[1]]["mu", "LME"])))

  # ----------------------------------------------------  
  # now can begin plotting
  # ----------------------------------------------------
  # PLOT 1st FX
  ylims <- c(.01, .07)
  plot(as.numeric(sel_n), modes$RM, 
       pch = 17, col = palette_choice,
       xaxt = "n", yaxt = "n", ylim = ylims,
       xlab = "N", ylab = xlabs[1], bty = "n")
  points(as.numeric(sel_n), meds$RM,
         pch = 19, col = palette_choice)
  points(as.numeric(sel_n), mus$RM,
         pch = 15, col = palette_choice)
  axis(1, at=as.numeric(sel_n), labels = sel_n)
  axis(2, at=ylims, labels=as.character(ylims))
  
  # add legend here
  legend(x=250, y=ylims[2],
         legend = c("mode", "med", "mean"),
         pch = c(17, 19, 15), bty="n")
  fig_label("C", cex = 2)
  
  # PLOT 2nd FX
  ylims <- c(.01, .14)
  plot(as.numeric(sel_n), modes$LME, 
       pch = 17, col = palette_choice,
       xaxt = "n", yaxt = "n", ylim = ylims,
       xlab = "N", ylab = xlabs[1], bty = "n")
  points(as.numeric(sel_n), meds$LME,
         pch = 19, col = palette_choice)
  points(as.numeric(sel_n), mus$LME,
         pch = 15, col = palette_choice)
  axis(1, at=as.numeric(sel_n), labels = sel_n)
  axis(2, at=ylims, labels=as.character(ylims))
}
