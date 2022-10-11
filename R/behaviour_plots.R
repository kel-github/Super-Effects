plot_SRT_results <- function(fname) {
  ## plot the random vs seq block effect for mean RTs
  # -----------------------------------------------------
  # read in data
  # -----------------------------------------------------
  dat <- read.csv(fname, header = TRUE)
  # -----------------------------------------------------
  # Create dataframes
  # ------------------------------------------------------
  min.RT <- 200 # in msec
  sd.crit <- 2.5
  
  ffx.dat <- dat %>% filter(Block.No > 2) %>%
    group_by(Subj.No, Block.No.Names) %>%
    filter(Accuracy == 1) %>%
    filter(RT.ms > min.RT) %>%
    filter(RT.ms < (mean(RT.ms) + sd.crit*sd(RT.ms))) %>%
    summarise(RT=mean(RT.ms))
  names(ffx.dat)[names(ffx.dat) == "Block.No.Names"] <- "trialtype"
  
  ffx.dat$RT <- ffx.dat$RT/1000
  
  # ---------------------------------------------------------------
  # plot data
  # ---------------------------------------------------------------
  par(bty = "n")
  vioplot::vioplot(ffx.dat$RT ~ ffx.dat$trialtype, 
                   col = adjustcolor(c(wes_palette("IsleofDogs1")[3], 
                                       wes_palette("IsleofDogs1")[4]), 
                                     alpha = 0.5),
                   xlab = "", ylab = "RT", yaxt = "n")
  axis(side = 1, at = c(1,2), labels = c("R", "S"))
  axis(side = 2, at = c(0.3, 0.8),
       labels = c("0.3", "0.8"),
       las = 1)
  # leg_cols <- wes_palette("IsleofDogs1")[c(3, 4)]
  # legend(2, 1000, legend = c("rand", "rep"),
  #        col = leg_cols, pch = 19, bty = "n", cex = 1)
}


plot_CC_results <- function(fname, type) {
  ## plot the trial type x block interaction for mean RTs
  # -----------------------------------------------------
  # read in data
  # -----------------------------------------------------
  dat <- read.csv(fname, header = TRUE)
  # ------------------------------------------------------
  # create summary dataframe
  # -------------------------------------------------------
  min_RT <- 200 # in msec
  sd_crit <- 2.5
  ffx_dat <- dat %>% mutate(Block.No = rep(c(1:12), each = 24, 
                            length(unique(dat$Subj.No)))) %>%
    group_by(Subj.No, Block.No, Trial.Type.Name) %>%
    filter(Accuracy == 1) %>%
    filter(RT.ms > min_RT) %>%
    filter(RT.ms < (mean(RT.ms) + sd_crit * sd(RT.ms))) %>%
    summarise(RT = mean(RT.ms)/1000)
  subs  <- unique(ffx_dat$Subj.No)
  names(ffx_dat) <- c("sub", "block", "type", "RT")

  # ---------------------------------------------------------------
  # plot data
  # ---------------------------------------------------------------
  if (type == "dist"){
    par(bty = "n")
    vioplot::vioplot(ffx_dat$RT ~ ffx_dat$type*ffx_dat$block, 
                     col = adjustcolor(c(wes_palette("IsleofDogs1")[3], 
                                         wes_palette("IsleofDogs1")[4]), 
                                     alpha = 0.5),
                     xlab = "", ylab = "RT", yaxt = "n",
                     pchMed = 20,
                     colMed = rep(c(wes_palette("IsleofDogs1")[3], wes_palette("IsleofDogs1")[4]), times=12))
    axis(side = 1, at = seq(3.5, 23.5, by = 4), labels = c("2", "4", "6", "8", "10", "12"))
    axis(side = 2, at = c(0.3, 3.8),
         labels = c("0.3", "3.8"),
         las = 1)
    # add legend 
    leg_cols <- wes_palette("IsleofDogs1")[c(3, 4)]
    legend(12, 4, legend = c("novel", "repeat"),
           col = leg_cols, pch = 19, bty = "n", cex = 1)
  } else if (type == "mean"){
    
    ffx_dat <- ffx_dat %>% group_by(sub, block, type) %>%
      summarise(RT = mean(RT)) %>%
      group_by(block, type) %>%
      summarise(mu = mean(RT),
                se = sd(RT) / sqrt(length(RT))) %>% 
      ungroup()
    par(las=1)
    with(ffx_dat, plot(x = block[type == "Novel"],
                       y = mu[type == "Novel"],
                       bty = "n",
                       pch = 20,
                       cex = 1,
                       col = wes_palette("IsleofDogs1")[3],
                       ylim = c(0.95, 1.3),
                       ylab = expression(italic(paste(mu, "RT", sep = " "))),
                       xlab = expression(italic("block")),
                       cex.lab = 1,
                       cex.axis = 1))
     with(ffx_dat, points(x = block[type == "Repeated"],
                          y = mu[type == "Repeated"],
                          pch = 20,
                          cex = 1,
                          col = wes_palette("IsleofDogs1")[4]))
   # now add error bars
     with(ffx_dat, arrows(x0 = block[type == "Novel"], 
                          y0 = mu[type == "Novel"] - se[type == "Novel"],
                          x1 = block[type == "Novel"],
                          y1 = mu[type == "Novel"] + se[type == "Novel"],
                          code = 3,
                          col = wes_palette("IsleofDogs1")[3],
                          angle = 90,
                          length = .025))
     with(ffx_dat, arrows(x0 = block[type == "Repeated"], 
                          y0 = mu[type == "Repeated"] - se[type == "Repeated"],
                          x1 = block[type == "Repeated"],
                          y1 = mu[type == "Repeated"] + se[type == "Repeated"],
                          code = 3,
                          col = wes_palette("IsleofDogs1")[4],
                          angle = 90,
                          length = .025))

    # leg_cols <- wes_palette("IsleofDogs1")[c(3, 4)]
    # legend(2, 1000, legend = c("novel", "repeat"),
    #        col = leg_cols, pch = 19, bty = "n", cex = 1)
  }
}

plot_MT_results <- function(fname) {

  library(vioplot)
  # ----------------------------------------------------------------------------------------------------
  # load data and wrangle into tidy form 
  # ----------------------------------------------------------------------------------------------------
  
  dat <- read.csv(fname, header=TRUE)
  # ----------------------------------------------------------------------------------------------------
  # Create dataframes 
  # ----------------------------------------------------------------------------------------------------
  # Create a summary of the data for ffx and rfx modelling
  min.RT <- .200 # in sec
  sd.crit <- 2.5
  
  rfx.dat <- dat %>% filter(Overall.Accuracy == 1) %>%
    select(c('Subj.No', 'Trial.Type.Name', 'Task.1.RT.Sec', 'Task.2.RT.Sec', 'Task.1.Response', 'Task.2.Response'))  %>%
    pivot_longer(c('Task.1.RT.Sec', 'Task.2.RT.Sec'), names_to = "task", values_to="RT") %>%
    drop_na()
  rfx.dat$task[rfx.dat$Trial.Type.Name == 'single_auditory'] = 'sound'
  rfx.dat$task[rfx.dat$Trial.Type.Name == 'single_visual'] = 'vis'
  rfx.dat$task[rfx.dat$Trial.Type.Name == 'dual_task' & rfx.dat$task == 'Task.1.RT.Sec'] = 'vis'
  rfx.dat$task[rfx.dat$Trial.Type.Name == 'dual_task' & rfx.dat$task == 'Task.2.RT.Sec'] = 'sound'
  
  rfx.dat <- rfx.dat %>% mutate(trialtype = fct_recode(Trial.Type.Name,
                                                       'single' = 'single_auditory',
                                                       'single' = 'single_visual',
                                                       'dual' = 'dual_task')) 
  rfx.dat$task.stim <- NA
  rfx.dat$task.stim[rfx.dat$trialtype == "single"] = rfx.dat$Task.1.Response[rfx.dat$trialtype == "single"]
  rfx.dat$task.stim[rfx.dat$trialtype == "dual" & rfx.dat$task == "vis"] = rfx.dat$Task.1.Response[rfx.dat$trialtype == "dual" & rfx.dat$task == "vis"]
  rfx.dat$task.stim[rfx.dat$trialtype == "dual" & rfx.dat$task == "sound"] = rfx.dat$Task.2.Response[rfx.dat$trialtype == "dual" & rfx.dat$task == "sound"]
  
  ffx.dat <- rfx.dat %>% select(-c("Trial.Type.Name", "Task.1.Response", "Task.2.Response")) %>%
    group_by(Subj.No, task, trialtype) %>%
    filter(RT > min.RT) %>%
    filter(RT < (mean(RT)+sd.crit*sd(RT))) %>%
    summarise(RT = mean(RT))
  names(ffx.dat)[names(ffx.dat) == "Subj.No"] = "sub"
  
  # reorder factors
  ffx.dat$task <- factor(ffx.dat$task, levels = c("vis", "sound"))
  ffx.dat$trialtype <- factor(ffx.dat$trialtype, levels = c("single", "dual"))
  
  
  # ---------------------------------------------------------------
  # plot data
  # ---------------------------------------------------------------
  par(bty = "n")
  vioplot(ffx.dat$RT ~ ffx.dat$task*ffx.dat$trialtype, 
          col = adjustcolor(rep(c(wes_palette("IsleofDogs1")[3], wes_palette("IsleofDogs1")[4]), times = 2), alpha = 0.5),
          xlab = "", ylab = "RT", yaxt = "n")
  axis(side = 1, at = 1:4, labels = c("VS", "SoS", "VM", "SoM"))
  axis(side = 2, at = c(min(ffx.dat$RT), max(ffx.dat$RT)),
                 labels = c(sprintf("%1.2f", min(ffx.dat$RT)),
                            sprintf("%1.2f", max(ffx.dat$RT))),
       las = 1)
}

plot_AB_results <- function(fname) {
  # ----------------------------------------------------------------------------------------------------
  # read in data
  # ----------------------------------------------------------------------------------------------------  
  dat <- read.csv(fname, header=TRUE)
  
  # ----------------------------------------------------------------------------------------------------
  # create summary dataframe
  # ----------------------------------------------------------------------------------------------------
  
  # Create a summary of the data for fixed fx modelling
  ffx_dat <- dat %>% group_by(Subj.No, Trial.Type.Name) %>%
    summarise(T1.Accuracy=mean(T1.Accuracy),
              T2T1.Accuracy=mean(T2T1.Accuracy)) %>%
    pivot_longer(cols=c("T1.Accuracy", "T2T1.Accuracy"),
                 names_to = "condition",
                 values_to = "accuracy") %>%
    ungroup()
  names(ffx_dat)[names(ffx_dat)=="Trial.Type.Name"] = "lag"
  ffx_dat$lag <- fct_recode(ffx_dat$lag, "2" = "lag_2", "3" = "lag_3", "5" = "lag_5", "7" = "lag_7")
  ffx_dat$condition <- fct_recode(ffx_dat$condition, "T1" = "T1.Accuracy", "T2gT1" = "T2T1.Accuracy")
  
  # ---------------------------------------------------------------
  # plot data
  # ---------------------------------------------------------------
  par(bty = "n")
  vioplot::vioplot(ffx_dat$accuracy ~ ffx_dat$condition*ffx_dat$lag, 
                   col = adjustcolor(rep(c(wes_palette("IsleofDogs1")[3], wes_palette("IsleofDogs1")[4]), 
                                           times = 4), alpha = 0.5),
                   xlab = "lag", ylab = "acc", yaxt = "n")
  axis(side = 1, at = c(1.5, 3.5, 5.5, 7.5), labels = c("2", "3", "5", "7"))
  axis(side = 2, at = c(0, 1),
       labels = c("0", "1"),
       las = 1)
  
  leg_cols <- wes_palette("IsleofDogs1")[c(3, 4)]
  legend(1.9, .3, legend = c("T1", "T2gT1"),
         col = leg_cols, pch = 19, bty = "n", cex = 1)
}

plot_VSL_results <- function(fname){
  
  # ----------------------------------------------------------------------------------------------------
  # read in data and wrangle
  # ----------------------------------------------------------------------------------------------------  
  dat <- read.csv(fname, header=TRUE) 
  
  
  # data frame contains TRUE ordering
  prev.dat <- dat %>% select(Subj.No, Trial.No, Response, Target.Order, Accuracy)
  prev.dat$Response <- as.factor(prev.dat$Response)
  prev.dat$Target.Order <- as.factor(prev.dat$Target.Order)
  
  prev.dat <- prev.dat %>% mutate(Response = recode(Response,
                                                    "122" = "Novel",
                                                    "109" = "Repeat"),
                                  Target.Order = recode(Target.Order,
                                                        "1" = "Novel",
                                                        "2" = "Repeat")) %>%
                          group_by(Subj.No) %>%
                          summarise(acc = mean(Response==Target.Order)) %>%
                          ungroup() 
  
  
  with(prev.dat, hist(x = acc, freq=FALSE, breaks = 20,
                      bty = "n", 
                      col = wes_palette("IsleofDogs1")[3],
                      xlim = c(0, 1),
                      ylab = expression(italic(paste("d"))),
                      xlab = expression(italic("acc")),
                      cex.lab = 1,
                      cex.axis = 1,
                      xaxt = "n",
                      yaxt = "n",
                      main = ""))
  axis(side=1, at = c(0, .5, 1), labels = c("0", ".5", "1"))
  axis(side=2, at = c(0, 4), labels = c("0", "4"))
  abline(v = 0.5, lty = 2, col = "grey48")
}