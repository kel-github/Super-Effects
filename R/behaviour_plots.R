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
  
  grp_dat <- ffx.dat %>% group_by(trialtype) %>%
                         mutate(RT = RT/1000) %>%
                         summarise(mu = mean(RT),
                                   N = length(RT),
                                   se = sd(RT)/sqrt(N)) %>%
                         ungroup()
  
  # ---------------------------------------------------------------
  # plot data
  # ---------------------------------------------------------------
  with(grp_dat, plot(x = c(1,2),
                     y = mu,
                     bty = "n",
                     pch = 20,
                     cex = 1,
                     col = c(wes_palette("IsleofDogs1")[3], wes_palette("IsleofDogs1")[4]),
                     ylim = c(.2, .450),
                     xlim = c(0.5, 2.5),
                     ylab = expression(italic(paste(mu, "RT", sep = " "))),
                     xlab = expression(italic("block")),
                     cex.lab = 1,
                     cex.axis = 1,
                     xaxt = "n"))
  axis(side = 1, at = c(1, 2), labels = c("Ran", "Rep"))
  
  with(grp_dat, arrows(x0 = c(1,2), 
                       y0 = c(mu[trialtype == "Random Block"] - (1.96*se[trialtype == "Random Block"]),
                              mu[trialtype == "Sequence Block"] - (1.96*se[trialtype == "Sequence Block"])),  
                       x1 = c(1,2), 
                       y1 = c(mu[trialtype == "Random Block"] + (1.96*se[trialtype == "Random Block"]),
                              mu[trialtype == "Sequence Block"] + (1.96*se[trialtype == "Sequence Block"])),
                       code = 3,
                       col = c(wes_palette("IsleofDogs1")[3], wes_palette("IsleofDogs1")[4]),
                       angle = 90,
                       length = .025))
  
  # leg_cols <- wes_palette("IsleofDogs1")[c(3, 4)]
  # legend(2, 1000, legend = c("rand", "rep"),
  #        col = leg_cols, pch = 19, bty = "n", cex = 1)
}


plot_CC_results <- function(fname) {
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
    summarise(RT = mean(RT.ms))
  subs  <- unique(ffx_dat$Subj.No)
  names(ffx_dat) <- c("sub", "block", "type", "RT")

  # Create a summary of the data for fixed fx depiction
  ffx_dat <- ffx_dat %>% group_by(sub, block, type) %>%
                 summarise(RT = mean(RT)) %>%
                 group_by(block, type) %>%
                 summarise(mu = mean(RT),
                 se = sd(RT) / sqrt(length(RT))) %>% 
                 ungroup()

  # ---------------------------------------------------------------
  # plot data
  # ---------------------------------------------------------------
  with(ffx_dat, plot(x = block[type == "Novel"],
                     y = mu[type == "Novel"],
                     bty = "n",
                     pch = 20,
                     cex = 1,
                     col = wes_palette("IsleofDogs1")[3],
                     ylim = c(800, max(mu)),
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

    leg_cols <- wes_palette("IsleofDogs1")[c(3, 4)]
    legend(2, 1000, legend = c("novel", "repeat"),
           col = leg_cols, pch = 19, bty = "n", cex = 1)
}

plot_MT_results <- function(fname) {

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
  
  grp_sum <- ffx.dat %>% group_by(task, trialtype) %>%
                         summarise(mu = mean(RT),
                                   N = length(RT),
                                   se = sd(RT)/N) %>% 
                         ungroup()
  
  # ---------------------------------------------------------------
  # plot data
  # ---------------------------------------------------------------
  ys <- with(grp_sum, c(mu[task == "sound" & trialtype == "single"], 
                        mu[task == "vis" & trialtype == "single"],
                        mu[task == "sound" & trialtype == "dual"],
                        mu[task == "vis" & trialtype == "dual"]))
  upper = ys + 1.96*with(grp_sum, c(se[task == "sound" & trialtype == "single"], 
                                    se[task == "vis" & trialtype == "single"],
                                    se[task == "sound" & trialtype == "dual"],
                                    se[task == "vis" & trialtype == "dual"]))
  lower = ys - 1.96*with(grp_sum, c(se[task == "sound" & trialtype == "single"], 
                                    se[task == "vis" & trialtype == "single"],
                                    se[task == "sound" & trialtype == "dual"],
                                    se[task == "vis" & trialtype == "dual"]))
  plot(x = c(1,2,3,4),
       y = ys,
       bty = "n",
       pch = 20,
       cex = 1,
       col = rep(c(wes_palette("IsleofDogs1")[3], wes_palette("IsleofDogs1")[4]), times = 2),
       xlim = c(.5, 4.5),
       ylim = c(0.4, 1.2),
       ylab = expression(italic(paste(mu, "RT", sep = " "))),
       xlab = expression(italic("task")),
       cex.lab = 1,
       cex.axis = 1,
       xaxt = "n")
  axis(side = 1, at = c(1, 2, 3, 4), labels = c("S", "S", "M", "M"))
  # now add error bars
  arrows(x0 = c(1,2,3,4), 
         y0 = lower,
         x1 = c(1,2,3,4), 
         y1 = upper,
         code = 3,
         col = rep(c(wes_palette("IsleofDogs1")[3], wes_palette("IsleofDogs1")[4]), times = 2),
         angle = 90,
         length = .05)
  leg_cols <- wes_palette("IsleofDogs1")[c(3, 4)]
  legend(3, .65, legend = c("A", "V"),
         col = leg_cols, pch = 19, bty = "n", cex = 1)
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
    group_by(Trial.Type.Name, condition) %>%
    summarise(acc=mean(accuracy),
              se=sd(accuracy)/sqrt(length(accuracy))) %>% 
    ungroup()
  names(ffx_dat)[names(ffx_dat)=="Trial.Type.Name"] = "lag"
  ffx_dat$lag <- fct_recode(ffx_dat$lag, "2" = "lag_2", "3" = "lag_3", "5" = "lag_5", "7" = "lag_7")
  ffx_dat$condition <- fct_recode(ffx_dat$condition, "T1" = "T1.Accuracy", "T2gT1" = "T2T1.Accuracy")
  
  # ---------------------------------------------------------------
  # plot data
  # ---------------------------------------------------------------
  with(ffx_dat, plot(x = c(2,3,5,7),
                     y = acc[condition == "T1"],
                     bty = "n",
                     pch = 20,
                     cex = 1,
                     col = wes_palette("IsleofDogs1")[3],
                     xlim = c(1.5, 7.5),
                     ylim = c(0.4, 1),
                     ylab = expression(italic(paste(mu, "acc", sep = " "))),
                     xlab = expression(italic("lag")),
                     cex.lab = 1,
                     cex.axis = 1))
  with(ffx_dat, points(x = c(2,3,5,7),
                       y = acc[condition == "T2gT1"],
                       pch = 20,
                       cex = 1,
                       col = wes_palette("IsleofDogs1")[4]))
  # now add error bars
  with(ffx_dat, arrows(x0 = c(2,3,5,7), 
                       y0 = acc[condition == "T1"] - se[condition == "T1"],
                       x1 = c(2,3,5,7),
                       y1 = acc[condition == "T1"] + se[condition == "T1"],
                       code = 3,
                       col = wes_palette("IsleofDogs1")[3],
                       angle = 90,
                       length = .05))
  with(ffx_dat, arrows(x0 = c(2,3,5,7),
                       y0 = acc[condition == "T2gT1"] - se[condition == "T2gT1"],
                       x1 = c(2,3,5,7),
                       y1 = acc[condition == "T2gT1"] + se[condition == "T2gT1"],
                       code = 3,
                       col = wes_palette("IsleofDogs1")[4],
                       angle = 90,
                       length = .05))
  
  leg_cols <- wes_palette("IsleofDogs1")[c(3, 4)]
  legend(5, .6, legend = c("T1", "T2gT1"),
         col = leg_cols, pch = 19, bty = "n", cex = 1)
}