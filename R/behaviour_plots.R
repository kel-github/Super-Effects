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