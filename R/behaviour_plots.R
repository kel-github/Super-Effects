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