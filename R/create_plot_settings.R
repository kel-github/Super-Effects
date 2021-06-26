task <- "AB"

# ----------------------------------------------------
# behavioural data
# ----------------------------------------------------
fname <- "../data/total_of_313_subs_AB_task_trial_level_data.csv"
sub_Ns <- paste(round(exp(seq(log(13), log(313), length.out = 20))))
# ----------------------------------------------------
# density variables
# ----------------------------------------------------
fx <- list(datpath = "../data/",
           task = "AB",
           jmax = 2,
           dv = "dens_fx",
           sel_n = paste(c(25, 59, 136, 313)),
           w = 1.96,
           h = 2.36 * 2,
           xlabs = c(expression(eta[p]^2), expression("r"^2)),
           xl = c(0, 1),
           max_idx = c(20, 20),
           leg_id = 2,
           leg_locs = c(0.5, 20),
           figlabel = "B",
           figlabelon = TRUE)

p <- list(datpath = "../data/",
          task = "AB",
          jmax = 2,
          dv = "dens_p",
          sel_n = paste(c(25, 59, 136, 313)),
          w = 1.96,
          h = 2.36 * 2,
          xlabs = c("p", "p"),
          xl = c(-60, 0),
          max_idx = c(20, 20),
          leg_id = 1,
          leg_locs = c(-60, 0.4),
          figlabel = "A",
          figlabelon = TRUE)

# ----------------------------------------------------
# KL divergence
# ----------------------------------------------------

kl <- list(datpath = "../data/",
           task = "AB",
           dv = "dens_fx",
           ratio_type = "KL",
           origin = "313",
           sub_Ns = sub_Ns,
           w = 1.96,
           h = 2.36,
           leg_id = FALSE,
           leg_locs = c(5, 20),
           leg_txt = c("RM-AN", "LME"),
           ylabel = expression(italic("KL p||q")),
           yl = NULL)

# ----------------------------------------------------
# fx sz ratio between models
# ----------------------------------------------------
           
model_rats <- list(datpath = "../data/",
              task = "AB",
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
              mods = c("RM-AN", "LME"),
              yl = NULL)

# ----------------------------------------------------
# meta-analytic vs observed fx sz ratio
# ----------------------------------------------------
sig <- list(datpath = "../data/",
            task = "AB",
            dv = "stats_sig",
            ratio_type = "stats_sig",
            origin = "",
            sub_Ns = sub_Ns,
            w = 1.96,
            h = 2.36,
            leg_id = TRUE,
            leg_locs = c(5, 0.5),
            leg_txt = c("RM-AN", "LME"),
            ylabel = expression(italic("meta / sim")),
            mods = c("RM-AN", "LME"),
            yl = c(0, 2))

save(task, fname, fx, p, kl, model_rats, sig, 
     file = paste("../data/", task, "/",
                  task, "_plot_settings.RData", sep = ""))
