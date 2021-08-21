rm(list=ls())

task <- "SRT"
datpath <- "../data/"
medN <- "36"
mods <- c("t", "LME")
# ----------------------------------------------------
# behavioural data
# ----------------------------------------------------
fname <- "../data/total_of_313_subs_SRT_task_trial_level_data.csv"
sub_Ns <- paste(round(exp(seq(log(13), log(313), length.out = 20))))
# ----------------------------------------------------
# density variables
# ----------------------------------------------------
fx <- list(datpath = datpath,
           task = task,
           jmax = 2,
           dv = "dens_fx",
           sel_n = paste(c(25, 59, 136, 313)),
           w = 1.96,
           h = 2.36 * 2,
           xlabs = c(expression(eta[p]^2), expression("r"^2)),
           xl = c(0, 1),
           max_idx = c(20, 20),
           leg_id = 2,
           leg_locs = c(0.4, 20),
           figlabel = "B",
           figlabelon = TRUE)

p <- list(datpath = datpath,
          task = task,
          jmax = 2,
          dv = "dens_p",
          sel_n = paste(c(25, 59, 136, 313)),
          w = 1.96,
          h = 2.36 * 2,
          xlabs = c("p", "p"),
          xl = c(-50, 2),
          max_idx = c(5, 20),
          leg_id = 1,
          leg_locs = c(-45, 0.25),
          figlabel = "A",
          figlabelon = TRUE)

# ----------------------------------------------------
# KL divergence
# ----------------------------------------------------

kl <- list(datpath = datpath,
           task = task,
           dv = "dens_fx",
           ratio_type = "KL",
           origin = "313",
           sub_Ns = sub_Ns,
           w = 1.96,
           h = 2.36,
           leg_id = FALSE,
           leg_locs = c(5, 20),
           leg_txt = mods,
           y_label = expression(italic("KL p||q")),
           yl = NULL)

# ----------------------------------------------------
# fx sz ratio between models
# ----------------------------------------------------
           
model_rats <- list(datpath = datpath,
              task = task,
              dv = "stats_fx",
              ratio_type = "model",
              origin = "",
              sub_Ns = sub_Ns,
              w = 1.96,
              h = 2.36,
              leg_id = TRUE,
              leg_locs = c(5, 20),
              leg_txt = "",
              ylabel = expression(italic(paste(mods[1], "/", mods[2], sep=" "))),
              mods = mods,
              yl = NULL)

# ----------------------------------------------------
# meta-analytic vs observed mus
# ----------------------------------------------------
meta_mu <- list(datpath = datpath,
                task = task,
                mods = mods,
                sub_Ns = paste(round(exp(seq(log(13), log(313), length.out = 20)))),
                yl = c(-0.5, .5),
                leg_locs = c(2, -.05),
                leg_id = TRUE,
                sig_lines = NULL,
                sig_y = NULL)

# ----------------------------------------------------
# model mu difference
# ----------------------------------------------------
model_mu_diff <- list(datpath = datpath,
                      task = task,
                      mods = mods,
                      sub_Ns = paste(round(exp(seq(log(13), log(313), length.out = 20)))),
                      yl = c(-0.3, .9),
                      leg_locs = NULL,
                      leg_id = FALSE,
                      sig_lines = NULL,
                      sig_y = NULL)

# ----------------------------------------------------
# p ratios
# ----------------------------------------------------
p_rat <- list(datpath = datpath,
              task = task,
              dv = "stats_p",
              ratio_type = "origin",
              origin = medN,
              y_label = expression(italic("q ratio")),
              sub_Ns = sub_Ns,
              leg_id = TRUE,
              leg_locs = c(1, 0.75),
              leg_txt = mods,
              yl = c(0, 2))

# ----------------------------------------------------
# meta-analytic vs observed fx sz ratio
# ----------------------------------------------------
sig <- list(datpath = datpath,
            task = task,
            dv = "stats_sig",
            ratio_type = "stats_sig",
            origin = "",
            sub_Ns = sub_Ns,
            w = 1.96,
            h = 2.36,
            leg_id = TRUE,
            leg_locs = c(5, 0.5),
            leg_txt = mods,
            ylabel = expression(italic("meta / sim")),
            mods = mods,
            yl = c(0, 2))

# ----------------------------------------------------
# qqplot settings
# ----------------------------------------------------
qq_inputs <- list(datpath = datpath,
                  task = task,
                  sub_Ns = sub_Ns,
                  median_N = as.numeric(medN),
                  leg_id = TRUE,
                  leg_locs = c(0.01, .9),
                  leg_txt = mods,
                  xl = c(0, 1),
                  yl = c(0, 1))

save(task, fname, fx, p, kl, meta_mu, model_mu_diff, model_rats, sig, p_rat,
     qq_inputs,
     file = paste("../data/", task, "/",
                  task, "_plot_settings.RData", sep = ""))


