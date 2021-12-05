rm(list=ls())

task <- "CC"
datpath <- "../data/"
medN <- "25"
mods <- c("RM-AN", "LME") #c("t", "p")
imm = TRUE
# ----------------------------------------------------
# behavioural data
# ----------------------------------------------------
fname <- "../data/total_of_313_subs_CC_task_trial_level_data.csv"
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
           xlabs = expression(eta[p]^2, eta[p]^2), #expression("d"), 
           xl =  c(0.0, 0.4, 0.0, 1),
           max_idx = c(20, 20),
           leg_id = 1,
           leg_locs = c(0.25, 49.0),
           figlabel = "A",
           figlabelon = TRUE,
           imm = imm)

p <- list(datpath = datpath,
          task = task,
          jmax = 2,
          dv = "dens_p",
          sel_n = paste(c(25, 59, 136, 313)),
          w = 1.96,
          h = 2.36 * 2,
          xlabs = c("p", "p"),
          xl = c(-40, 2, -40, 2),
          max_idx = c(1, 1),
          leg_id = 1,
          leg_locs = c(-40, 0.35),
          leg_txt = c("b x tt", "tt"),
          figlabel = "A",
          figlabelon = TRUE,
          imm = imm)

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
           leg_txt = c("b x tt", "tt"),
           y_label = expression(italic("KL p||q")),
           yl = NULL,
           xvals = c(0, 1, 0, 1),
           mods = c("RM-AN", "LME"),
           imm = imm)

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
                yl = c(-0.25, .15),
                leg_locs = c(2, 0.145),
                leg_id = TRUE,
                leg_txt = c("b x tt", "tt"),
                sig_lines = NULL,
                sig_y = NULL,
                imm = imm)

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
              leg_locs = c(1, 0.95),
              leg_txt = mods,
              yl = c(0, 2.5),
              imm = imm)

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
                  leg_locs = c(0.01, .79),
                  leg_txt = c("b x tt", "tt"),
                  xl = c(0.0, 0.8),
                  yl = c(0.0, 0.8),
                  mods = c("RM-AN", "LME"),
                  modn = 2,
                  imm = imm)

if (imm){
  save(task, fname, fx, p, kl, meta_mu, model_mu_diff, model_rats, sig, p_rat,
       qq_inputs,
       file = paste("../data/", task, "/", "IMM",
                    task, "_plot_settings.RData", sep = ""))
} else {
  save(task, fname, fx, p, kl, meta_mu, model_mu_diff, model_rats, sig, p_rat,
       qq_inputs,
       file = paste("../data/", task, "/",
                    task, "_plot_settings.RData", sep = ""))
}
