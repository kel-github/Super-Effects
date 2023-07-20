##############################################################################
## functions
############################################################################

###########################################################################
# BASICS TO GET DATA AND PERFORM STEPWISE - TASK GENERIC
##########################################################################
get_dat_4corrs <- function(datpath, task, zipnm, N){
  # for a given N, compile fx size and distributional data
  # datpath [str] - where is the data - e.g. "../data"
  # task [str] - which task are you getting data for?
  # zipnm [str] - what is the name of the zipped folder?
  # N [int] - which sample size are we working on?
  
  # first, list the files in the zip rxv
  tmp = unzip(zipfile=paste(datpath, task, zipnm, sep="/"), list = TRUE)
  tmp = tmp$Name[!is.na(str_extract(tmp$Name, paste("N-", N, "_", sep = "")))]
  dst_idx <- do.call(cbind, lapply(tmp, str_detect, pattern = "diststats"))
  fnms <- tmp[!dst_idx]
  # now get the fx data
  fx_dat <- do.call(rbind, lapply(fnms, unzp_fx_n_spew, datpath=datpath, task=task, zipnm=zipnm))
  # and reshape/relabel it if it has 2 effects in it
  if(task == "SD"){
    fx_dat$mod[fx_dat$mod == "RM-AN"] = "ME"
    fx_dat$mod[fx_dat$mod == "LME"] = "int"
    fx_dat <- fx_dat %>% pivot_wider(id_expand = FALSE, names_from = mod, values_from = c(esz, p))
  }
  # and the mu diff and var stats
  dst_nms <- tmp[dst_idx]
  dst_dat <- do.call(rbind, lapply(dst_nms, unzp_mu_stats_n_spew, 
                                   datpath=datpath,
                                   task=task,
                                   zipnm=zipnm))
  
  # join up fx_dat and dst_dat and return
  inner_join(fx_dat, dst_dat, by = "parent")
}


unzp_fx_n_spew <- function(datpath, task, zipnm, fnm){
  # given a list of filenames, unextract
  # each one, open, and add to dataframe
  # output the key variables from the 
  # dataframe
  f <- unzip(paste(datpath, task, zipnm, sep="/"), 
             file=fnm,
             exdir=paste(datpath, task, sep="/"))
  load(f)
  row_idx <- !is.na(out$p)
  out <- out[row_idx, ]
  # add parent #
  pN <- strsplit(fnm, split = "-")
  pNidx <- !is.na(sapply(pN, str_extract, pattern = ".RData"))
  nustr <- pN[[1]][pNidx]
  parent <- as.numeric(strsplit(nustr, split="R")[[1]][1])
  out$parent <- parent
  
  # will need to pivot wider here for tasks with 2 fx
  
  # output
  out %>% select(n, p, esz, mod, parent)
}

unzp_mu_stats_n_spew <- function(datpath, task, zipnm, fnm){
  # given a list of filenames, unextract
  # each one, open, and add to dataframe
  # output the key variables from the 
  # dataframe
  f <- unzip(paste(datpath, task, zipnm, sep="/"), 
             file=fnm,
             exdir=paste(datpath, task, sep="/"))
  load(f)
  
  if (task == "AB" | task == "SRT") {
    dst_dat <- cbind(dist_out[[1]][1], dist_out[[2]])
  } else if (task == "SD"){
    dst_dat <- cbind(dist_out[[1]]$ME, dist_out[[1]]$int, dist_out[[2]])
    names(dst_dat) <- c("ME_mu", "ME_sigma", "ME_skew", "ME_k", "ME_r",
                        "int_mu", "int_sigma", "int_skew", "int_k", "int_r",
                        "sigma_mu", "skew_mu", "k_mu")
  } 
  # now get parent number
  
  # add parent #
  pN <- strsplit(fnm, split = "-")
  pNidx <- !is.na(sapply(pN, str_extract, pattern = "_diststats"))
  nustr <- pN[[1]][pNidx]
  parent <- as.numeric(strsplit(nustr, split="_")[[1]][1])
  dst_dat$parent <- parent
  
  # output
  dst_dat
}

get_ev_non_norm <- function(x){
  # given a random variable x, compute the expected value
  # note, check with collaborators they're happy with this
  d <- density(x)
  mu <- sum(d$x*d$y)*diff(d$x)[1]
  mu
}

do_stp_n_prs <- function(df, dv, task, datpath = "../data/", prs_fnm, rsd_fm){
  # given the dataframe df, perform an AIC stepwise
  # regression process, with dv as the dv
  # dv [str], corresponding to a named variable in df
  # the remaining columns in df will serve as the predictor 
  # variables
  # task [str] - which task are you performing analysis on?
  # datpath [str] - where is the parent data folder?
  # prs_fnm [str] - what base do you want the plot fnms to have?
  # rsd_fnm [str] - same as above but for residuals
  # first do pairs plot and save
  prs <- ggpairs(df)
  ggsave(paste(prs_fnm, "pairs.pdf", sep="_"), plot=prs, 
         path=paste(datpath, task, sep="/"),
         width=30,
         height=30,
         units = "cm")
  
  # perform linear regression and test the VIF factor
  full_model <- lm(as.formula(paste(dv, "~ .")), data=df)
  vif_out <- vif(full_model)
  step_model <- stepAIC(full_model, direction = "both", 
                        trace = FALSE)
  pdf(paste(datpath, task, paste(rsd_fm, "resid.pdf", sep="_"), sep="/"))
  plot(step_model, which = c(1, 2))
  dev.off()
  
  pdf(paste(datpath, task, paste(rsd_fm, "aV.pdf", sep="_"), sep="/"))
  avPlots(step_model)
  dev.off()
  
  list(vif_out, step_model)
}

r2z <- function(r){
  # Fisher R to Z transform https://en.wikipedia.org/wiki/Fisher_transformation
  0.5*log((1+r)/(1-r))
}

###########################################################################
# TASK SPECIFIC - COLLATE MODEL OVER N
##########################################################################
run_models <- function(f, Ns, datpath, task, zipnm, maxN){
  
  betas <- lapply(Ns, eval(f), datpath=datpath,
                  task=task,
                  zipnm=zipnm,
                  maxN=maxN)
  names(betas) <- as.character(Ns)
  save(betas, file=paste(datpath, task, paste(task, "moments_betas.RData", 
                                              sep="_"), sep="/"))
}

run_AB_models <- function(datpath = "../data", task = "AB", zipnm = "EPSAB_wv.zip", N, maxN = 313){
  # for one N, run the model on the AB data
  # and save the results
  # datpath [str] - where is the data - e.g. "../data"
  # task [str] - which task are you getting data for?
  # zipnm [str] - what is the name of the zipped folder?
  # N [int] - which sample size are we working on?
  dat <- get_dat_4corrs(datpath,"AB", zipnm, N)
  Nmaxdat <- get_dat_4corrs(datpath,"AB", zipnm, maxN)
  
  mu <- get_ev_non_norm(Nmaxdat$esz)
  # now compute the distance of each value of the fx size
  dat$esz_dist <- dat$esz - mu
  
  AB_model <- lm(esz_dist ~ AB.skew + AB.k + AB.r + sigma_mu, data = dat)
  
  prs <- ggpairs(dat %>% select(esz_dist, AB.skew, AB.k, AB.r, sigma_mu))
  ggsave(paste(task, N, "pairs.pdf", sep="_"), plot=prs, 
         path=paste(datpath, task, sep="/"),
         width=30,
         height=30,
         units = "cm") 
  
  pdf(paste(datpath, task, paste(task, N, "resid.pdf", sep="_"), sep="/"))
  plot(AB_model, which = c(1, 2))
  dev.off()
  
  pdf(paste(datpath, task, paste(task, N, "aV.pdf", sep="_"), sep="/"))
  avPlots(AB_model)
  dev.off()
  
  # return coefficients
  AB_model$coefficients
}

run_SRT_models <- function(datpath = "../data", task = "SRT", 
                           zipnm = "SRT_wv.zip", N, maxN = 313){
  
  # for one N, run the model and save the results
  # datpath [str] - where is the data - e.g. "../data"
  # task [str] - which task are you getting data for?
  # zipnm [str] - what is the name of the zipped folder?
  # N [int] - which sample size are we working on?
  dat <- get_dat_4corrs(datpath, task, zipnm, N)
  Nmaxdat <- get_dat_4corrs(datpath, task, zipnm, maxN)
  mu <- get_ev_non_norm(Nmaxdat$esz)
  # now compute the distance of each value of the fx size
  dat$esz_dist <- dat$esz - mu
  
  # now transform the data
  x = 0.01 - min(dat$esz_dist)
  dat <- dat %>% mutate(esz_dist_t = log(esz_dist + x),
                          SRT.r_t = scale(r2z(SRT.r), scale=F),
                          SRT.k_t = log(SRT.k), # same
                          sigma_mu_t = log(sigma_mu)) %>%
                  select(esz_dist_t, SRT.skew, SRT.k_t, SRT.r_t, sigma_mu)
  SRT_model <- lm(esz_dist_t ~ SRT.skew + SRT.k_t + SRT.r_t + sigma_mu, data = dat)
  
  prs <- ggpairs(dat)
  ggsave(paste(task, N, "pairs.pdf", sep="_"), plot=prs, 
         path=paste(datpath, task, sep="/"),
         width=30,
         height=30,
         units = "cm") 
  
  pdf(paste(datpath, task, paste(task, N, "resid.pdf", sep="_"), sep="/"))
    plot(SRT_model, which = c(1, 2))
  dev.off()
  
  pdf(paste(datpath, task, paste(task, N, "aV.pdf", sep="_"), sep="/"))
    avPlots(SRT_model)
  dev.off()
  
  # return coefficients
  SRT_model$coefficients
}

run_SD_ME <- function(datpath = "../data", task = "SD", 
                      zipnm = "SD_wv.zip", N, maxN = 313){
  
  dat <- get_dat_4corrs(datpath, task, zipnm, N)
  Nmaxdat <- get_dat_4corrs(datpath, task, zipnm, maxN)
  mu <- get_ev_non_norm(Nmaxdat$esz_ME)
  # now compute the distance of each value of the fx size
  dat$esz_dist <- dat$esz_ME - mu
  
  # now transform
  dat <- dat %>% select(esz_dist, ME_skew, ME_k,
                          ME_r, sigma_mu) %>%
                 mutate(ME_k_t = log(ME_k),
                          ME_r_t = r2z(ME_r)) %>%
                 select(esz_dist, ME_skew, ME_k_t,
                          ME_r_t, sigma_mu)

  SD_ME_model <- lm(esz_dist ~ ME_skew + ME_k_t + ME_r_t + sigma_mu, data = dat)
  
  prs <- ggpairs(dat)
  ggsave(paste(task, "ME", N, "pairs.pdf", sep="_"), plot=prs, 
         path=paste(datpath, task, sep="/"),
         width=30,
         height=30,
         units = "cm") 
  
  pdf(paste(datpath, task, paste(task, "ME", N, "resid.pdf", sep="_"), sep="/"))
    plot(SD_ME_model, which = c(1, 2))
  dev.off()
  
  pdf(paste(datpath, task, paste(task, "ME", N, "aV.pdf", sep="_"), sep="/"))
    avPlots(SD_ME_model)
  dev.off()
  
  # return coefficients
  SD_ME_model$coefficients
}

run_SD_int <- function(datpath = "../data", task = "SD", 
                      zipnm = "SD_wv.zip", N, maxN = 313){
  
    dat <- get_dat_4corrs(datpath, task, zipnm, N)
    Nmaxdat <- get_dat_4corrs(datpath, task, zipnm, maxN)
    mu <- get_ev_non_norm(Nmaxdat$esz_int)
    # now compute the distance of each value of the fx size
    dat$esz_dist <- dat$esz_int - mu
   
    # now transform
    x <- 0.01 - min(dat$esz_dist)
    dat <- dat %>% mutate(esz_dist_t = log(esz_dist + x),
                                           int_k_t = log(int_k),
                                           k_mu_t = log(k_mu)) %>%
                  select(esz_dist_t, int_skew, int_k_t,
                          int_r, sigma_mu)
    
    # run model
    SD_int_model <- lm(esz_dist_t ~ int_skew + int_k_t + int_r + sigma_mu, data = dat)
    
    prs <- ggpairs(dat)
    ggsave(paste(task, "int", N, "pairs.pdf", sep="_"), plot=prs, 
           path=paste(datpath, task, sep="/"),
           width=30,
           height=30,
           units = "cm") 
    
    pdf(paste(datpath, task, paste(task, "int", N, "resid.pdf", sep="_"), sep="/"))
      plot(SD_int_model, which = c(1, 2))
    dev.off()
    
    pdf(paste(datpath, task, paste(task, "int", N, "aV.pdf", sep="_"), sep="/"))
      avPlots(SD_int_model)
    dev.off()

    # return coefficients
    SD_int_model$coefficients    
}    