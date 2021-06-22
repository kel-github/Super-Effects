# ----------------------------------------------------------------------------------------------------
# AB specific analysis
# ---------------------------------------------------------------------------------------------------- 
plot.AB.results <- function(fname){
  ## plot the mean T1 and T2|T1 group level data
  
  # ----------------------------------------------------------------------------------------------------
  # read in data
  # ----------------------------------------------------------------------------------------------------  
  dat <- read.csv(fname, header=TRUE)
  
  # ----------------------------------------------------------------------------------------------------
  # create summary dataframe
  # ----------------------------------------------------------------------------------------------------
  
  # Create a summary of the data for fixed fx modelling
  ffx.dat <- dat %>% group_by(Subj.No, Trial.Type.Name) %>%
                              summarise(T1.Accuracy=mean(T1.Accuracy),
                                        T2T1.Accuracy=mean(T2T1.Accuracy)) %>%
                     pivot_longer(cols=c("T1.Accuracy", "T2T1.Accuracy"),
                                  names_to = "condition",
                                  values_to = "accuracy") %>%
                     group_by(Trial.Type.Name, condition) %>%
                              summarise(acc=mean(accuracy),
                                        se=sd(accuracy)/sqrt(length(accuracy))) %>% 
                     ungroup()
  names(ffx.dat)[names(ffx.dat)=="Trial.Type.Name"] = "lag"
  ffx.dat$lag <- fct_recode(ffx.dat$lag, "2" = "lag_2", "3" = "lag_3", "5" = "lag_5", "7" = "lag_7")
  ffx.dat$condition <- fct_recode(ffx.dat$condition, "T1" = "T1.Accuracy", "T2gT1" = "T2T1.Accuracy")
  # ----------------------------------------------------------------------------------------------------
  # plot data
  # ----------------------------------------------------------------------------------------------------
  ffx.dat %>% ggplot(aes(x=lag, y=acc, group=condition, fill=condition)) +
    geom_line(aes(colour=condition)) +
    geom_errorbar(aes(ymin=acc-(1.96*se), ymax=acc+(1.96*se), colour=condition), width=0.2) +
    scale_colour_manual(values=wes_palette("IsleofDogs1")[c(3,2)]) +
    theme_cowplot() + ylim(c(0,1)) + ylab("accuracy")
  
}

do.AB.analysis <- function(fname){
  
  # ----------------------------------------------------------------------------------------------------
  # load data and wrangle into tidy form 
  # ----------------------------------------------------------------------------------------------------
  
  dat <- read.csv(fname, header=TRUE)
  # ----------------------------------------------------------------------------------------------------
  # Create dataframes 
  # ----------------------------------------------------------------------------------------------------
  
  ffx.dat <- dat %>% group_by(Subj.No, Trial.Type.Name) %>%
    summarise(T1=mean(T1.Accuracy),
              T2gT1=mean(T2T1.Accuracy)) %>%
              ungroup()
  names(ffx.dat)[c(1:2)] <- c("sub", "lag")
  # ----------------------------------------------------------------------------------------------------
  # run anova
  # ----------------------------------------------------------------------------------------------------
  an <- get_anova_table(anova_test(data=ffx.dat, dv=T2gT1, wid=sub, within=lag, effect.size="pes", type=3))
  
  # ----------------------------------------------------------------------------------------------------
  # do post-hoc t - tests
  # ----------------------------------------------------------------------------------------------------
  a = c("lag_2", "lag_2", "lag_2", "lag_3", "lag_3", "lag_5")
  b = c("lag_3", "lag_5", "lag_7", "lag_5", "lag_7", "lag_7")
  ts = mapply(function(x,y) t.test(ffx.dat$T2gT1[ffx.dat$lag == x],
                                   ffx.dat$T2gT1[ffx.dat$lag == y],
                                   paired=TRUE,
                                   var.equal=FALSE), a, b)
 
  # ----------------------------------------------------------------------------------------------------
  # return
  # ----------------------------------------------------------------------------------------------------
  list(an, ts)
  
}


# ----------------------------------------------------------------------------------------------------
# CC specific analysis
# ---------------------------------------------------------------------------------------------------- 
plot.CC.results <- function(fname){
  ## plot the trial type x block interaction for mean RTs

  # ----------------------------------------------------------------------------------------------------
  # read in data
  # ----------------------------------------------------------------------------------------------------  
  dat <- read.csv(fname, header=TRUE)
  
  # ----------------------------------------------------------------------------------------------------
  # create summary dataframe
  # ----------------------------------------------------------------------------------------------------
  min.RT <- 200 # in msec
  sd.crit <- 2.5
  
  ffx.dat <- dat %>% mutate(Block.No = rep(c(1:12), each = 24, length(unique(dat$Subj.No)))) %>%
    group_by(Subj.No, Block.No, Trial.Type.Name) %>%
    filter(Accuracy == 1) %>%
    filter(RT.ms > min.RT) %>%
    filter(RT.ms < (mean(RT.ms) + sd.crit*sd(RT.ms))) %>%
    summarise(RT=mean(RT.ms))
  subs  <- unique(ffx.dat$Subj.No)
  names(ffx.dat) <- c("sub", "block", "type", "RT")
  
  # Create a summary of the data for fixed fx depiction
  ffx.dat <- ffx.dat %>% group_by(sub, block, type) %>%
    summarise(RT=mean(RT)) %>%
    group_by(block, type) %>%
    summarise(mu=mean(RT),
              se=sd(RT)/sqrt(length(RT))) %>% 
    ungroup()
  
  # ----------------------------------------------------------------------------------------------------
  # plot data
  # ----------------------------------------------------------------------------------------------------
  ffx.dat %>% ggplot(aes(x=as.factor(block), y=mu, group=type, fill=type)) +
    geom_line(aes(colour=type)) +
    geom_errorbar(aes(ymin=mu-(1.96*se), ymax=mu+(1.96*se), colour=type), width=0.2) +
    scale_colour_manual(values=wes_palette("IsleofDogs1")[c(3,2)]) +
    theme_cowplot() + ylab(paste(bquote(mu), "RT", sep=" ")) + xlab("block")
  
}

do.CC.analysis <- function(fname){
  
  # ----------------------------------------------------------------------------------------------------
  # load data and wrangle into tidy form 
  # ----------------------------------------------------------------------------------------------------
  
  dat <- read.csv(fname, header=TRUE)
  # ----------------------------------------------------------------------------------------------------
  # create summary dataframe for CC data
  # ----------------------------------------------------------------------------------------------------
  min.RT <- 200 # in msec
  sd.crit <- 2.5
  ffx.dat <- dat %>% mutate(Block.No = rep(c(1:12), each = 24, length(unique(dat$Subj.No)))) %>%
    group_by(Subj.No, Block.No, Trial.Type.Name) %>%
    filter(Accuracy == 1) %>%
    filter(RT.ms > min.RT) %>%
    filter(RT.ms < (mean(RT.ms) + sd.crit*sd(RT.ms))) %>%
    summarise(RT=mean(RT.ms))
  subs  <- unique(ffx.dat$Subj.No)
  names(ffx.dat) <- c("sub", "block", "type", "RT")

  # ----------------------------------------------------------------------------------------------------
  # run anova
  # ----------------------------------------------------------------------------------------------------
  an <- rstatix::get_anova_table(rstatix::anova_test(data=ffx.dat%>%ungroup(), dv=RT, wid=sub, within=c(block,type), effect.size="pes", type=3))
  
  # ----------------------------------------------------------------------------------------------------
  # do post-hoc t - tests on 1st block and 12th block
  # ----------------------------------------------------------------------------------------------------
  a = c(1, 12)
  ts = lapply(a, function(x) t.test(ffx.dat$RT[ffx.dat$block == x & ffx.dat$type == "Novel"],
                                    ffx.dat$RT[ffx.dat$block == x & ffx.dat$type == "Repeated"],
                                    paired=TRUE,
                                    var.equal=FALSE))
  
  # ----------------------------------------------------------------------------------------------------
  # return
  # ----------------------------------------------------------------------------------------------------
  list(an, ts)
  
}

CC.sum.dat <- function(data){
  # ----------------------------------------------------------------------------------------------------
  # create summary dataframe for CC data
  # ----------------------------------------------------------------------------------------------------
  min.RT <- 200 # in msec
  sd.crit <- 2.5
  
  ffx.dat <- dat %>% mutate(Block.No = rep(c(1:12), each = 24, length(unique(dat$Subj.No)))) %>%
    group_by(Subj.No, Block.No, Trial.Type.Name) %>%
    filter(Accuracy == 1) %>%
    filter(RT.ms > min.RT) %>%
    filter(RT.ms < (mean(RT.ms) + sd.crit*sd(RT.ms))) %>%
    summarise(RT=mean(RT.ms))
  subs  <- unique(ffx.dat$Subj.No)
  names(ffx.dat) <- c("sub", "block", "type", "RT")
  ffx.dat
}


plot_fx_szs <- function(d, med, d_scale_ffx, d_scale_rfx, px_rng_d){
  ### Call plotting functions for observed effect sizes 
  # kwargs
  # -- d: dataframe
  # -- med: median sample size
  # -- d_scale_ffx: scaling factor (ffx)
  # -- d_scale_rfx: scaling factor (lme/prev)
  # -- px_rng_d: range of values to plot
  
  # ----------------------------------------------------------------------------------------------------
  # define factors and plot
  # ----------------------------------------------------------------------------------------------------
  d$Nsz <- as.factor(d$Nsz)
  d$mod <- as.factor(d$mod)
  plts <- lapply(print(levels(d$mod)), function(x) plot.d(d, x, px_rng_d, d_scale_ffx, med))
  # ----------------------------------------------------------------------------------------------------
  # return!
  # ----------------------------------------------------------------------------------------------------  
  plts
}
  
plot_ps <- function(d){  
  # ----------------------------------------------------------------------------------------------------
  # load p, define factors and plot
  # ----------------------------------------------------------------------------------------------------
  
  load(fnames[2])
  
  
  d$Nsz <- as.factor(d$Nsz)
  d$mod <- as.factor(d$mod)
  
  if (task == "SD") d <- d %>% filter(as.numeric(Nsz) < 17) # remove sample sizes saturated at 0
  
  ffx.p <- plot.p(d, "ffx", px_rng_p_ffx, p_scale_ffx, med)
  rfx.p <- plot.p(d, "rfx", px_rng_p_rfx, p_scale_rfx, med)
  
  list(ffx.p, rfx.p)
  # p = plot_grid(ffx.d, rfx.d, ffx.p, rfx.p, labels=c('A', 'B', 'C', 'D'), label_size = 12, align="v")
  # p
  # #p # print out the plot so you can see it
  #p = p + ggsave(paste("../images/", task, ".png", sep=""), width = width, height = height, units="in")
}


do.stats <- function(subfol, task){
  
  # ----------------------------------------------------------------------------------------------------
  # load d values
  # ----------------------------------------------------------------------------------------------------
  
  fnames = c(paste("../data/", subfol, "/", task, "_esz", "_d.RData", sep=""), paste("../data/", subfol, "/", task, "_p", "_d.RData", sep=""))
  load(fnames[1])
  
  # ----------------------------------------------------------------------------------------------------
  # define factors 
  # ----------------------------------------------------------------------------------------------------
  d$Nsz <- as.factor(d$Nsz)
  d$mod <- as.factor(d$mod)
  
  total_p <- d %>% group_by(mod, Nsz) %>% summarise(sum=sum(d)) 
  d <- d %>% inner_join(total_p, by=c("mod", "Nsz")) %>% mutate(dp = d/sum) 
  
  # ----------------------------------------------------------------------------------------------------
  # get ratios
  # ----------------------------------------------------------------------------------------------------
  d <- d %>% group_by(mod, Nsz) %>% mutate(c=cumsum(dp))
  
  get.ratios <- function(dist_a, dist_b){
    # get the ratio of the 95% iqr of two distributions
    dist_a_IQR <- abs(min(dist_a$x[which(abs(dist_a$c-.025) == min(abs(dist_a$c-.025)))]) - 
                      min(dist_a$x[which(abs(dist_a$c-.975) == min(abs(dist_a$c-.975)))]))
    
    dist_b_IQR <- abs(min(dist_b$x[which(abs(dist_b$c-.025) == min(abs(dist_b$c-.025)))]) - 
                        min(dist_b$x[which(abs(dist_b$c-.975) == min(abs(dist_b$c-.975)))]))
    
    dist_a_IQR/dist_b_IQR
  }
  
  ffx_ratios <- lapply(levels(d$Nsz)[1:(length(levels(d$Nsz))-1)], 
                       function(x) get.ratios(dist_a = d[d$mod == "ffx" & d$Nsz == x,], dist_b = d[d$mod == "ffx" & d$Nsz == "313",]))
  
  rfx_ratios <- lapply(levels(d$Nsz)[1:length(levels(d$Nsz))], 
                       function(x) get.ratios(dist_a = d[d$mod == "ffx" & d$Nsz == x,], dist_b = d[d$mod == "rfx" & d$Nsz == x,]))
  

  # ----------------------------------------------------------------------------------------------------
  # load p values
  # ----------------------------------------------------------------------------------------------------
  

  load(fnames[2])
  # ----------------------------------------------------------------------------------------------------
  # define factors 
  # ----------------------------------------------------------------------------------------------------
  d$Nsz <- as.factor(d$Nsz)
  d$mod <- as.factor(d$mod)
  
  total_p <- d %>% group_by(mod, Nsz) %>% summarise(sum=sum(d)) 
  d <- d %>% inner_join(total_p, by=c("mod", "Nsz")) %>% mutate(dp = d/sum)   
  
  get.propor <- function(x, frq){
    sum(frq[x<log(.05)])/sum(frq)
  }
  
  pwr.ffx <- do.call(cbind, lapply(levels(d$Nsz), function(y) get.propor(x=d$x[d$Nsz == y & d$mod == "ffx"], frq=d$d[d$Nsz == y & d$mod == "ffx" ])))
  pwr.rfx <- do.call(cbind, lapply(levels(d$Nsz), function(y) get.propor(x=d$x[d$Nsz == y & d$mod == "rfx"], frq=d$d[d$Nsz == y & d$mod == "rfx" ])))
  
  list(ffx_ratios = unlist(ffx_ratios), ffx_vs_rfx = unlist(rfx_ratios), pwr_ffx = pwr.ffx, pwr_rfx = pwr.rfx)
}

