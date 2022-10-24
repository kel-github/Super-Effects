
do.AB.analysis <- function(fname){
  # ----------------------------------------------------------------------------------------------------
  # AB specific analysis
  # ---------------------------------------------------------------------------------------------------- 
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
  eps <- sapply(unique(an$Effect), function(x) (an$F[an$Effect == x] - 1) / 
                                      (an$F[an$Effect == x] + (an$DFd[an$Effect == x]/an$DFn[an$Effect == x])))
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
  list(an, ts, eps)
  
}

do.MT.analysis <- function(fname){
  # ----------------------------------------------------------------------------------------------------
  # MT specific analysis
  # ---------------------------------------------------------------------------------------------------- 
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

  options(contrasts = c("contr.sum", "contr.poly")) # set options  
  
  # RUN RM ANOVA
  # -----------------------------------------------------------------------------
  an <- rstatix::get_anova_table(rstatix::anova_test(data=ffx.dat %>% ungroup(), dv=RT, wid=sub, within=c(task,trialtype), effect.size="pes", type=3))

  eps <- sapply(unique(an$Effect), function(x) (an$F[an$Effect == x] - 1) / 
                  (an$F[an$Effect == x] + (an$DFd[an$Effect == x]/an$DFn[an$Effect == x])))
  # ----------------------------------------------------------------------------------------------------
  # do post-hoc t - tests
  # ----------------------------------------------------------------------------------------------------
  post_hoc_dat <- ffx.dat %>% pivot_wider(id_cols=sub, names_from = c(task, trialtype), values_from = RT) %>% 
                              mutate(sound_MT = sound_dual - sound_single) %>% 
                              mutate(vis_MT = vis_dual - vis_single) %>%
                              mutate(mu_sing = mean(sound_single, vis_single)) %>%
                              mutate(mu_dual = mean(sound_dual, vis_dual))
  ts = list(with(post_hoc_dat, t.test(sound_MT,
                                      vis_MT,
                                      paired=TRUE,
                                      var.equal=FALSE)),
            with(post_hoc_dat, t.test(mu_sing,
                                      mu_dual, 
                                      paired = TRUE,
                                      var.equal = FALSE)))
  names(ts) <- c("interaction", "me_MT")
  
  # ----------------------------------------------------------------------------------------------------
  # return
  # ----------------------------------------------------------------------------------------------------
  list(an, ts, eps)
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
  eps <- sapply(unique(an$Effect), function(x) (an$F[an$Effect == x] - 1) / 
                  (an$F[an$Effect == x] + (an$DFd[an$Effect == x]/an$DFn[an$Effect == x])))
  
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
  list(an, ts, eps)
  
}

do.SRT.analysis <- function(fname){
  # ----------------------------------------------------------------------------------------------------
  # SRT specific analysis
  # ----------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------
  # load data and wrangle into tidy form
  # ----------------------------------------------------------------------------------------------------

  dat <- read.csv(fname, header=TRUE)
  # ----------------------------------------------------------------------------------------------------
  # Create dataframes
  # ----------------------------------------------------------------------------------------------------
  # Create a summary of the data for fixed fx modelling
  min.RT <- 200 # in msec
  sd.crit <- 2.5

  ffx.dat <- dat %>% filter(Block.No > 2) %>%
    group_by(Subj.No, Block.No.Names) %>%
    filter(Accuracy == 1) %>%
    filter(RT.ms > min.RT) %>%
    filter(RT.ms < (mean(RT.ms) + sd.crit*sd(RT.ms))) %>%
    summarise(RT=mean(RT.ms))
  names(ffx.dat)[names(ffx.dat) == "Block.No.Names"] <- "trialtype"

  t <- with(ffx.dat, t.test(x = RT[trialtype == "Random Block"],
                            y = RT[trialtype == "Sequence Block"],
                            paired = TRUE,
                            var.equal = FALSE))
  d <- get.cohens.d(ffx.dat, iv = "trialtype", dv = "RT", x = "Random Block")
  r <- (d / (sqrt((d^2) + 4)))^2
  list(t, d, r)
}

d2r <- function(dat, m) {
  # apply d2r transform for data in dat corresponding to model m
  dat %>% filter(mod == m) %>%
    group_by(n) %>%
    mutate(esz = (esz / (sqrt((esz^2) + 4)))^2) %>% # Cohen 1988 equation 2.2.6
    ungroup()
}


