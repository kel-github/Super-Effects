### written by K. Garner, April 2020
### edited by Z. Nott, September 2020
### much enhanced by C. Nolan, Jan+ 2021
### for the project 'On the detectability of effects in executive
### function and implicit learning tasks'

### custom functions for running statistical analyses and plotting 

# ------------------------------------------------------------------------------
# select a data subset
# ------------------------------------------------------------------------------
sample.N <- function(subs, N, k, replace) {
  # get k x N subject numbers and assign to a dataframe for future filtering  
  this.sample <- data.frame(sub=unlist(lapply(N, sample, x=subs, replace=replace)),
                            Nsz=unlist(lapply(N,  function(x) rep(x, each=x))),
                            perm=k)
  this.sample
}

# ------------------------------------------------------------------------------
###### t-test functions (paired samples)
###### -------------------------------------------------------------------------

get.ps.t.test <- function(data, iv, dv, x){
  # run t.test on the dv between the variables x & y 
  # return the p value
  # data = dataframe for testing
  # iv = name of iv
  # dv = name of dv
  # x = first level of iv
  # y = second level of iv
  t.dat = data[eval(dv)]
  t.idx = data[eval(iv)] == x
  t = t.test(t.dat[t.idx == TRUE], t.dat[t.idx==FALSE], paired = TRUE)
  t$p.value
}

get.cohens.d <- function(data, iv, dv, x){
  # get Cohen's d measure for paired samples
  # uses Cohen's d_z (equation 6 in https://www.frontiersin.org/articles/10.3389/fpsyg.2013.00863/full)
  # data = dataframe for testing
  # iv = name of iv
  # dv = name of dv
  # x = first level of iv
  # y = second level of iv  
  d.dat = data[eval(dv)]
  d.idx = data[eval(iv)] == x
  mu_diff = (mean(d.dat[d.idx == TRUE]) - mean(d.dat[d.idx == FALSE]))
  sd.pool = sum(((d.dat[d.idx == TRUE] - d.dat[d.idx == FALSE]) - mu_diff)^2)
  sd.pool = sqrt(sd.pool / (length(d.dat[d.idx == TRUE])-1))
#  sd.pool = sqrt( ( sd( d.dat[d.idx == TRUE]  )^2 + sd(  d.dat[d.idx == FALSE]  )^2 ) / 2 )
  d =  mu_diff/ sd.pool
  d
}

run.t.test.sim <- function(data, iv="trialtype", dv="RT", x="Random Block", subs){
  
  out = data.frame( p    = get.ps.t.test(data, iv, dv, x),
                    d    = get.cohens.d(data, iv, dv, x))
  out
}

###### -------------------------------------------------------------------------
###### t-test functions (one sample t test, used for VSL)
###### -------------------------------------------------------------------------

get.os.t.test <- function(sum.data, dv){
  # run one sample t.test 
  # return the p value
  # data = dataframe for testing
  # iv = name of iv
  # dv = name of dv
  # x = iv
  t.dat = sum.data[eval(dv)]
  #  t.idx = data[eval(iv)] == x
  #  t = t.test(t.dat[t.idx == TRUE], alt = "greater", mu = 0.5)
  t = t.test(t.dat, alt = "greater", mu = 0.5)
  t
}

get.os.cohens.d <- function(sum.data, dv){
  # get Cohen's d measure for one sample t test
  # sum.data = dataframe for testing
  # dv = dv of interest
  
  d.dat = sum.data[eval(dv)]
  meanH0 = 0.5
  sd = sd( d.dat$acc )
  d = (mean(d.dat$acc) - meanH0) / sd
  d
}

run.os.t.test.sim <- function(data, dv = "acc"){
  
  # first label subs
  u <- length(unique(data$sub))
  data$sub <- rep(1:(length(data$sub)/max(data$Trial.No)), each = max(data$Trial.No))
  # get acc per sub
  sum.data = data %>% group_by(sub) %>%
    summarise(acc=mean(Response==Target.Order))
  # do tests
  t = get.os.t.test(sum.data, dv)
  out = data.frame( p    = t$p.value,
                    df   = t$parameter,
                    t    = t$statistic,
                    d    = get.os.cohens.d(sum.data, dv),
                    u    = u)
  out
}

###### -------------------------------------------------------------------------
###### prevalence statistics functions
###### -------------------------------------------------------------------------

# run.mont.frst.lvl <- function(data, N){
#   # data = 1 participants VSL data! N = number of montecarlo simulations
#   # https://arxiv.org/pdf/1512.00810.pdf
#   # see Algorithm section
#   # As 24! is in the millions, going to do a monte carlo sampling for the first level permutation
#   sub.data <- data.frame(sub=rep(data$sub[1], times=N),
#                          acc=NA)
#   sub.data$acc[1] = with(data, mean(Target.Order==Response))
#   sub.data$acc[2:N]=replicate(N-1, with(data, mean(sample(Target.Order)==Response)))
#   sub.data
# }

run.sing.mont.frst.lvl <- function(data){
  data %>% group_by(sub) %>%
    summarise(acc = mean(sample(Target.Order) == Response))
}

run.mont.frst.lvl.over.subs <- function(data,NfL){
  # feed in all VSL data and the number of monte carlo perms to run (N)
  # will apply the run.mont.frst.lvl over each subject and return a dataframe
  # note: the first permutation is the preserved orderings, as they occurred in the experiment
  nsubs = data$Nsz[1]
  idx <- gen.flvl.perms(NfL, nsubs, max(data$Trial.No))
  names(data)[names(data) == "Trial.No"] = "trial"
  
  perms <- rbind(data %>% group_by(sub) %>%
                   summarise(acc=mean(Response==Target.Order)) %>%
                   mutate(p=1),
                 inner_join(idx, data, by=c("sub", "trial")) %>% group_by(sub, p) %>%
                   summarise(acc=mean(Response==Target.Order)))
  perms %>% arrange(sub)
}

gen.flvl.perms <- function(NfL, nsubs, ntrials, neut=FALSE){
  if (neut == FALSE){
    idx <- data.frame(trial = unlist(replicate(nsubs*(NfL-1), sample(ntrials, replace=FALSE), simplify=FALSE)),
                      p = rep(2:NfL, each=ntrials*nsubs),
                      sub = rep(1:nsubs, each=ntrials, times=NfL-1))
  } else {
    idx <- data.frame(trial = rep(c(1:ntrials), times=nsubs*(NfL-1)),
                      p = rep(2:NfL, each=ntrials*nsubs),
                      sub = rep(1:nsubs, each=ntrials, times=NfL-1))
  }
  idx
}

calc.flvl.nways <- function(a, b, x, y) {
  accuracy <- (x + (b - y)) / (a + b)
  nways <- choose(a, x) * choose(b, y)
  list(accuracy, nways)
}

calc.flvl.perms <- function(answers_a, answers_b, responses_a) {
  n_answers <- answers_a + answers_b
  responses_b <- n_answers - responses_a
  if (responses_a <= answers_a) {
    a_min <- max(responses_a - answers_b, 0)
    ratio_counts <- lapply(a_min:responses_a, function (x) calc.flvl.nways(answers_a, answers_b, x, responses_a-x))
  } else {
    b_min <- max(responses_b - answers_a, 0)
    ratio_counts <- lapply(0:responses_b, function (x) calc.flvl.nways(answers_b, answers_a, x, responses_b-x))
  }
  rows <- do.call(rbind, ratio_counts)
  data.frame(acc=do.call(rbind, rows[,1]), counts=do.call(rbind, rows[,2]))
}

calc.flvl.dist <- function(answers_a, answers_b, responses_a) {
  ratio_counts <- calc.flvl.perms(answers_a, answers_b, responses_a)
  ratio_counts$p <- ratio_counts$counts / sum(ratio_counts$counts)
  ratio_counts
}

gen.samps <- function(k, nsubs) {
  data <- data.frame(shuffle=rep(2:k, nsubs),
                     sub=rep(1:nsubs, each=k-1),
                     p=unlist(lapply(1:nsubs, function (x) sample(2:k, k-1, replace=T))))
  data
}

get.perm.mins <- function(slvl.idx, flvl.idx, flvl.neut, samp.data) {   
  inner_join(slvl.idx, cbind(inner_join(flvl.idx, samp.data, by=c("sub","trial")) %>% select(-Response), inner_join(flvl.neut, samp.data, by=c("sub", "trial")) %>% select(Response)), by=c("sub", "p")) %>%
    group_by(shuffle, sub) %>%
    summarise(acc=mean(Target.Order == Response)) %>%
    group_by(shuffle) %>%
    summarise(min.acc=min(acc)) %>%
    select(min.acc)
}

get.flvl.sub.dists <- function(data, nsubs) {
  responses <- data %>% group_by(sub) %>% count(Response) %>% spread(Response, n)
  answers <- data %>% group_by(sub) %>% count(Target.Order) %>% spread(Target.Order, n)
  sar <- inner_join(answers, responses, by="sub", suffix=c(".ans", ".resp"))
  mapply(calc.flvl.dist, sar$Novel.ans, sar$Repeat.ans, sar$Novel.resp, SIMPLIFY=FALSE)
}

gen.slvl <- function(flvldists){
  unlist(lapply(flvldists, function(x) sample(x$acc, 1, prob=x$p, replace=TRUE)))
}

prev.test <- function(samp.data, alpha, k, NfL) {
  # samp.data = sample of data
  # k = the number of 2nd level perms
  # NfL = the number of first level perms
  # alpha = criteria for significance
  
  # first, run first level permutations on input dataset
  # assign a unique subject number to each data entry and
  # get the data from the first level perms
  ntrials = max(samp.data$Trial.No)
  nsubs <- length(samp.data$Trial.No)/ntrials
  samp.data$sub <- rep(1:nsubs, each=ntrials)
  names(samp.data)[names(samp.data) == "Trial.No"] = "trial"
  # compute flvl idx
  flvldists <- get.flvl.sub.dists(samp.data, nsubs)
  # Now generate a neutral idx as long and as appropriate as flvl.idx, for a subsequent cbind, prior to 
  # summarising and joining to slvl idx, summarising and then taking minimum stat.
  # Now the data is prepared generate indexing for minimum stat
  # slvl.idx = rbind(gen.samps(k, nsubs)) %>% arrange(shuffle, sub)
  # computes prevalence statistic, given a set of second level permutations (and original scores)
  # Based on: https://github.com/allefeld/prevalence-permutation/blob/master/prevalenceCore.m - lines 160-168, also
  # k = the number second level permutation you want to extract from the data
  # first select the minimum statistic from the neutral permutation
  # sort out this one 
  neut_m <- min(unlist(lapply(unique(samp.data$sub), function(x) with(samp.data[samp.data$sub == x, ], mean(Response==Target.Order)))))
  # now compute the probability of the minimum value (equation 24 of 10.1016/j.neuroimage.2016.07.040)
  puGN <- sum(c(neut_m, unlist(replicate(k, min(gen.slvl(flvldists))))) >= neut_m)/(k+1)
  # puGN <- sum(unlist(replicate(k, min(gen.slvl(flvldists)))) >= neut_m)/k
  # probability uncorrected of global null (puGN)
  # the above gives the statement of existance, next step is to evaluate against the prevalence null
  # if this is below alpha, then we would say its significant
  # so puGN is the probability that the random perms are greater than the minimum true accuracy
  # the more that are, the less likely there is an effect
  
  ####### attain upper bound on the gamma null that can be rejected (see equation 20 of 10.1016/j.neuroimage.2016.07.040)
  gamma_zero = (alpha^(1/nsubs) - puGN^(1/nsubs)) / (1 - puGN^(1/nsubs))
  # if gamma_zero = - inf, it means that all the perms were greater than the min accuracy
  # thus, the prevalence should be 0 in this case
  if (puGN > alpha) gamma_zero = 0 # not significant so not defined with a prevalence value
  
  ####### the below completes step 5a in the algorithm section of 10.1016/j.neuroimage.2016.07.040
  # # first define prevalence nulls (line 196 or equation 19)
  # null_gammas <- seq(0.01, 1, by = .01)
  # # probability uncorrected of prevalence null
  # puPN <- ((1 - null_gammas) * puGN ^ (1/nsubs) + null_gammas) ^ nsubs
  # sigMN = round(puPN,2) <= alpha
  # # return info in a dataframe
  # if (length(puPN[sigMN]) < 1){
  #   gamma_zero <- 0
  # } else {
  #   gamma_zero <- max(null_gammas[sigMN])
  # }
  
  results <- data.frame(d = gamma_zero,
                        p = puGN)
  results
}

run.prev.test <- function(data, alpha=.05, k=1000000, Np=1000){
  # data = the input data for that N and sample
  # N = the desired sample size
  # subs = the list of subjects
  # alpha = the alpha level against which to assess significance
  # k = the number of 2nd level perms
  # Np = the number of 1st level perms
  results <- prev.test(data, alpha, k, Np)
  results
}

get.ps.vsl <- function(data) {
  # run t-test and prevalence test for VSL data
  t <- run.os.t.test.sim(data)
  prev <- run.prev.test(data)

  ##### tmp variable as re-running without prev test
  #prev <- data.frame(p = NA, d = NA)
  # COLLATE OUTPUT VARIABLES
  # ----------------------------------------------------------------------------
  out <- list()

  # p 
  p <- c(t$p, prev$p)
  out$p <- p
  
  # get ef szs
  efs <- c(t$d, prev$d) 
  out$esz <- efs
  
  # get residuals from lme model
  out$esub = c(NA, NA)
  out$eRes = c(NA, NA)
  # out$df = c(t$df, NA)
  # out$t = c(t$t, NA)
  # out$u = c(t$u, NA)
  out
}

# ------------------------------------------------------------------------------
###### LME and sim functions for SRT data
#### ---------------------------------------------------------------------------
get.ps.srt <- function(data) {
  # run t-test and linear mixed effects model on SRT data
  
  options(contrasts = c("contr.sum", "contr.poly")) # set options   
  
  names(data) <- c("sub", "Nsz", "perm", "trialtype", "RT")  
  
  # number subjects
  data$trialtype <- as.factor(data$trialtype)
  nsubs = length(data$sub)/length(levels(data$trialtype))
  data$sub <- rep(1:nsubs, each = length(levels(data$trialtype)))
  
  # GET DATA STATS
  # -----------------------------------------------------------------------------
  stat_data <- data %>% select(sub, trialtype, RT) %>% 
    pivot_wider(id_cols = sub, names_from = trialtype, values_from = RT)
  stat_data <- get_RTdist_stats_dfs(stat_data)
  
  # RUN t-test
  # ----------------------------------------------------------------------------
  t <- run.t.test.sim(data, subs = unique(data$sub))
  
  # RUN LME VERSION
  # ----------------------------------------------------------------------------
  # mod <- lmer(RT ~ trialtype + (1|sub), data=data)
  # lme.an <- Anova(mod)
  # 
  # # get effect size for random effects
  # lme.peta <- summary(mod)$coefficients["trialtype1", "Estimate"]/
  #   sqrt(sum(as.data.frame(VarCorr(mod))$sdcor^2)) # get the variance of the random effects
  # 
  # COLLATE OUTPUT VARIABLES
  # -----------------------------------------------------------------------------
  out <- list()
  
  # p 
  p <- c(t$p, NA)
  out$p <- p
  
  
  # get ef szs
  efs <- c(t$d, NA) 
  out$esz <- efs
  
  # get residuals from lme model
  #df = as.data.frame(VarCorr(mod))
  out$esub = c(NA, NA)
  out$eRes = c(NA, NA)
  
  # put the distribution data into the dataframe
  list(out, stat_data)
}

# ----------------------------------------------------------------------------------------------------
###### aov and LME functions for CC data
#### -------------------------------------------------------------------------------------------------

get.ps.CC <- function(data){
  # run aov on the contextual cueing data
  # return the p value, and the conversion of peta to d for both RM Anova and LME
  # data = dataframe for testing
  # dv = name of dv (typically RT)
  
  options(contrasts = c("contr.sum", "contr.poly")) # set options  
  
  names(data) <- c("sub", "Nsz", "perm", "block", "trialtype", "RT")
  
  data$block <- as.factor(data$block)
  data$trialtype <- as.factor(data$trialtype)
  
  # number subjects and get number of reps to set SS type
  ntrials = 24
  data$trial = 1:ntrials
  type = set.SS.type(data)
  
  # now get n subs
  nsubs <- length(data$sub)/(length(levels(data$block)) * length(levels(data$trialtype)))
  # set subject identifiers as unique
  data$sub <- as.factor(rep(1:nsubs, each=length(levels(data$block)) * length(levels(data$trialtype))))

  # GET DATA STATS
  # -----------------------------------------------------------------------------
  stat_data <- data %>% select(sub, block, trialtype, RT) %>% 
    pivot_wider(id_cols = sub, names_from = c(block, trialtype), values_from = RT)
  stat_data <- get_RTdist_stats_dfs(stat_data)
  
  # RUN RM ANOVA (see anova function notes for choice)
  # -----------------------------------------------------------------------------
  # running type 3 SS to match with spss - however, data are balanced so the results are the same either way
  # however SS calcs crash when 1 sub appears > 3 times, so at theat time use SS type 1

  # adding a further tryCatch to the situ for times where type 3 fails
  an <- tryCatch({
          get_anova_table(anova_test(data=data%>%ungroup(), dv=RT, wid=sub, within=c(block,trialtype), effect.size="pes", type=type), correction='GG')
  }, error=function(cond) {
          get_anova_table(anova_test(data=data%>%ungroup(), dv=RT, wid=sub, within=c(block,trialtype), effect.size="pes", type=1), correction='GG')
  })

  # compute partial epsilon sq:
  peps <- compute_partial_epsilon_sq(an)
  peps <- peps[names(peps) %in% c("block:trialtype", "block")]
  # RUN LME VERSION 
  # Notes: this is commented out as I wound up ditching the LME side of the analyses
  # However, I kept the infrastructure so that I could send the main effect of 
  # trial type to the 'lme' structure, to save both for the analyses
  # -----------------------------------------------------------------------------
  # mod <- lmer( RT ~ block*trialtype + (1|sub), data=data )
  # lme.an <- Anova(mod)
  # 
  # # get effect size for random effects
  # lme.peta <- (summary(mod)$coefficients["block11:trialtype1", "Estimate"]-summary(mod)$coefficients["block1:trialtype1", "Estimate"])/
  #   sqrt(sum(as.data.frame(VarCorr(mod))$sdcor^2)) # get the variance of the random effects
  
  # COLLATE OUTPUT VARIABLES
  # -----------------------------------------------------------------------------
  out <- list()
  
  # p 
  p <- c(an$p[an$Effect == 'block:trialtype'], an$p[an$Effect == 'trialtype'])
  out$p <- p
  
  #peta <- c(an$pes[an$Effect == 'block:trialtype'], an$pes[an$Effect == 'trialtype']) 
  peta <- c(peps[names(peps) == "block:trialtype"], peps[names(peps) == "block"])
  out$esz <- peta
  # get residuals from lme model
  #df = as.data.frame(VarCorr(mod))
  # out$esub = c(NA, df$sdcor[df$grp=="sub"])
  # out$eRes = c(NA, df$sdcor[df$grp=="Residual"])
  out$esub = c(NA, NA)
  out$eRes = c(NA, NA)
  out
  
  list(out, stat_data)
}

# ----------------------------------------------------------------------------------------------------
###### aov and LME functions for SD data
#### -------------------------------------------------------------------------------------------------
get.ps.SD <- function(data){
  # run aov on the contextual cueing data
  # return the p value, and the conversion of peta to d
  # data = dataframe for testing
  
  options(contrasts = c("contr.sum", "contr.poly")) # set options  
  
  names(data) <- c("sub", "Nsz", "perm", "task", "trialtype", "RT")
  
  # org factors, number subjects and set SS type
  data$task <- as.factor(data$task)
  data$trialtype <- as.factor(data$trialtype)
  data$trial <- c(1:4)
  type <- set.SS.type(data)
  
  nsubs <- length(data$sub)/(length(levels(data$task)) * length(levels(data$trialtype)))
  data$sub <- as.factor(rep(1:nsubs, each=length(levels(data$task)) * length(levels(data$trialtype))))
  
  # RUN RM ANOVA
  # -----------------------------------------------------------------------------
  an <- tryCatch({
          get_anova_table(anova_test(data=data%>%ungroup(), dv=RT, wid=sub, within=c(task,trialtype), effect.size="pes", type=type))
  }, error=function(cond) {
          get_anova_table(anova_test(data=data%>%ungroup(), dv=RT, wid=sub, within=c(task,trialtype), effect.size="pes", type=1))
  })
  
  # compute partial epsilon sq:
  peps <- compute_partial_epsilon_sq(an)
  peps <- peps[names(peps) %in% c("trialtype", "task:trialtype")]

  # RUN LME VERSION
  # -----------------------------------------------------------------------------
#  mod <- lmer( RT ~ task*trialtype + (1|sub), data=data )
#  lme.an <- Anova(mod)
  
  # get effect size for random effects
#  lme.peta <- summary(mod)$coefficients["trialtype1", "Estimate"]/sqrt(sum(as.data.frame(VarCorr(mod))$sdcor^2)) # get the variance of the random effects
  
  # COLLATE OUTPUT VARIABLES
  # -----------------------------------------------------------------------------
  out <- list()
  
  # p 
  p <- c(an$p[an$Effect == "trialtype"], an$p[an$Effect == "task:trialtype"])
  out$p <- p
  
  # compute partial eta squared
  peta <- c(peps[names(peps) == "trialtype"], peps[names(peps) == "task:trialtype"])
  out$esz <- peta
  
  # get residuals from lme model
  #df = as.data.frame(VarCorr(mod))
  out$esub = c(NA, NA)
  out$eRes = c(NA, NA)
  out
}

# ----------------------------------------------------------------------------------------------------
###### aov and LME functions for AB data
#### -------------------------------------------------------------------------------------------------
# compute partial epsilon sq:
get.ps.aov.AB <- function(data){
  # run lme and aov on the dv with lag as the iv
  # return the p value
  # data = dataframe for testing - has 4 columns - sub, lag, T1, T2gT1
  # dv = name of dv (T1 or T2gT1)
  
  options(contrasts = c("contr.sum", "contr.poly")) # set options
  names(data) <- c("sub", "Nsz", "perm", "lag", "T1", "T2gT1")
  
  data$lag <- as.factor(data$lag)
  # get SS type - WARNING - HARD CODING!
  data$trial <- c(1:4)
  type <- set.SS.type(data)
  
  # number subs
  nsubs <- length(data$sub)/length(levels(data$lag))
  data$sub <- as.factor(rep(1:nsubs, each=length(levels(data$lag))))
  
  # GET DATA STATS
  # -----------------------------------------------------------------------------
  stat_data <- data %>% select(sub, lag, T2gT1) %>% 
                pivot_wider(id_cols = sub, names_from = lag, values_from = T2gT1)
  stat_data <- get_RTdist_stats_dfs(stat_data)
  
  
  # RUN RM ANOVA
  # -----------------------------------------------------------------------------
  # aov doesn't do sum of squares 3, ezANOVA = ~500 ms slower than get_anova_table
  an <- tryCatch({
          get_anova_table(anova_test(data=data%>%ungroup(), dv=T2gT1, wid=sub, within=lag, effect.size="pes", type=type))
  }, error=function(cond) {
          get_anova_table(anova_test(data=data%>%ungroup(), dv=T2gT1, wid=sub, within=lag, effect.size="pes", type=1))
  })
  
  # compute partial epsilon
  # -----------------------------------------------------------------------------
  peps <- compute_partial_epsilon_sq(an)

  # RUN LME VERSION
  # -----------------------------------------------------------------------------
  # mod <- lmer( T2gT1 ~ lag + (1|sub), data=data )
  # lme.an <- Anova(mod)
  # 
  # # get effect size for random effects
  # lme.peta <- summary(mod)$coefficients["lag1", "Estimate"]/sqrt(sum(as.data.frame(VarCorr(mod))$sdcor^2)) # get the variance of the random effects
  # 
  # COLLATE OUTPUT VARIABLES
  # -----------------------------------------------------------------------------
  out <- list()
  
  # p 
  p <- c(an$p, NA)
  out$p <- p
  
  # compute partial eta squared
  peta <- c(peps, NA) 
  out$esz <- peta
  
  # get residuals from lme model
#  df = as.data.frame(VarCorr(mod))
  out$esub = c(NA, NA)
  out$eRes = c(NA, NA)
  
  # put the distribution data into the dataframe
  
  list(out, stat_data)
}


# SET SS TYPE
set.SS.type <- function(data){
  # note! data must have the trial numbers added
  sub_list = data %>% filter(trial == 1) %>% select(sub) 
  reps=unique(sapply(sub_list$sub, function(x)(sum(x==sub_list$sub))))
  if (max(reps) > 2){
    type=1
  } else {
    type=3
  }
  type
}

# ----------------------------------------------------------------------------------------------------
###### functions to run sims ACROSS ALL TASKS
#### -------------------------------------------------------------------------------------------------
run.outer <- function(in.data, subs, N, k, j, outer_index, cores, fstem,  f, samp, seeds) {
  # this function runs the outer permutation loop
  # for 1:j permutations, for a given N
  # arguments:
  # -- in.data = dataframe for relevant analysis
  # -- subs = the list of all subject numbers
  # -- N = the number to sample from subs
  # -- k = the total number of inner permutations, 1 if using immediate sampling
  # -- j = how many outer loops? 1000 for both types of sampling
  # -- outer_index = which outer_index to run (NA == all)
  # -- cores = how many cores do you want to use?
  # -- fstem = what do you want the output files from the inner loop to be called?
  # -- f: is the reference to the function you want to run for this task
  # -- samp: sampling type: 'int' for intermediate or 'imm' for immediate
  # -- seeds: random seeds for reproducible results

  # settings for the sampling technique
  if (samp == "int") {
    replace <- FALSE
  } else if (samp == "imm") {
    replace <- TRUE
  }
  if (!is.na(outer_index)) {
    # reseed
    set.seed(seeds[outer_index])
  }
  # select an idx of N unique if int, or with replacement if imm
  sub.idx = lapply(1:j, function(x) sample.N(subs, N, x, replace=replace))
  sub.idx = do.call(rbind, sub.idx) # index for all the outerloops
  # here I mclapply over each 'parent sample' to pass into run.inner
  names(in.data)[names(in.data)=="Subj.No"] <- "sub"
  
  # so pass in all the data that will be used only on this outer loop
  outer_filter <- function(x) inner_join(sub.idx, in.data, by="sub") %>%
                              filter(perm == x) %>%
                              select(-c("Nsz", "perm"))

  if (!is.na(outer_index)) {
    run.inner(in.data=outer_filter(outer_index),
              parent.subs=sub.idx$sub[sub.idx$perm == outer_index],
              N=N,
              k=k,
              j=outer_index,
              fstem=fstem,
              f=f,
              samp=samp)
  } else {
    mclapply(1:j,
            function(x) run.inner(in.data=outer_filter(x),
                                  parent.subs=sub.idx$sub[sub.idx$perm == x],
                                  N=N,
                                  k=k,
                                  j=x,
                                  fstem=fstem,
                                  f=f,
                                  samp=samp),
            mc.cores=cores)
  }
}

run.inner <- function(in.data, parent.subs, N, k, j, fstem, f, samp) {
  # this function runs the inner permutation loop
  # for 1:k permutations, run.AB.models is run, and the output
  # collated. The results are saved in a binary file in the
  # local directory
  # arguments:
  # --in.data = dataframe (see notes of get.ps.aov.AB)
  # --parent.subs = the list of subject numbers for this parent sample
  # --N = the number to sample from parent.subs
  # --k = the total number of inner permutations, 1 if using immediate sampling
  # --j = which outer permutation are we on?
  # -- f: the function to be run
  # -- fstem: see run.outer
  # -- samp: sampling type: 'int' for intermediate or 'imm' for immediate
  
  tmp = replicate(k, run.models(in.data=in.data, subs=parent.subs, N=N, f=f, samp=samp), simplify=FALSE)
  
  # first get the standard data frame of effect size results
  tmp.fx = tmp[[1]][[1]] # note, a multi inner loop will no longer work!
                         # but I have decided I'll worry about this if the inner loop
                         # ever becomes a thing again
  out = data.frame( n    = c(N, N),
                    p    = tmp.fx$p,
                    esz  = tmp.fx$esz,
                    esub = tmp.fx$esub,
                    eRes = tmp.fx$eRes,
                    mod = c("RM-AN", "LME") )
  out$k = rep(1:k, each=2)
  # now save the out data
  fname = sprintf(fstem, N, j)
  save(out, file=fname)
  
  # now to get/save distributional data
  dist_summary <- tmp[[1]][[2]][[1]] # this code is gross, sorry world
  dist_out <- data.frame(dist_summary)
  dist_out$n <- N
  dist_out$k <- k
  dist_fstem <- str_split(fstem, pattern = ".RD", simplify = TRUE)[1,1]
  dist_fstem <- paste(dist_fstem, "_diststats.RData", sep="")
  dist_fname = sprintf(dist_fstem, N, j)
  save(dist_out, file=dist_fname)
  
  # now get and save density function data
  dens_func_dat <- tmp[[1]][[2]][2]
  names(dens_func_dat) <- sprintf("N%d_k%d", N, k)
  dens_fstem <- str_split(fstem, pattern = ".RD", simplify = TRUE)[1,1]
  dens_fstem <- paste(dens_fstem, "_DV_dens_funcs.RData", sep="")
  dens_fname <- sprintf(dens_fstem, N, j)
  save(dens_func_dat, file=dens_fname)
}

run.models <- function(in.data, subs, N, f, samp){
  # this function runs 1 simulation 
  # this function will sample the requested N with the replacement mode set by samp
  # and then apply the analysis and output to a dataframe
  # in.data = dataframe 
  # subs = the list of subject numbers to be sampled from
  # N = the number to sample from subs
  # -- f: the analysis function to be run 
  # samp: sampling type: 'int' for intermediate or 'imm' for immediate
  
  # settings for the sampling technique 
  if (samp == "int"){
    replace=TRUE
  } else if (samp == "imm"){
    replace=FALSE
  }
  
  idx = sample.N(subs=subs, N=N, k=1, replace=replace) # get the samples for this one permutation
  
  f <- eval(f) # get the function
  tmp.fx = f(inner_join(idx, in.data, by="sub"))
  tmp = tmp.fx # some re-labelling due to legacy issues
  
  tmp.fx = tmp[[1]]
  
  if (length(tmp.fx) > 4){
    out = data.frame( n    = c(N, N),
                      p    = tmp.fx$p,
                      esz  = tmp.fx$esz,
                      df = tmp.fx$df,
                      t = tmp.fx$t,
                      u = tmp.fx$u,
                      mod = c("RM-AN", "LME") )
  } else {
    out = data.frame( n    = c(N, N),
                      p    = tmp.fx$p,
                      esz  = tmp.fx$esz,
                      esub = tmp.fx$esub,
                      eRes = tmp.fx$eRes,
                      mod = c("RM-AN", "LME") )
  }
  
  list(out, tmp[[2]]) # now output the dataframe of effect results
                      # and the data distribution results
}

# ----------------------------------------------------------------------------------------------------
# Density generating functions for plotting, and for computing 95 % CI for the FFX/RFX ratio measure
# ----------------------------------------------------------------------------------------------------

dens.across.N <- function(fstem, Ns, j, min, max, spacer, dv, savekey, task, datpath, rxvnme, rxvsub, cores){
  # grab density functions for dv of choice, across all N sizes
  # save to a binary file output
  # INPUTS
  # -- fstem: see add.dens
  # -- Ns: vector of subject Ns for which there are data to grab
  # -- j: number of outer permutations
  # -- min: see add.dens (and beyond)
  # -- max: see add.dens
  # -- spacer: see get.dens
  # -- dv: which dv are we pulling data out for?
  # -- savekey: usually the initials of the paradigm, for saving output
  # -- datpath: filepath to data rxv
  # -- rxvnme: name of rxv
  # -- rxvsub: name of the sub folder in rxv
  tmp = mclapply(Ns, add.dens, fstem=fstem, j=j, min=min, max=max, 
                 dv=dv, spacer=spacer, datpath=datpath, rxvnme=rxvnme,
                 rxvsub=rxvsub, task=task, mc.cores=cores)
  d = do.call(rbind, tmp)
  fname = paste(savekey, dv, "d.RData", sep="_")
  save(d, file=paste(datpath, fname, sep=""))
  unlink(paste(datpath, "/", rxvsub, sep=""), recursive=TRUE)
}


add.dens <- function(fstem, N, j, min, max, spacer, dv, datpath, rxvnme, rxvsub, task){
  # add density functions for the ffx or rfx based values
  # for a given N, across all outer loops. 
  # output is two density functions.
  # INPUTS
  # -- fstem: fstem = filestem to be sprintf'd with N and j
  # -- N: sub sample size
  # -- j: outer loop size
  # -- min: see get.dens
  # -- max: see get.dens
  # -- spacer: see get.dens
  # -- dv: which dv do you want to know about?
  
  # first get data extracted for that N
  unzp(datpath, rxvnme, rxvsub, task, j, N)
  ds <- lapply(1:j, get.dens, fstem=paste(datpath, rxvsub, "/", task, fstem, sep=""), N=N, min=min, max=max, dv=dv, spacer=spacer)
  ds <- lapply(1:j, function(x) do.call(rbind, ds[[x]])) 
  ds <- Reduce('+', ds)
  out <- data.frame(Nsz = rep(N, each=ncol(ds)*2),
                    mod = rep(c("ffx", "rfx"), each=ncol(ds)),
                    d = c(ds["ffx",], ds["rfx",]),
                    x = seq(min, max, by=abs(max-min)/spacer))
  # remove data files
  #  unlink(paste(datpath, "/Rdat", sep=""), recursive=TRUE)
  out
}


get.dens <- function(fstem, N, j, min, max, dv, spacer){
  # this function loads the relevant data based on fstem, and calls
  # gen.dens to compute a density function
  # for the given dv, for the given N, and j, using min and max as 
  # the density range
  # INPUTS
  # -- fstem = filestem to be sprintf'd with N and j
  # -- N = subject sample size
  # -- j = outer loop permutation number
  # -- min value for density range
  # -- max value for density range
  # -- dv = which d for which to compute density
  
  load(sprintf(fstem, N, j))
  ffx = out[out$mod=="RM-AN",]
  rfx = out[out$mod=="LME",]
  if (dv == "p"){
    # if (dens_dtype == "probit"){
    #   # see https://www.rdocumentation.org/packages/gtools/versions/3.5.0/topics/logit
    #   ffx.d <- gen.dens(min, max, spacer=spacer, log( ffx[,eval(dv)] / (1- ffx[,eval(dv)] )))
    #   rfx.d <- gen.dens(min, max, spacer=spacer, log( rfx[,eval(dv)] / (1- rfx[,eval(dv)] )))
    # } else {
    ffx.d <- gen.dens(min, max, spacer=spacer, log(ffx[,eval(dv)]))
    rfx.d <- gen.dens(min, max, spacer=spacer, log(rfx[,eval(dv)]))
    # }
  } else if (dv == "esz") {
    ffx.d <- gen.dens(min, max, spacer=spacer, ffx[,eval(dv)])
    rfx.d <- gen.dens(min, max, spacer=spacer, rfx[,eval(dv)])
  } else {
    ffx.d <- rep(NA, length(seq(min, max, by=abs(max-min)/spacer)))
    rfx.d <- gen.dens(min, max, spacer=spacer, rfx[,eval(dv)])
  }
  list(ffx = ffx.d, rfx=rfx.d)
}

gen.dens <- function(min, max, spacer = 10000, data){
  # generate an UNORMALISED density function between min and max, 
  # for the vector 'data'
  # INPUTS
  # -- min: the minimum range on the x of the density function
  # -- max: same, but the max
  # -- data: the vector over which you wish to compute density
  x = seq(min, max, by=abs(max-min)/spacer)
  idx <- sapply(data, function(y) which.min(abs(x-y)))
  idx.idx = seq(1, length(idx), by=1)
  coords = t(rbind(idx, idx.idx))
  d = matrix(0, nrow=length(x), ncol=length(idx))
  # http://eamoncaddigan.net/r/programming/2015/10/22/indexing-matrices/
  d[idx + nrow(d) * (idx.idx-1)] = 1
  apply(d, 1, sum)
}


# ----------------------------------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------------------------------

plot.d <- function(d, m, yl, sc, med){
  # this function is for the mega sample
  # inputs:
  #--d: data
  #--m: model - ffx or rfx
  #--yl: ylim - c(0,3) - for example
  #--sc: how much to scale
  #--med: median sample size
  
  total_p <- d %>% group_by(mod, Nsz) %>% summarise(sum=sum(d)) 
  p <- d %>% inner_join(total_p, by=c("mod", "Nsz")) %>% mutate(dp = d/sum) %>% 
        filter(mod == eval(m)) %>% 
        ggplot(aes(x=x, y=Nsz, height=dp, group=Nsz)) +
        geom_density_ridges(stat="identity", scale=sc, rel_min_height=.01, fill=wes_palette("IsleofDogs1")[1], color=wes_palette("IsleofDogs1")[5]) +
        theme_ridges() +
        geom_hline(yintercept=which(abs(as.numeric(levels(d$Nsz))-med) == min(abs(as.numeric(levels(d$Nsz))-med)))[1],
                    linetype="dashed", color=wes_palette("IsleofDogs1")[3]) +
        ylab('N') + theme_cowplot() + xlim(yl) +
        guides(fill = FALSE, colour = FALSE) +
        ggtitle(eval(m)) 
  if (m == "RM-AN") {
    p + xlab(expression(paste(eta[p]^{2}))) +
          theme(axis.title.x = element_text(face = "italic"))
  } else if (m=="LME"){
    p + xlab(expression(tilde(paste(eta[p]^{2})))) +
          theme(axis.title.x = element_text(face = "italic"))
  } else if (m == "t") {
    p + xlab("d") +
          theme(axis.title.x = element_text(face = "italic"))
 } else if (m=="p") {
  p + xlab(bquote(gamma)) +
    theme(axis.title.x = element_text(face = "italic"))
 }
}

calc.crit.d <- function(df, N){
  max(qt(c(.025, .975), df))/sqrt(N)
}


plot.d.by.samp <- function(d, yl, sc, m){
  # this function is to plot the distributions
  # attained from a select few N, by sampling strategy
  # inputs:
  #--d: data
  #--yl: ylim - c(0,3) - for example
  #--sc: how much to scale
  total_p <- d %>% group_by(samp, mod, Nsz) %>% summarise(sum=sum(d)) 
  d %>% inner_join(total_p, by=c("mod", "Nsz", "samp")) %>% mutate(dp = d) %>% # amend if want to make a normalised distribution
    filter(mod == eval(m)) %>% 
    ggplot(aes(x=x, y=as.factor(Nsz), height=dp, group=as.factor(Nsz), fill=as.factor(samp))) +
    geom_density_ridges(stat="identity", scale=sc, rel_min_height=.0001, fill=wes_palette("IsleofDogs1")[1], color=wes_palette("IsleofDogs1")[1]) +
    theme_ridges() + facet_wrap(~as.factor(samp)) +
    xlab('effect') + ylab('N') + theme_cowplot() + xlim(yl) +
    scale_color_manual(values="white") +
    scale_fill_manual(values=wes_palette("IsleofDogs1")[1]) +
    theme(axis.title.x = element_text(face = "italic"))
}

plot.p <- function(d, m, yl, sc, med){
  # inputs:
  #--d: data
  #--m: model - ffx or rfx
  #--yl: ylim - c(0,3) - for example
  #--sc: how much to scale
  total_p <- d %>% group_by(mod, Nsz) %>% summarise(sum=sum(d)) 
  d %>% inner_join(total_p, by=c("mod", "Nsz")) %>% mutate(dp = d/sum) %>% 
    filter(mod == eval(m)) %>% 
    ggplot(aes(x=x, y=Nsz, height=dp, group=Nsz)) +
    geom_density_ridges(stat="identity", scale=sc, rel_min_height=.01, fill=wes_palette("IsleofDogs1")[1], color=wes_palette("IsleofDogs1")[5]) +
    theme_ridges()  +
    geom_hline(yintercept=which(abs(as.numeric(levels(d$Nsz))-med) == min(abs(as.numeric(levels(d$Nsz))-med)))[1],
               linetype="dashed", color=wes_palette("IsleofDogs1")[3]) +
    xlab('p') + ylab('N') + theme_cowplot() + xlim(yl) +
    guides(fill = FALSE, colour = FALSE) +
    ggtitle(eval(m)) +
    theme(axis.title.x = element_text(face = "italic")) +
    geom_vline(aes(xintercept=log(.05)), linetype="dashed", colour="black")
}


plt.fx.sz <- function(data, ylims){
  # plot effect size, given dataframe of 'n', 'measure', and 'value'
  data %>% filter(measure=="d") %>%
    ggplot(mapping=aes(x=value, y=n)) + #, fill=stat(x))) +
    geom_density_ridges(scale=2, rel_min_height=.01, fill=wes_palette("IsleofDogs1")[1], color=wes_palette("IsleofDogs1")[5]) +
    #    geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, gradient_lwd = 1.) +
    theme_ridges() +
    #    scale_fill_viridis_c(name = "value", option = "C") +
    xlab('d') + ylab('N') + theme_cowplot() + xlim(ylims) +
    #   scale_y_discrete(breaks = seq(23, 303, by = 20), labels=as.character(seq(23, 303, by = 20))) +
    guides(fill = FALSE, colour = FALSE) +
    ggtitle(paste(data$model[1])) +
    theme(axis.title.x = element_text(face = "italic"))
}


plt.ps <- function(data, xlims, rel_min_height){
  # same as plt.fx.sz but for p values.
  data %>% filter(measure=="p") %>%
    ggplot(mapping=aes(x=value, y=n)) + #, fill=stat(x))) +
    geom_density_ridges(scale=2, rel_min_height=rel_min_height, fill=wes_palette("IsleofDogs1")[1], color=wes_palette("IsleofDogs1")[5]) + # 
    theme_ridges() +
    xlab('p') + ylab('N') + theme_cowplot() + xlim(xlims) +
    geom_vline(aes(xintercept=log(.05)), linetype="dashed") +
    guides(fill = FALSE, colour = FALSE) +
    ggtitle(paste(data$model[1])) +
    theme(axis.title.x = element_text(face = "italic"))
}


RmType <- function(string) { # remove 1st label from facet_wrap
  sub("._", "", string)
}

plt.rfx <- function(data, xlims){
  # same as plt.fx.sz but for p values.
  data %>% pivot_longer(names(data)[!names(data) %in% c("n","model")], names_to = "rfx", values_to="var") %>%
    drop_na() %>%
    ggplot(mapping=aes(x=var, y=n)) + #, fill=stat(x))) +
    geom_density_ridges(scale=2, rel_min_height=.01, fill=wes_palette("IsleofDogs1")[5], color=wes_palette("IsleofDogs1")[4]) +
    #    geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, gradient_lwd = 1.) +
    theme_ridges() +
    #    scale_fill_viridis_c(name = "value", option = "C") +
    xlab(expression(sigma)) + ylab('N') + theme_cowplot() + xlim(xlims) + 
    #    scale_y_discrete(breaks = seq(23, 303, by = 20), labels=as.character(seq(23, 303, by = 20))) +
    facet_wrap(~model*rfx) +
    guides(fill = FALSE, colour = FALSE) +
    theme(axis.title.x = element_text(face = "italic"))
}

# ----------------------------------------------------------------------------------------------------
# More data wrangles
# ----------------------------------------------------------------------------------------------------
unzp <- function(datpath, rxvnme, rxvsub, task, j, subN) {
  # function to unzip specific files from the data archive
  # :: datpath = where is the data?
  # :: rxnme = name of zipped (rxiv) file e.g. IMMCC
  # :: rxvsub = sub folder name in rxv
  # :: task = which task do you want to extract data for? 
  # e.g. "CC" or "imm_CC"
  # :: j = total number of permutations/parent sets
  # :: subs = the sub Ns used in the perms
  if (length(j) == 1){
    fnums <- 1:j
  } else {
    fnums <- j
  }
  # amended to get exact filenames for each sub num
  tmp = unzip(zipfile=paste(datpath, rxvnme, sep=""), list = TRUE)
  tmp = tmp$Name[!is.na(str_extract(tmp$Name, paste("N-", subN, "_", sep = "")))]
  fnums = tmp
  
  print_and_unzip <- function(y, datpath, rxvnme, sep = ""){
    
    tryCatch({unzip(paste(datpath, rxvnme, sep=""), 
              file=y,
              exdir=datpath)
      },
      warning = function(cond){
        message(paste(y))
      })
  }
  lapply(fnums, function(x) print_and_unzip(x, datpath=datpath, rxvnme = rxvnme))
}

d2r <- function(d, model_name){
  # function to convert d's to r^2
  tmp = d %>% filter(mod == model_name) %>%
         group_by(Nsz, mod) %>%
         mutate(x = (x/(sqrt((x^2) + 4)))^2) %>% # Cohen 1988 equation 2.2.6
         group_by(Nsz, mod, as.factor(x)) %>%
         summarise(d=sum(d)) %>%
         ungroup()
  names(tmp)[names(tmp) == "as.factor(x)"] = 'x'
  tmp$x <- as.numeric(paste(tmp$x))
  rbind(tmp,
        d %>% filter(mod != model_name) %>%
        filter(x>0))
  }

compute_partial_epsilon_sq <- function(an){
  # given 'an' [output from get_anova_table(anova_test(...))], compute
  # partial epsilon sq for each effect
  # uses equation of appendix A of Albers & Lakens 2018, see: 
  # https://www.sciencedirect.com/science/article/pii/S002210311630230X
  ## inputs:
  # -- an [output from get_anova_table(anova_test(...))]
  ## outputs:
  # -- 
  
  sapply(unique(an$Effect), function(x) (an$F[an$Effect == x] - 1) / 
           (an$F[an$Effect == x] + (an$DFd[an$Effect == x]/an$DFn[an$Effect == x])))
}

# ----------------------------------------------------------------------------------------------------
###### functions to get distributional info of RT/accuracy data
#### -------------------------------------------------------------------------------------------------

get_RTdist_stats <- function(x, cond){
  # given a vector x, return distribution
  # moments
  # Kwargs:
  #   -- x : a vector of real values
  #   -- cond : str, condition name e.g. "lag2"

  mu <- mean(x)
  med <- median(x)
  sigma <- var(x)
  skew <- (mu - med)/sqrt(sigma) # see https://en.wikipedia.org/wiki/Skewness
  
  # now kurtosis - see https://en.wikipedia.org/wiki/Kurtosis - Moore's interpretation
  Z <- (x - mu) / sqrt(sigma)[1]
  k <- var(Z^2)+1 # excessive kurtosis < -3 ok 3 > excessive 
  
  # now make the ecdf so we can recreate the distribution later
  d <- density(x)
  
  list(data.frame(cond = cond,
                  mu = mu,
                  med = med, 
                  sigma = sigma,
                  skew = skew,
                  k = k),
                  d = d)
}

get_RTdist_stats_dfs <- function(stat_data){
  # return a data frame containing the stats info for the dependent variable (dv) distributions
  # return a list of the density info
  # Kwargs
  # -- stat_data - a dataframe of the condition data in WIDEFORM, containing only the columns 'sub' and the
  #                  other variables of interest
  # Returns
  # --  lists, one element contains the dataframe
  #            the other contains the densities
  
  stat_conds <- colnames(stat_data)[colnames(stat_data)!="sub"]
  stats <- lapply(stat_conds, function(y) get_RTdist_stats(unlist(stat_data[,y]), y))
  stats_df <- do.call(rbind, lapply(c(1:length(stats)), function(x) stats[[x]][[1]]))
  stats_dens <- do.call(rbind, lapply(c(1:length(stats)), function(x) stats[[x]][[2]]))
  rownames(stats_dens) <- stat_conds
  
  list(stats_df, stats_dens)
}













