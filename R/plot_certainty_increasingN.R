## K. Garner 2022
## plotting bits of information attained at each N for each task
###############################################################
rm(list=ls())

############ LOAD PACKAGES ##################
library(tidyverse)
library(gridExtra)

########### DEFINE VARIABLES ################
Ns <- round(exp(seq(log(13), log(313), length.out = 20)))
tNs <- length(Ns)
tfx <- 6 # how many fx are you plotting (across tasks)

bw <- 30
vers <- "defrngon13"

############ LOAD H'S AND TIDY ########################
load(file= sprintf("../data/AB/AB_H_bn%d_v%s.RData", bw, vers))
load(file=sprintf("../data/SRT/SRT_H_bn%d_v%s.RData", bw, vers))
load(file=sprintf("../data/SD/SD_ME_H_bn%d_v%s.RData", bw, vers))
load(file=sprintf("../data/SD/SD_INT_H_bn%d_v%s.RData", bw, vers))
load(file=sprintf("../data/CC/CC_ME_H_bn%d_v%s.RData", bw, vers))
load(file=sprintf("../data/CC/CC_INT_H_bn%d_v%s.RData", bw, vers))

allH <- tibble(H = c(AB_H, SRT_H, SD_ME_H, SD_INT_H, CC_ME_H, CC_INT_H),
               cond = rep(c("AB", "SRT", "MT_me", "ME_int", "CC_ME", "CC_int"), each=tNs),
               N = rep(Ns, times=tfx))

allH %>% ggplot(aes(x=N, y=H, group=cond, colour=cond)) +
  geom_line()

########### now get the certainty gain by adding an extra participant ##############
CG <- allH %>% mutate(one_p = H/N) %>% 
               mutate(Hp = one_p*(N+1),
                      CG = (Hp - H)/1000)
CG %>% ggplot(aes(x=N, y=CG, group = cond, colour = cond)) +
  geom_line()
