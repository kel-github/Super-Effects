### written by K. Garner, Jan 2021
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

### Call plotting functions for observed effect sizes and p values

rm(list=ls())

# ----------------------------------------------------------------------------------------------------
# load packages and source function files
# ----------------------------------------------------------------------------------------------------

library(tidyverse) # for data wrangling
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(lme4) # for mixed effects modelling
library(ggridges)
library(car)
library(parallel)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting

# ----------------------------------------------------------------------------------------------------
# define session variables
# ----------------------------------------------------------------------------------------------------

task = "SD"
med = 24
d_scale_ffx = 2
d_scale_rfx = 2
p_scale_ffx = 2
p_scale_rfx = 2
px_rng_d = c(0,3)
px_rng_p_ffx = c(-1000,0)
px_rng_p_rfx = c(-1000,0)
width = 8
height = 8

# ----------------------------------------------------------------------------------------------------
# LIST OF SETTINGS
# ----------------------------------------------------------------------------------------------------
# AB
# med=24, d_scale_ffx = 2, d_scale_rfx = 2, p_scale_ffx = 2, p_scale_rfx = 2, p_rng_d = c(0,5), px_rng_p_ffx = c(-800,0), px_rng_p_rfx = c(-800,0)

# CC 
# med = 23, d_scale_ffx = 2, d_scale_rfx = 2, p_scale_ffx = 2, p_scale_rfx = 2, p_rng_d = c(0,1), px_rng_p_ffx = c(-50,0), px_rng_p_rfx = c(-50,0)

# SRT
# med = 39, d_scale_ffx = 2, d_scale_rfx = 2, p_scale_ffx = 2, p_scale_rfx = 2, p_rng_d = c(0,3), px_rng_p_ffx = c(-400,0), px_rng_p_rfx = c(-400,0)

# SD
# med = 24, d_scale_ffx = 2, d_scale_rfx = 2, p_scale_ffx = 2, p_scale_rfx = 2, p_rng_d = c(0,3), px_rng_p_ffx = c(-1000,0), px_rng_p_rfx = c(-1000,0)

# ----------------------------------------------------------------------------------------------------
# define datas and load ds
# ----------------------------------------------------------------------------------------------------

fnames = c(paste("../data/", task, "_d", "_d.RData", sep=""), paste("../data/", task, "_p", "_d.RData", sep=""))
load(fnames[1])

# ----------------------------------------------------------------------------------------------------
# define factors and plot
# ----------------------------------------------------------------------------------------------------

d$Nsz <- as.factor(d$Nsz)
d$mod <- as.factor(d$mod)


ffx.d <- plot.d(d, "ffx", px_rng_d, d_scale_ffx, med)
rfx.d <- plot.d(d, "rfx", px_rng_d, d_scale_rfx, med)

# ----------------------------------------------------------------------------------------------------
# load p, define factors and plot
# ----------------------------------------------------------------------------------------------------

load(fnames[2])


d$Nsz <- as.factor(d$Nsz)
d$mod <- as.factor(d$mod)

if (task == "SD") d <- d %>% filter(as.numeric(Nsz) < 17) # remove sample sizes saturated at 0

ffx.p <- plot.p(d, "ffx", px_rng_p_ffx, p_scale_ffx, med)
rfx.p <- plot.p(d, "rfx", px_rng_p_rfx, p_scale_rfx, med)

p = plot_grid(ffx.d, rfx.d, ffx.p, rfx.p, labels=c('A', 'B', 'C', 'D'), label_size = 12, align="v")
# #p # print out the plot so you can see it
p = p + ggsave(paste("../images/", task, ".png", sep=""), width = width, height = height, units="in")