## Use this script to plot a metric across N, grouped by task
## Note: I wrote an ran this code to create the two figs as I had the document 
## open.
## K. Garner 2022

# ----------------------------------------------------
## common settings
# ----------------------------------------------------
library(RColorBrewer)
fx <- c("AB", "tc", "tc*m", "SRT", "b*c", "c")
fig_palette <- brewer.pal(length(fx), "Dark2")
# now make these colours transparent

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

fig_palette <- unlist(lapply(fig_palette, t_col))
names(fig_palette) <- fx

w <- 10/2.54 # to get inches
h <- w

# ----------------------------------------------------
## qq-ratio plot
# ----------------------------------------------------
# check data is loaded
if (!exists("qqrats_acrss_tsks")){
  warning("load the data dummy")
} else {
  
  ylims <- c(0, 7)
  n <- unique(qqrats_acrss_tsks$N)
  
  pdf(paste("../images/", "EPS", "_", "all_tasks", "_", "qqratio", ".pdf", sep = ""),
      width = w, height = h) 

  plot(n, qqrats_acrss_tsks$qqrat[qqrats_acrss_tsks$task == fx[1]], 
       type="l", lty=1, lwd=3, col = fig_palette[1],
       xaxt = "n", yaxt = "n", ylim = ylims,
       xlab = "N", ylab = "qq-ratio", bty = "n")
  for (i in 2:length(fx)){
    points(n, qqrats_acrss_tsks$qqrat[qqrats_acrss_tsks$task == fx[i]],
           type="l", lty=1, lwd=3, col=fig_palette[i])
  }
  axis(1, at=n, labels = paste(n))
  axis(2, at=ylims, labels=as.character(ylims))
  
  # add legend here
  legend(x=150, y=ylims[2],
         legend = fx, col = fig_palette,
         lty=1, lwd=3, bty="n")
  dev.off()

}
      

# ----------------------------------------------------
## p-hit plot
# ----------------------------------------------------
tasks <- c("AB", "CC", "SRT", "SD")
ftmplt = "../data/%s/%s_phitbest.RData"

########### LOAD DATASETS ################
get_prct_dats <- function(ftmplt, task){
  load(sprintf(ftmplt, task, task))
  if(task == "AB" | task == "SRT"){ 
    out$mod <- NA
    out$mod <- task
  }
  out
}

plt_dat <- do.call(rbind, lapply(tasks, get_prct_dats, ftmplt=ftmplt))
plt_dat$mod[plt_dat$mod == "tt"] <- "tc"
plt_dat$mod[plt_dat$mod == "ta*tt"] <- "tc*m"

if (!exists("plt_dat")){
  warning("get the p data dum dum")
} else {
  
  ylims <- c(0, 1)
  n <- unique(plt_dat$n)
  
  pdf(paste("../images/", "EPS", "_", "all_tasks", "_", "phit", ".pdf", sep = ""),
      width = w, height = h) 
  
  plot(n, plt_dat$p_upper[plt_dat$mod == fx[1]], 
       type="l", lty=1, lwd=3, col = fig_palette[1],
       xaxt = "n", yaxt = "n", ylim = ylims,
       xlab = "N", ylab = "p(hit|N)", bty = "n")
  for (i in 2:length(fx)){
    points(n, plt_dat$p_upper[plt_dat$mod == fx[i]],
           type="l", lty=1, lwd=3, col=fig_palette[i])
  }
  axis(1, at=n, labels = paste(n))
  axis(2, at=c(0, 0.5, 1), labels=paste(c(0, 0.5, 1)))
  
  # add legend here
  legend(x=189, y=0.7,
         legend = fx, col = fig_palette,
         lty=1, lwd=3, bty="n")
  dev.off()
  
}