# written by K. Garner 2021
# perform inferential stats, given 'res' list/structure

do_KL_and_spearmans <- function(kl_inputs, res){
  # compute the correlation between N and the divergence
  # ----------------------------------------------------
  # kwargs
  # -- kl_inputs: a list comprising of the following fields
  # origin <- e.g. "313"
  # dv <- "stats_fx", "stats_p", "stats_sig"
  # sub_Ns <- names of subject groups (Ns) - pasted as string
  # -- res: stored results from get_statistics.R
  kls <- calc_kl_sing_origin(kl_inputs, res)
  kls <- lapply(colnames(kls), function(x) do.call(rbind, kls[,x]))
  lapply(kls, function(z) cor.test(x=z, y=as.numeric(kl_inputs$sub_Ns)))
}

do_z_tests <- function(z_inputs, res){
  # compute z-tests between each density of interest
  # ----------------------------------------------------  
  # kwards:
  # -- z_inputs: a list comprising
  # ---- x: a string determining the z-tests to be 
  # conducted - "meta" or "model"
  # ----- sub_Ns: names of subject groups (Ns) - pasted as string
  # -- res: stored results from get_statistics.R
  x <- z_inputs$x
  sub_Ns <- z_inputs$sub_Ns
  
  # define z function
  comp_z <- function(mu_x, sd_x, mu_y, sd_y, n=1000^2){
    (mu_x - mu_y) / sqrt( (sd_x/n) + (sd_y/n) )
  }
  
  if (x == "meta"){
    
    zs <- do.call(rbind, lapply(sub_Ns, function(z)
                                lapply(c("RM-AN", "LME"), function(k)
                                comp_z(mu_x = res[,"stats_fx"][[z]][[1, k]],
                                       sd_x = res[,"stats_fx"][[z]][[2, k]],
                                       mu_y = res[,"stats_sig"][[z]][[1, k]],
                                       sd_y = res[,"stats_sig"][[z]][[2, k]]))))
    colnames(zs) <- c("RM-AN", "LME")
    
    ps <- do.call(cbind, lapply(c("RM-AN", "LME"), function(x)
                   2*pnorm(q=abs(do.call(rbind, zs[,x])), lower.tail=FALSE)))
    ps <- matrix(apply(ps, 2, p.adjust, method="fdr"), nrow=20)
  } else if (x == "model"){
    
    zs <- do.call(rbind, lapply(sub_Ns, function(z)
                            comp_z(mu_x = res[,"stats_fx"][[z]][[1, "RM-AN"]],
                                 sd_x = res[,"stats_fx"][[z]][[2, "RM-AN"]],
                                 mu_y = res[,"stats_fx"][[z]][[1, "LME"]],
                                 sd_y = res[,"stats_fx"][[z]][[2, "LME"]])))
    ps <- 2*pnorm(q=abs(zs), lower.tail=FALSE)
    ps <- matrix(p.adjust(ps, method = "fdr"))
  }
 
   list(zs = zs, ps = ps)
 
}


