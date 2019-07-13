#' Apply \code{fisher_enrichment_test} to simulated data analyzed by \code{hitman}
#' 
#' Apply \code{fisher_enrichment_test} to simulated data analyzed by \code{hitman}.
#' 
#' 
#' @param effect.v Numeric vector of log fold-changes or percent of phentotypes to add.
#' @inheritParams roast_contrasts
#' @inheritParams hitman
#' @inheritParams sim_barfield


# test mult pwys for efficiency
sim_fisher_hitman <- function(G, E, effect.v=c(0, 0.2), alpha=0.05, nsim=10**3, seed=1, verbose=TRUE){
  
  all.feats <- rownames(feat.tab)
  prop.sig.mat <- matrix(NA, nrow=nsim, ncol=length(effect.v), 
                         dimnames=list(paste0("sim", 1:nsim), paste0("eff_", effect.v)))
  
  set.seed(seed)
  
  err <- stats::rnorm(n=length(E), sd=stats::sd(E)/2)
  phenotype <- E + err
  
  for (sim in 1:nsim){
    for (ev in effect.v){
      obj.test <- matrix(stats::rnorm(n=length(all.feats)*length(grp)), ncol=length(grp), nrow=length(all.feats),
                         dimnames=list(all.feats, names(grp)))
      if (ev > 0){
        obj.test <- t(t(obj.test) + phenotype + ev*err)
      }
      
      feat.tab <- hitman(M=obj.test, E=E, Y=phenotype)[, 1, drop=FALSE]
      gset <- list(hm=rownames(feat.tab)[feat.tab$EMY.p < alpha])
      fet <- fisher_enrichment_test(sig.sets = gset, G=G, feat.tab = feat.tab, min.nfeats = 3)
      prop.sig.mat[sim, paste0("eff_", ev)] <- mean(fet$hm.p < alpha)
    }
    if (sim %% 100 == 0) cat("sim: ", sim, "\n")
  }
  (prop.sig.mat <- rbind(avg=colMeans(prop.sig.mat), prop.sig.mat))
}