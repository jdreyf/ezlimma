#' Apply \code{fisher_enrichment} to simulated data
#' 
#' Apply \code{fisher_enrichment} to simulated data.
#' 
#' @param effect.v Numeric vector of log fold-changes or percent of phentotypes to add.
#' @param alpha Alpha level.
#' @param nsim Number of simulations.
#' @param test.ft Compare to \code{\link[stats]{fisher.test}}.
#' @param seed Random seed for reprodubility.
#' @param verbose Logical; should the number of simulations be printed every 100 simulations?
#' @inheritParams roast_contrasts

# competitive test, so can only test one pathway at a time
# need high alpha for less variance in estimate
sim_fisher <- function(G, feat.tab, grp, effect.v=c(0, 0.2), alpha=0.2, nsim=99, test.ft=TRUE, seed=1, verbose=FALSE){
  stopifnot(nrow(feat.tab) > 10)
  all.feats <- rownames(feat.tab)
  g1 <- G[[1]]$genes
  if (is.null(names(grp))) names(grp) <- paste0("s", 1:length(grp))
  
  prop.sig.mat <- matrix(NA, nrow=nsim, ncol=length(effect.v), 
                         dimnames=list(paste0("sim", 1:nsim), paste0("eff_", effect.v)))
  contrast.v <- c(vs=paste(unique(grp)[2], unique(grp)[1], sep="-"))
  
  set.seed(seed)
  for (sim in 1:nsim){
    for (ev in effect.v){
      obj.test <- matrix(stats::rnorm(n=length(all.feats)*length(grp)), ncol=length(grp), nrow=length(all.feats),
                         dimnames=list(all.feats, names(grp)))
      if (ev > 0){
        obj.test[g1, grp == unique(grp)[2]] <- obj.test[g1, grp == unique(grp)[2]] + ev
      }
      
      feat.tab.tmp <- limma_contrasts(object=obj.test, grp=grp, contrast.v=contrast.v)
      # rownames are maintained
      feat.tab.tmp$vs.p <- two2one_tailed(feat.tab.tmp, alternative = "greater")
      feat.tab.tmp <- feat.tab.tmp[order(feat.tab.tmp$vs.p),]
      
      sig.set <- list(top=rownames(feat.tab.tmp)[1:10])
      
      fet <- fisher_enrichment(sig.sets = sig.set, G=G, feat.tab = feat.tab, name=NA)
      #test
      if (test.ft){
        ftp <- stats::fisher.test(all.feats %in% sig.set[[1]], all.feats %in% g1, alternative = "greater")$p.value
        stopifnot(all.equal(fet["pwy1", "top.p"], ftp))
      }
      prop.sig.mat[sim, paste0("eff_", ev)] <- fet["pwy1", "top.p"] < alpha
    }
    if (verbose & sim %% 100 == 0) cat("sim: ", sim, "\n")
  }
  (prop.sig.mat <- rbind(avg=colMeans(prop.sig.mat), prop.sig.mat))
}