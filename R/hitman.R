#' High-throughput mediation analysis
#'
#' High-throughput mediation analysis to test if rows of \code{M} mediate the effect of exposure \code{E} on outcome
#' \code{Y}.
#' 
#' @param E A vector representing exposure; can be a character with groups or numeric.
#' @param M A numeric matrix-like data object with one row per feature and one column per sample representing the mediator.
#' @param Y A numeric vector of \code{length(E)} of outcomes.
#' @param covar A matrix-like data object with rows corresponding to samples and columns to covariates. If \code{covar}
#' has rownames, they should match \code{names(E)}.
#' @export
#' @import stats

hitman <- function(E, M, Y, covar=NULL){
  stopifnot(length(Y)==ncol(M), is.numeric(Y), names(Y)==colnames(M), length(E)==ncol(M), names(E)==colnames(M))
  
  n <- length(E)

  tt.em <- limma_dep(object=M, Y=E, covar=covar, prefix="EM")

  if (!is.null(covar)){
    covar2 <- data.frame(E, covar)
  } else {
    covar2 <- data.frame(E)
  }
  
  tt.my <- limma_pcor(object=M, phenotype=Y, covar=covar2, prefix="MY")
  
  ret <- cbind(tt.my, tt.em[rownames(tt.my),])
  comb.p <- apply(ret[,c("EM.p", "MY.p")], MARGIN=1, FUN=max)
  comb.FDR <- stats::p.adjust(comb.p, method="BH")
  ret <- cbind(comb.p, comb.FDR, ret)
  ret <- ret[order(ret$comb.p),]
  return(ret)
}
