#' High-throughput mediation analysis
#'
#' High-throughput mediation analysis to test if rows of \code{M} mediate the effect of exposure \code{E} on outcome
#' \code{Y}.
#' 
#' @param E A vector representing exposure; can be a character with groups or numeric.
#' @param M A numeric matrix-like data object with one row per feature and one column per sample representing the mediator.
#' @param Y A numeric vector of \code{length(E)} of outcomes.
#' @param covariates Numeric vector or matrix of covariates.
#' @param verbose Logical indicating if messages should be reported.
#' @export
#' @import stats

#can add covariates in future
hitman <- function(E, M, Y, covariates=NULL, verbose=FALSE){
  stopifnot(length(Y)==ncol(M), is.numeric(Y), names(Y)==colnames(M), length(E)==ncol(M), names(E)==colnames(M))
  
  if (!is.numeric(E)){
    if (verbose) cat("E being treated as a nominal group variable.\n")
    #leave covariates
    batch <- E
    my.covar <- covariates
  } else {
    if (verbose) cat("E being treated as a numeric variable.\n")
    batch <- NULL
    if (is.null(covariates)){
      my.covar <- E
    } else {
      my.covar <- cbind(E, covariates)
    }
  }

  tt.em <- limma_dep(object=M, Y=E, covariates=covariates, prefix="EM")
  
  tt.my <- limma_pcor(object=M, phenotype=Y, batch=batch, covariates=my.covar, prefix="MY")
  
  ret <- cbind(tt.my, tt.em[rownames(tt.my),])
  comb.p <- apply(ret[,c("EM.p", "MY.p")], MARGIN=1, FUN=function(v){
    max(v)^2
  })
  comb.FDR <- stats::p.adjust(comb.p, method="BH")
  ret <- cbind(comb.p, comb.FDR, ret)
  ret <- ret[order(ret$comb.p),]
  return(ret)
}
