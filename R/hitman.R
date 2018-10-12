#' High-throughput mediation analysis
#'
#' High-throughput mediation analysis to test if rows of \code{M} mediate the effect of exposure \code{E} on outcome
#' \code{Y}.
#' 
#' @param E A vector representing exposure; can be numeric or a character/factor with nominal groups. 
#' @param M A numeric matrix-like data object with one row per feature and one column per sample representing the mediator.
#' Must have more than one feature.
#' @param Y A numeric vector of \code{length(E)} of outcomes.
#' @param covariates Numeric vector or matrix of covariates.
#' @param verbose Logical indicating if messages should be reported.
#' @details If \code{E} and \code{Y} have names, and \code{M} has colnames, they should all match.
#' @export

#can add covariates in future
hitman <- function(E, M, Y, covariates=NULL, verbose=FALSE){
  stopifnot(length(Y)==ncol(M), is.numeric(Y), names(Y)==colnames(M), length(E)==ncol(M), names(E)==colnames(M), 
            length(E) > 0, nrow(M) > 1)
  
  if (is.numeric(E)){
    if (verbose) message("E treated as continuous numeric vector.")
    batch <- NULL
    if (stats::var(E, na.rm=TRUE)==0) stop("E treated as numeric, but has no variance.")
    #ok if covariates is NULL
    my.covar <- cbind(E, covariates)

  } else {
    if (verbose) message("E treated as an unordered factor.")
    ngrps <- length(unique(E))
    if (ngrps==1) stop("E is not numeric, but has only one group.")
    if (any(is.na(E))) stop("E is not numeric, but has an NA.")
    #leave covariates
    batch <- E
    my.covar <- covariates
  }
  
  #if mult grps, get an F-stat -> no ey.sign
  ey.sign <- NA
  if (is.numeric(E) || ngrps==2){
    #Y treated as gene expression -> dependent variable
    tt.ey <- limma_dep(object=Y, Y=E, covariates=covariates, prefix="EY")
    if (tt.ey$EY.t != 0){
      ey.sign <- sign(tt.ey$EY.t)
      if (abs(tt.ey$EY.t) < 1){
        warning("E and Y are weakly associated, so mediation may not be meaningful.")
      }
    } else {
      #leave ey.sign as NA
      warning("E and Y are not associated, so not accounting for direction of association.")
    }
  }
  
  tt.em <- limma_dep(object=M, Y=E, covariates=covariates, prefix="EM")
  tt.my <- limma_pcor(object=M, phenotype=Y, batch=batch, covariates=my.covar, prefix="MY")
  ret <- cbind(tt.em[rownames(tt.my),], tt.my)
  
  #modify separate columns, to keep stats of two-sided tests for inspection.
  if (!is.na(ey.sign)){
    ret <- cbind(EM_dir.p=ret$EM.p, MY_dir.p=ret$MY.p, ret)
    p.cols <- c("EM_dir.p", "MY_dir.p")
    ret <- modify_hitman_pvalues(tab=ret, overall.sign = ey.sign, p.cols=p.cols)
  } else {
    p.cols <- c("EM.p", "MY.p")
  }
  
  EMY.p <- apply(ret[,p.cols], MARGIN=1, FUN=function(v){
    max(v)^2
  })
  EMY.FDR <- stats::p.adjust(EMY.p, method="BH")
  ret <- cbind(EMY.p, EMY.FDR, ret)
  ret <- ret[order(ret$EMY.p),]
  return(ret)
}
