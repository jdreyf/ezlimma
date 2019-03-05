#' High-throughput mediation analysis
#'
#' High-throughput mediation analysis to test if rows of \code{M} mediate the effect of exposure \code{E} on outcome
#' \code{Y}.
#' 
#' @param E A numeric vector or matrix of exposures. Can be a design matrix returned by 
#' \code{\link[stats]{model.matrix}} or \code{\link[ezlimma]{batch2design}}.
#' @param M A numeric matrix-like data object with one row per feature and one column per sample of mediators.
#' Must have more than one feature.
#' @param Y A numeric vector of \code{length(E)} of outcomes.
#' @param covariates Numeric vector or matrix of covariates.
#' @return Data frame with columns
#' \describe{
#' \item{EMY.p}{Overall p-value for mediation}
#' \item{EMY.FDR}{Overall FDR for mediation}
#' \item{EM_dir.p}{p-value for E-->M accounting for direction of mediation}
#' \item{MY_dir.p}{p-value for M-->Y accounting for direction of mediation}
#' \item{EM.p}{p-value for E-->M, not accounting for direction}
#' \item{EM.t or EM.F}{F-statistic or t-statistic for E-->M, not accounting for direction}
#' \item{MY.p}{p-value for M-->Y, not accounting for direction}
#' \item{MY.slope}{slope of regression for M-->Y, not accounting for direction}
#' }
#' and annotation.
#' @details If \code{E} and \code{Y} have names, and \code{M} has colnames, they should all match. \code{E} and \code{Y}
#' cannot have \code{NA}s.
#' 
#' If \code{E} is a matrix, the direction of E-->Y is not evaluated, which reduces power and may inflate the false
#' positive rate. It is preferable to test a column of \code{E} with other columns as covariates.
#' @export

#can add covariates in future
hitman <- function(E, M, Y, covariates=NULL){
  stopifnot(limma::isNumeric(E), limma::isNumeric(M), limma::isNumeric(Y), !is.na(E), !is.na(Y), length(E) > 0, 
            nrow(M) > 1, nrow(as.matrix(E))==ncol(M), length(Y)==ncol(M), names(Y)==colnames(M))
  
  if (ncol(as.matrix(E))==1){
    stopifnot(colnames(M)==names(E))
  } else {
    stopifnot(colnames(M)==rownames(E))
  }
  
  if (any(apply(X=as.matrix(E), MARGIN=2, FUN=stats::var, na.rm=TRUE) == 0)){
    stop("E treated as numeric, but has one or more columns with no variance.")
  }
  #ok if covariates is NULL
  my.covar <- cbind(E, covariates)
  
  #test EY; return ey.sign & weak assoc warning
  #Y treated as gene expression -> dependent variable
  tt.ey <- limma_dep(object=Y, Y=E, covariates=covariates, prefix="EY")
  if (tt.ey$EY.p > 0.99){
    stop("E and Y are not associated.")
  }
  if (tt.ey$EY.p > 0.1){
    warning("E and Y are not associated, so mediation may not be meaningful.")
  }
  
  if (ncol(as.matrix(E)) == 1){
    ey.sign <- sign(tt.ey$EY.t)
  } else {
    #if mult grps, get an F-stat -> no ey.sign
    ey.sign <- NA
  }
  
  #change order of columns so it's consistent with c("MY.p", "MY.slope")
  tt.em <- limma_dep(object=M, Y=E, covariates=covariates, prefix="EM")[,2:1]
  tt.my <- limma_pcor(object=M, phenotype=Y, covariates=my.covar, prefix="MY")
  tt.my <- tt.my[,setdiff(colnames(tt.my), "MY.FDR")]
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
