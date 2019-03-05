#' High-throughput mediation analysis
#'
#' High-throughput mediation analysis to test if rows of \code{M} mediate the effect of exposure \code{E} on outcome
#' \code{Y}.
#' 
#' @param E A numeric vector of exposures.
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
  stopifnot(limma::isNumeric(E), limma::isNumeric(M), limma::isNumeric(Y), !is.na(E), !is.na(Y), var(E) > 0, 
            nrow(M) > 1, length(E)==ncol(M), length(Y)==ncol(M), colnames(M)==names(E), names(Y)==colnames(M))
  #ok if covariates is NULL
  my.covar <- cbind(E, covariates)
  
  #test EY; return ey.sign & weak assoc warning
  #Y treated as gene expression -> dependent variable
  #not sure what this does when object is a vector
  tt.ey <- limma_dep(object=Y, Y=E, covariates=covariates, prefix="EY")
  if (tt.ey$EY.p > 0.99){
    stop("E and Y are not associated.")
  }
  if (tt.ey$EY.p > 0.1){
    warning("E and Y are not associated, so mediation may not be meaningful.")
  }
  ey.sign <- sign(tt.ey$EY.t)
  
  #change order of columns so it's consistent with c("MY.p", "MY.slope")
  tt.em <- limma_dep(object=M, Y=E, covariates=covariates, prefix="EM")[,2:1]
  tt.my <- limma_pcor(object=M, phenotype=Y, covariates=my.covar, prefix="MY")
  tt.my <- tt.my[,setdiff(colnames(tt.my), "MY.FDR")]
  ret <- cbind(tt.em[rownames(tt.my),], tt.my)
  
  #modify separate columns, to keep stats of two-sided tests for inspection.
  ret <- cbind(EM_dir.p=ret$EM.p, MY_dir.p=ret$MY.p, ret)
  p.cols <- c("EM_dir.p", "MY_dir.p")
  ret <- modify_hitman_pvalues(tab=ret, overall.sign = ey.sign, p.cols=p.cols)
  
  EMY.p <- apply(ret[,p.cols], MARGIN=1, FUN=function(v){
    max(v)^2
  })
  EMY.FDR <- stats::p.adjust(EMY.p, method="BH")
  ret <- cbind(EMY.p, EMY.FDR, ret)
  ret <- ret[order(ret$EMY.p),]
  return(ret)
}