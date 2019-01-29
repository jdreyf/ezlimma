#' Test (conditional) dependence of two variables
#'
#' Test (conditional) dependence of \code{Y} with rows of \code{object} via \code{\link[ezlimma]{limmaF}}. Let a row of
#' \code{object} be \eqn{o}, then we test the coefficient(s) of \code{Y} in the model \eqn{o = 1 + Y + covariates}. All
#' inputs should be numeric, unlike \code{\link[ezlimma]{limma_contrasts}}.
#' 
#' @param Y A numeric vector or matrix.
#' @param covariates A numeric matrix-like data object with rows corresponding to samples and columns to covariates
#' to be adjusted for.
#' @inheritParams limma_contrasts
#' @inheritParams limma_cor
#' @return Data frame with \code{t-statistic} or \code{F-statistic} and \code{p-value} per row of \code{object}.
#' @export

limma_dep <- function(object, Y, covariates=NULL, prefix=NULL){
  stopifnot(ncol(object)==nrow(as.matrix(Y)), !is.null(covariates)||all(nrow(covariates)==nrow(as.matrix(Y))),
            limma::isNumeric(Y))
  
  p.col <- "p"
  
  if (ncol(as.matrix(Y))==1){
    stopifnot(colnames(object)==names(Y))
  } else {
    stopifnot(colnames(object)==rownames(Y))
  }
  
  if (any(apply(X=as.matrix(Y), MARGIN=2, FUN=stats::var, na.rm=TRUE) == 0)){
    stop("Y treated as numeric, but has one or more columns with no variance.")
  }

  #dat must be data.frame for model.matrix
  dat <- data.frame(cbind(Y, covariates))
  # include intercept in the design matrix
  design <- stats::model.matrix(~., data=dat)
  # exclude intercept in coef
  # if ncol(Y) == 1, then coef = 2
  coef <- 1:ncol(as.matrix(Y)) + 1
  #limmaF == limma_cor Y is vector, will return 't' instead of 'F'
  tt <- limmaF(object=object, design = design, coef=coef, prefix=prefix)
  return(tt)
}