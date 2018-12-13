#' Test (conditional) dependence of two variables
#'
#' Test (conditional) dependence of \code{Y} with rows of \code{object} via \code{\link[ezlimma]{limmaF}} or 
#' \code{\link[ezlimma]{limma_cor}}.
#' 
#' @param object A matrix-like data object with rows corresponding to features and columns to samples.
#' @param Y A numeric vector or matrix.
#' @param covariates A matrix-like data object with rows corresponding to samples and columns to covariates.
#' @param prefix Character string to add to beginning of column names.
#' @return Dataframe with \code{t-statistic} or \code{F-statistic} and \code{p-value} per row of \code{object}.
#' @export

limma_dep <- function(object, Y, covariates=NULL, prefix=NULL){
  stopifnot(ncol(object)==nrow(as.matrix(Y)), !is.null(covariates)||all(nrow(covariates)==nrow(as.matrix(Y))),
            is.numeric(Y))
  
  if (ncol(as.matrix(Y))==1){
    stopifnot(colnames(object)==names(Y))
  } else {
    stopifnot(colnames(object)==rownames(Y))
  }
  
  p.col <- "p"
  
  if (any(apply(X=as.matrix(Y), MARGIN=2, FUN=stats::var, na.rm=TRUE) == 0)){
    stop("Y treated as numeric, but has one or more columns with no variance.")
  }

  if (ncol(as.matrix(Y)) == 1){
    if (is.null(covariates)){
      tt <- limma_cor(object=object, phenotype = Y, cols=c('t', 'P.Value'))
    } else {
      dat <- data.frame(Y, covariates)
      design <- stats::model.matrix(~., data=dat)
      #design includes an intercept, so coefficient of Y is tested
      tt <- limma_cor(object = object, design = design, cols=c('t', 'P.Value'))
    }
  } else {
    coef <- 1:ncol(Y)
    design <- cbind(Y, covariates)
    tt <- limmaF(object, design = design, coef=coef, cols=c('F', 'P.Value'))
  }
  
  if (!is.null(prefix)){ colnames(tt) <- paste(prefix, colnames(tt), sep='.') }
  return(tt)
}