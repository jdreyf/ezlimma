#' Test (conditional) dependence of two variables
#'
#' Test (conditional) dependence of \code{Y} with rows of \code{object}.
#' 
#' @param object A matrix-like data object with rows corresponding to features and columns to samples.
#' @param Y A vector. If not numeric, treated as a nominal group variable. To ensure that \code{object} matches \code{Y},
#' it is checked that \code{colnames(object)==names(Y)}. This check will pass if one of these is \code{NULL}.
#' @param covariates A matrix-like data object with rows corresponding to samples and columns to covariates. If \code{covariates}
#' has rownames, they should match \code{names(Y)}.
#' @param prefix Character string to add to beginning of column names.
#' @param verbose If \code{TRUE} reports messages.
#' @return Dataframe with \code{t-statistic} or \code{F-statistic} and \code{p-value} per row of \code{object}.
#' @export
#' @import stats

#transform to |z-score| using sqrt of chisq w/ df=1
limma_dep <- function(object, Y, covariates=NULL, prefix=NULL, verbose=FALSE){
  stopifnot(ncol(object)==length(Y), colnames(object)==names(Y), rownames(covariates)==names(Y),
            !is.null(covariates)||all(nrow(covariates)==length(Y)))
  p.col <- "p"
  if (is.numeric(Y)){
    if (verbose) message("'Y' treated as continuous numeric vector")
    if (stats::var(Y, na.rm=TRUE)==0) stop("'Y' treated as numeric, but has no variance.")
    if (is.null(covariates)){
      tt <- limma_cor(object=object, phenotype = Y, cols=c('t', 'P.Value'))
    } else {
      dat <- data.frame(Y, covariates)
      design <- stats::model.matrix(~., data=dat)
      tt <- limma_cor(object = object, design = design, cols=c('t', 'P.Value'))
    }
  } else {
    if (verbose) message("'Y' treated as an unordered factor")
    ngrps <- length(unique(Y))
    if (ngrps==1) stop("'Y' is not numeric, but has only one group.")
    if (any(is.na(Y))) stop("'Y' is not numeric, but has an NA.")
    if (!is.null(covariates)){
      design <- stats::model.matrix(~1+Y+covariates)
      covar.ncol <- ncol(as.matrix(covariates))
      coef <- 2:(ncol(design)-covar.ncol)
    } else {
      design <- stats::model.matrix(~1+Y)
      coef <- 2:ncol(design)
    }

    if (length(coef)==1){
      tt <- limmaF(object, design = design, coef=coef, cols=c('t', 'P.Value'))
    } else {
      tt <- limmaF(object, design = design, coef=coef, cols=c('F', 'P.Value'))
    }
  }#end else not numeric
  if (!is.null(prefix)){ colnames(tt) <- paste(prefix, colnames(tt), sep='.') }
  return(tt)
}