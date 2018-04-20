#' F-test of each row of object using design matrix
#'
#' F-test of each row of object using design matrix for selected coefficients.
#' 
#' @param object A matrix-like data object containing log-ratios or 
#'  log-expression values for a series of samples, with rows corresponding to 
#'  genes and columns to samples.
#' @param design the design matrix of the experiment, with rows corresponding to 
#'  samples and columns to coefficients to be estimated. Can be used to provide 
#'  covariates.
#' @param coef Vector of 2 or more coefficients to include in F-test. These should be in \code{colnames(design)}. If your design matrix has an intercept, you likely want to exclude it.
#' @param prefix character string to add to beginning of column names.
#' @param weights non-negative observation weights. Can be a numeric matrix of 
#'  individual weights, of same size as the object expression matrix, or a 
#'  numeric vector of array weights with length equal to \code{ncol} of the 
#'  expression matrix, or a numeric vector of gene weights with length equal to 
#'  \code{nrow} of the expression matrix.
#' @param trend logical, should an intensity-trend be allowed for the prior 
#'  variance? Default is that the prior variance is constant.
#' @param block vector or factor specifying a blocking variable on the arrays. 
#'   Has length equal to the number of arrays.
#' @param correlation the inter-duplicate or inter-technical replicate 
#'   correlation.
#' @param adjust.method method used to adjust the p-values for multiple testing.
#' @param cols columns of \code{topTable} output the user would like in the 
#'  result. Some column names, such as \code{adj.P.Val} are changed. If \code{logFC}
#'  is specified, \code{FC} will also be given.
#' @param reorder.rows logical, should rows be reordered by F-statistic from 
#'  \code{\link[limma]{toptable}} or be left in the same order as 
#'  \code{object}?
#' @return Dataframe.
#' @seealso \code{\link[limma]{lmFit}} and \code{\link[limma]{eBayes}}.
#' @examples 
#' object <- matrix(rnorm(60), nrow=10, ncol=6)
#' grp <- rep(c("A", "B", "ctrl"), each=2)
#' covar <- cbind(covar1=1:6, covar2=rnorm(6))
#' design <- model.matrix(~0+grp+covar)
#' lf <- limmaF(object=object, design=design, coef=colnames(design)[1:3])
#' @export

limmaF <- function(object, design=NULL, coef=colnames(design), prefix='', weights=NULL, trend=FALSE, block = NULL, 
                   correlation = NULL, adjust.method='BH', reorder.rows=TRUE, cols=c('F', 'P.Value')){
  
  stopifnot(dim(weights)==dim(object)|length(weights)==nrow(object)|length(weights)==ncol(object),
            length(coef)>1, coef %in% colnames(design) | coef %in% 1:ncol(design))
  
  int <- grep("intercept", coef, ignore.case = TRUE, value = TRUE)
  if (length(int)>0) message("You included the column ", int, " in 'coefs', which you may want to exclude.")
  
  if (!missing(weights)){
    if (!is.matrix(object) && !is.null(object$weights)){ warning('object$weights are being ignored') }
    fit <- limma::lmFit(object, design=design, block = block, correlation = correlation, weights=weights)
  } else {
    fit <- limma::lmFit(object, design=design, block = block, correlation = correlation)
  }
  fit2 <- limma::eBayes(fit, trend=trend)
  res.mat <- eztoptab(fit2, coef=coef, cols=cols, adjust.method=adjust.method)
  
  if (!reorder.rows){ res.mat <- res.mat[rownames(object),] }
  if (prefix!=''){ colnames(res.mat) <- paste(prefix, colnames(res.mat), sep='.') }
  return(res.mat)
}