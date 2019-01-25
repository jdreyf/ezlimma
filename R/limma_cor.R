#' Test correlation of each row of object to phenotype using moderated variance
#'
#' Test correlation of each row of object to phenotype. By default, it uses the model 
#' \code{design=model.matrix(~1+phenotype)} and tests 2nd coefficient.
#'
#' @param object A matrix-like data object containing log-ratios or 
#'  log-expression values for a series of samples, with rows corresponding to 
#'  genes and columns to samples.
#' @param phenotype Vector of phenotypes of the samples. Should be same length as
#'  \code{ncol(object)}. If the vector is named, names should match 
#'  \code{colnames(object)}.
#' @param design the design matrix of the experiment, with rows corresponding to 
#'  samples and columns to coefficients to be estimated. Can be used to provide 
#'  covariates.
#' @param prefix character string to add to beginning of column names.
#' @param weights non-negative observation weights. Can be a numeric matrix of 
#'  individual weights, of same size as the object expression matrix, or a 
#'  numeric vector of array weights with length equal to \code{ncol} of the 
#'  expression matrix, or a numeric vector of gene weights with length equal to 
#'  \code{nrow} of the expression matrix. Set to \code{NULL} to ignore \code{object$weights}. \code{weights=NA} 
#'   (with length one) doesn't pass weights to \code{limma}.
#' @param trend logical, should an intensity-trend be allowed for the prior 
#'  variance? Default is that the prior variance is constant.
#' @param adjust.method method used to adjust the p-values for multiple testing.
#' @param coef Integer coefficient of the linear model to test; passed to \code{\link[ezlimma]{eztoptab}}.
#' @param reorder.rows logical, should rows be reordered by F-statistic from 
#'  \code{\link[limma]{toptable}} or be left in the same order as 
#'  \code{object}?
#' @param reduce.df Number of degrees of freedom to subtract from residual. This may be necessary if
#' \code{\link[limma]{removeBatchEffect}} was previously applied to \code{object}. Must be <= df.residual returned by
#' \code{\link[limma]{lmFit}}.
#' @param check_names Logical indicating if \code{names(phenotype)=rownames(object)} should be checked.
#' @param cols columns of \code{topTable} output the user would like in the 
#'  result. Some column names, such as \code{adj.P.Val} are changed. If \code{logFC}
#'  is specified, \code{FC} will also be given.
#' @return Dataframe.
#' @details Exactly one of \code{design} or \code{phenotype} must be non-null. If \code{design} is \code{NULL} and \code{phenotype} 
#' is given, design will be calculated as \code{model.matrix(~0+phenotype)}. See further details in \code{\link[limma]{lmFit}}.
#' @seealso \code{\link[limma]{lmFit}} and \code{\link[limma]{eBayes}}.
#' @export

limma_cor <- function(object, phenotype=NULL, design=NULL, prefix=NULL, weights=NA, trend=FALSE, adjust.method='BH', coef=2,
                      reorder.rows=TRUE, reduce.df=0, check_names=TRUE,
                      cols=c('AveExpr', 'P.Value', 'adj.P.Val', 'logFC')){
   stopifnot(dim(weights)==dim(object)|length(weights)==nrow(object)|length(weights)==ncol(object), is.numeric(reduce.df), 
             reduce.df >= 0, is.null(phenotype)!=is.null(design))
  
  if (!is.null(phenotype)){
    stopifnot(length(phenotype)==ncol(object), limma::isNumeric(phenotype), !is.na(phenotype))
    if (check_names){
      stopifnot(names(phenotype)==colnames(object))
    }
    #if want to handle NAs in pheno, need to account for object, object$weights, and weights (as vector or matrix)
    design <- stats::model.matrix(~1+phenotype)
  } else {
    #if pheno is NULL, design was given
    stopifnot(is.numeric(design[,2]))
  }
  
  if (length(weights)!=1 || !is.na(weights)){
    if (!is.matrix(object) && !is.null(object$weights)){ warning('object$weights are being ignored') }
    fit <- limma::lmFit(object, design, weights=weights)
  } else {
    fit <- limma::lmFit(object, design)
  }
  
  if (reduce.df > 0){
    if (any(reduce.df >= fit$df.residual)){
      stop("reduce.df=", reduce.df, " >= df.residual=", min(fit$df.residual))
    }
    fit$df.residual <- fit$df.residual - reduce.df
  }
  
  fit2 <- limma::eBayes(fit, trend=trend)
  res.mat <- eztoptab(fit2, coef=coef, cols=cols, adjust.method=adjust.method)
  
  #change logFC to coeff and get rid of FC
  colnames(res.mat) <- gsub('logFC', 'slope', colnames(res.mat))
  res.mat <- res.mat[,setdiff(colnames(res.mat), 'FC')]
  
  if (!reorder.rows){ res.mat <- res.mat[rownames(object),] }
  if (!is.null(prefix)){ colnames(res.mat) <- paste(prefix, colnames(res.mat), sep='.') }
  return(res.mat)
}