#' Test correlation of each row of object to phenotype using moderated variance
#'
#' Test correlation of each row of object to phenotype. By default, it uses the model 
#' \code{design=model.matrix(~1+phenotype)} and tests 2nd coefficient. See examples in vignette.
#'
#' @param phenotype Numeric vector of sample characteristics (e.g. phenotypes or treatments). 
#' Should be same length as \code{ncol(object)}.
#' @param prefix Character string to add to beginning of column names. \code{NULL} does not add a prefix.
#' @param coef Column index or column name of the linear model to test, passed to \code{\link{eztoptab}}.
#' @param reorder.rows Logical, should rows be reordered by p-value?
#' @param reduce.df Number degrees of freedom to subtract from residual. This may be necessary if 
#' \code{\link[limma]{removeBatchEffect}} was previously applied to \code{object}. Must be <= \code{df.residual} 
#' returned by \code{\link[limma]{lmFit}}.
#' @param check.names Logical; should \code{names(phenotype)==rownames(object)} be checked?
#' @inheritParams limma_contrasts
#' @return Data frame.
#' @details Exactly one of \code{design} or \code{phenotype} must be non-null. If \code{design} is \code{NULL} and 
#' \code{phenotype} is given, design will be calculated as \code{model.matrix(~0+phenotype)}. See further details 
#' in \code{\link[limma]{lmFit}}.
#' 
#' When \code{moderated} is FALSE, an error is generated if \code{trend} is TRUE.
#' @seealso \code{\link[limma]{lmFit}}; \code{\link[limma]{eBayes}}; \code{\link[ezlimma]{ezcor}}
#' @export

limma_cor <- function(object, phenotype=NULL, design=NULL, prefix=NULL, weights=NA, trend=FALSE, block=NULL, correlation=NULL,
                      adjust.method="BH", coef=2, reorder.rows=TRUE, moderated=TRUE, reduce.df=0, check.names=TRUE,
                      cols=c("AveExpr", "P.Value", "adj.P.Val", "logFC")){
  
   stopifnot(all(is.na(weights)) || is.null(weights) || dim(weights)==dim(object) || length(weights)==nrow(object) || 
             length(weights)==ncol(object), is.numeric(reduce.df), reduce.df >= 0, is.null(phenotype)!=is.null(design), 
             moderated || !trend)
  if (!is.null(phenotype)){
    stopifnot(length(phenotype)==ncol(object), limma::isNumeric(phenotype), !is.na(phenotype))
    if (check.names){
      stopifnot(names(phenotype)==colnames(object))
    }
    #if want to handle NAs in pheno, need to account for object, object$weights, and weights (as vector or matrix)
    design <- stats::model.matrix(~1+phenotype)
  } else {
    #if pheno is NULL, design was given
    # pheno could be in first column if not intercept
    stopifnot(is.numeric(design))
  }
  
  if (length(weights)!=1 || !is.na(weights)){
    if (!is.matrix(object) && !is.null(object$weights)){ warning("object$weights are being ignored") }
    fit <- limma::lmFit(object, design=design, block = block, correlation = correlation, weights=weights)
  } else {
    fit <- limma::lmFit(object, design=design, block = block, correlation = correlation)
  }
  
  if (reduce.df > 0){
    if (any(reduce.df >= fit$df.residual)){
      stop("reduce.df=", reduce.df, " >= df.residual=", min(fit$df.residual))
    }
    fit$df.residual <- fit$df.residual - reduce.df
  }
  
  fit2 <- ezebayes(fit, moderated=moderated, trend=trend)
  res.mat <- eztoptab(fit2, coef=coef, cols=cols, adjust.method=adjust.method)
  
  #change logFC to slope and get rid of FC
  colnames(res.mat) <- gsub("logFC", "slope", colnames(res.mat))
  res.mat <- res.mat[, setdiff(colnames(res.mat), "FC"), drop=FALSE]
  
  if (!reorder.rows){ res.mat <- res.mat[rownames(object),, drop=FALSE] }
  if (!is.null(prefix)){ colnames(res.mat) <- paste(prefix, colnames(res.mat), sep=".") }
  return(res.mat)
}