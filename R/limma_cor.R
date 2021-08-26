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
#' The defaults of arguments \code{ndups} and \code{spacing} are set to \code{NULL}, 
#' which allows these arguments to be overridden by the elements \code{object$printer$ndups} and
#' \code{object$printer$spacing}, respectively, if these exist. Whereas, if an element does not exist,
#' the corresponding argument is treated as being its default in \code{\link[limma]{lmFit}}, 
#' i.e. \code{ndups=1} or \code{spacing=1}. 
#' If either of these arguments are specified, they would override the 
#' respective element of \code{object$printer}, if the element existed.
#' 
#' When \code{moderated} is FALSE, an error is generated if \code{trend} is TRUE.
#' @seealso \code{\link[limma]{lmFit}}; \code{\link[limma]{eBayes}}; \code{\link[ezlimma]{ezcor}}
#' @export

limma_cor <- function(object, phenotype=NULL, design=NULL, prefix=NULL, weights=NA, trend=FALSE, ndups=NULL, spacing=NULL,
                      block=NULL, correlation=NULL, adjust.method="BH", coef=2, reorder.rows=TRUE, moderated=TRUE, 
                      reduce.df=0, check.names=TRUE, cols=c("AveExpr", "P.Value", "adj.P.Val", "logFC")){
  
  stopifnot(!is.null(dim(object)), !is.null(rownames(object)), !is.null(colnames(object)), ncol(object)>1,
            length(weights)!=1 || is.na(weights), length(weights)<=1 || 
              (is.numeric(weights) && all(weights>=0) && !all(is.na(weights))), 
            length(weights)<=1 || all(dim(weights)==dim(object)) || 
              length(weights)==nrow(object) || length(weights)==ncol(object),
            is.null(ndups) || (is.numeric(ndups) && ndups >= 1), 
            is.null(spacing) || (is.numeric(spacing) && spacing>=1), 
            is.null(correlation) || (is.numeric(correlation) && abs(correlation)<=1),
            is.null(ndups) || ndups < 2 || is.null(block),  moderated || !trend,
            is.numeric(reduce.df), reduce.df >= 0, is.null(phenotype)!=is.null(design))
  
  if ((!is.null(block) || (!is.null(ndups) && ndups > 1)) && is.null(correlation))
    stop("!is.null(block) or ndups>1, so correlation must not be NULL.")
  
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
  
  args.lst <- list(object=object, design=design, block = block, correlation = correlation)
  
  if (length(weights)!=1 || !is.na(weights)) args.lst <- c(args.lst, list(weights=weights))
  if (!is.null(ndups)) args.lst <- c(args.lst, ndups=ndups)
  if (!is.null(spacing)) args.lst <- c(args.lst, spacing=spacing)
  
  # use do.call() so that can accommodate ndups/spacing
  fit <- do.call(limma::lmFit, args = args.lst)
  
  # if ndups>1, lmFit returns y$genes=uniquegenelist(y$probes,ndups=ndups,spacing=spacing), but does not fix row names
  # assign row names before multitoptab, since it reorders rows 
  if (!is.null(ndups) && ndups >= 2){
    if (is.null(spacing)) spacing <- 1
    rownames(fit$coefficients) <- limma::uniquegenelist(rownames(object), ndups=ndups, spacing=spacing)
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
  
  if (!reorder.rows){ 
    # if ndups > 1, nrow(res.mat) < nrow(object)
    rnms <- rownames(object)[rownames(object) %in% rownames(res.mat)]
    res.mat <- res.mat[rnms,, drop=FALSE]
  }
  if (!is.null(prefix)){ colnames(res.mat) <- paste(prefix, colnames(res.mat), sep=".") }
  return(res.mat)
}