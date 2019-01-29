#' F-test of each row of object using design matrix
#'
#' F-test of each row of object using design matrix for selected coefficients. If there's only one selected coefficient, 
#' a t-test is applied.
#' 
#' @inheritParams limma_contrasts
#' @inheritParams limma_cor
#' @return Data frame.
#' @seealso \code{\link[limma]{lmFit}} and \code{\link[limma]{eBayes}}.
#' @examples 
#' object <- matrix(rnorm(60), nrow=10, ncol=6)
#' grp <- rep(c("A", "B", "ctrl"), each=2)
#' covar <- cbind(covar1=1:6, covar2=rnorm(6))
#' design <- model.matrix(~0+grp+covar)
#' lf <- limmaF(object=object, design=design, coef=colnames(design)[1:3])
#' @export

limmaF <- function(object, design, coef=colnames(design), prefix=NULL, weights=NULL, trend=FALSE, block = NULL, 
                   correlation = NULL, adjust.method='BH', reorder.rows=TRUE, cols=c('F', 'P.Value')){
  
  stopifnot(dim(weights)==dim(object)|length(weights)==nrow(object)|length(weights)==ncol(object),
            coef %in% colnames(design) | coef %in% 1:ncol(design))
  #limma only retains "F" if length(coef)>1
  if (length(coef)==1 && "F" %in% cols) cols <- sub("^F$", "t", cols)
  
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
  if (!is.null(prefix)){ colnames(res.mat) <- paste(prefix, colnames(res.mat), sep='.') }
  return(res.mat)
}