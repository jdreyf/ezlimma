#' Transform factor vector into design matrix
#' 
#' Transform factor vector into design matrix, based on code from \code{\link[limma]{removeBatchEffect}}.
#' 
#' @param batch Vector that is a factor or can be coerced to one.
#' @return Numeric matrix.

batch2design <- function(batch){
  stopifnot(length(unique(batch)) > 1)
  batch <- as.factor(batch)
  stats::contrasts(batch) <- stats::contr.sum(levels(batch))
  batch <- stats::model.matrix(~batch)[, -1, drop = FALSE]
  return(batch)
}