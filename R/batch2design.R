#' Transform factor vector into design matrix
#' 
#' Transform factor vector into design matrix, based on code from \code{\link[limma]{removeBatchEffect}}.
#' 
#' @param batch Vector that is a factor or can be coerced to one.
#' @return Numeric matrix.
#' @export

batch2design <- function(batch){
  stopifnot(length(unique(batch)) > 1)
  
  batch <- as.factor(batch)
  contrasts(batch) <- contr.sum(levels(batch))
  batch <- model.matrix(~batch)[, -1, drop = FALSE]
  return(batch)
}