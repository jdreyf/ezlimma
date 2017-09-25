#'Combine p-values of a feature over multiple p-value columns of an object
#'
#'Combine p-values of a feature over multiple p-value columns of an object
#'
#'@param mat A matrix-like object with statistical columns, including some
#'  containing p-values.
#'@param pv.cols the column names or column indices with p-values. If
#'  \code{NULL}, the function searches for columns that end with \code{.p} or
#'  \code{.pval}.
#'@details This function uses the z-transform method to combine p-values, equivalently to using unweighted \code{method="z.transform"} in
#'\code{\link[survcomp]{combine.test}}.
#'@return p-value.
#'@seealso \code{\link[survcomp]{combine.test}}.
#'@export

#don't export
combine_pvalues <- function(mat, pv.cols=NULL){
  if (is.null(pv.cols)){
    pv.cols <- grep(paste0('\\.', '(p|pval)', '$'), colnames(mat), ignore.case=TRUE)
  }
  stopifnot(length(pv.cols)>0)
  comb.p <- apply(as.matrix(mat[,pv.cols]), 1, FUN=function(p){ 
    p.nona <- p[!is.na(p)]
    z <- qnorm(p.nona, lower.tail = FALSE)
    cp <- pnorm(sum(z)/sqrt(length(p.nona)), lower.tail = FALSE)
    return(cp)
  })
  return(comb.p)
}