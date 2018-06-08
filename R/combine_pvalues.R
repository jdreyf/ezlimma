#' Combine p-values of a feature over multiple p-value columns of an object
#' 
#' Combine p-values of a feature over multiple p-value columns of an object
#' 
#' @param mat A matrix-like object with statistical columns, including some
#'   containing p-values. Must have \code{nrow(mat)>1} & \code{ncol(mat)>1}.
#' @param pv.cols the column names or column indices with p-values. If
#'   \code{NULL}, the function searches for columns that end with \code{.p} or \code{.pval}.
#' @details This function uses the z-transform method to combine p-values across rows, equivalently to using unweighted 
#' \code{method="z.transform"} in \code{survcomp::combine.test}.
#' @return A vector of p-values.
#' @examples
#'  tab <- data.frame(foo.p=(1:9)/9, bar.p=(9:1)/9)
#'  combine_pvalues(tab)
#' @export

#don't export
combine_pvalues <- function(mat, pv.cols=NULL){
  stopifnot(ncol(mat) > 1, nrow(mat) > 1, !is.null(pv.cols) || !is.null(colnames(mat)))
  #if pv.cols not given, grep for them at end of column names
  if (is.null(pv.cols)){
    pv.cols <- grep(pattern=paste0('(\\.|^)', '(p|pval)', '$'), colnames(mat), ignore.case=TRUE)
  }
  stopifnot(length(pv.cols)>0)
  comb.p <- apply(as.matrix(mat[,pv.cols]), MARGIN=1, FUN=function(p){ 
    p.nona <- p[!is.na(p)]
    z <- stats::qnorm(p.nona, lower.tail = FALSE)
    cp <- stats::pnorm(sum(z)/sqrt(length(p.nona)), lower.tail = FALSE)
    return(cp)
  })
  return(comb.p)
}