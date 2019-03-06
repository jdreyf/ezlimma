#' Combine p-values of a feature over multiple p-value columns of an object
#' 
#' Combine p-values of a feature over multiple p-value columns of an object.
#' 
#' @param tab Matrix-like object with statistical columns, including some containing p-values. Must have 
#' \code{nrow(tab)>1} & \code{ncol(tab)>1}.
#' @param p.cols Column names or column indices with p-values. If \code{NULL}, the function searches for
#' columns that end with \code{.p} or \code{.pval}.
#' @details Z-transform method is used to combine p-values across rows, equivalently to using unweighted 
#' \code{method="z.transform"} in \code{survcomp::combine.test}.
#' @return Vector of p-values.
#' @examples
#'  tab <- data.frame(foo.p=(1:9)/9, bar.p=(9:1)/9)
#'  combine_pvalues(tab)
#' @export

combine_pvalues <- function(tab, p.cols=NULL){
  stopifnot(ncol(tab) >= 1, nrow(tab) >= 1, !is.null(p.cols) || !is.null(colnames(tab)))
  # if p.cols not given, grep for them at end of column names
  if (is.null(p.cols)){
    p.cols <- grep(pattern=paste0("(\\.|^)", "(p|pval)", "$"), colnames(tab), ignore.case=TRUE)
  }
  stopifnot(length(p.cols)>0)
  comb.p <- apply(as.matrix(tab[,p.cols]), MARGIN=1, FUN=function(p){ 
    p.nona <- p[!is.na(p)]
    z <- stats::qnorm(p.nona, lower.tail = FALSE)
    cp <- stats::pnorm(sum(z)/sqrt(length(p.nona)), lower.tail = FALSE)
    return(cp)
  })
  return(comb.p)
}