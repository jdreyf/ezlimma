#' Combine p-values of a feature over multiple p-value columns of an object
#' 
#' Combine p-values of a feature over multiple p-value columns of an object. If \code{alternative!="two.sided"},
#' uses the direction of \code{stat.cols} to transform \code{p.cols} into one-sided p-values; this assumes that
#' these p-values were originally calculated from two-sided tests, as is done in \pkg{limma}.
#' 
#' @param tab Matrix-like object with statistical columns, some containing p-values.
#' @param p.cols Indices or \code{\link{regexp}} with column names or column names suffix of numeric p-value columns.
#' @param stat.cols Indices or \code{\link{regexp}} with column names or column names suffix with numeric signed 
#' statistics, or with \code{"Up", "Down"} values. These should match \code{p.cols}.
#' @param only.p Logical; should only combined p-values be returned? If not, returns matrix with z-scores and FDRs also.
#' @param alternative Direction of change: \code{"two.sided"}; \code{"greater"} or \code{"less"}, or their synonyms  
#' \code{"Up"} or \code{"Down"}.
#' @details Z-transform method is used to combine p-values across rows, equivalently to using unweighted 
#' \code{method="z.transform"} in \code{survcomp::combine.test}.
#' 
#' \code{stat.cols} are ignored if \code{alternative} is \code{"two.sided"}.
#' @return Vector of p-values.
#' @examples
#'  tab <- data.frame(foo.p=(1:9)/9, bar.p=(9:1)/9)
#'  combine_pvalues(tab)
#' @export

combine_pvalues <- function(tab, p.cols="p|PValue", stat.cols="logFC|slope|cor|rho|Direction", only.p=TRUE,
                            alternative=c("two.sided", "greater", "less", "Up", "Down")){
  
  stopifnot(ncol(tab) >= 1, nrow(tab) >= 1, !is.null(colnames(tab)))
  alternative <- match.arg(alternative)
  p.colnms <- grep_cols(tab, p.cols=p.cols)
  tab.p <- data.matrix(tab[, p.colnms])
  if (any(tab.p == 0, na.rm = TRUE)){
    small.p <- max(10**(-15), min(tab.p[tab.p > 0])/2)
    warning("p-values should not be zero; these have been converted to ", small.p, ".", call. = FALSE)
    wh <- which(tab.p == 0)
    tab.p[wh] <- small.p
  }
  
  if (alternative != "two.sided"){
    stat.colnms <- grep_cols(tab, stat.cols=stat.cols)
    if (length(stat.colnms) != length(p.colnms)) stop("Lengths of p columns and stat columns must match.")
    for (col.ind in seq_along(length(p.colnms))){
      tab.p[, col.ind] <- two2one_tailed(tab, stat.cols=stat.colnms[col.ind], p.cols=p.colnms[col.ind], 
                                         alternative=alternative)
    }
    if (any(rowMeans(is.na(tab.p)) == 1)){
      stop("All rows of p-values, after accounting for stats, must not be all NA.")
    }
  }
  tab.z <- apply(as.matrix(tab.p), MARGIN=2, stats::qnorm, lower.tail = FALSE)
  # account for NAs
  combz.v <- apply(as.matrix(tab.z), MARGIN=1, FUN=function(zv){
    zv.nona <- zv[!is.na(zv)]
    sum(zv)/sqrt(length(zv))
  })
  combp.v <- stats::pnorm(combz.v, lower.tail = FALSE)
  if (!only.p){
    return(cbind(z=combz.v, p=combp.v, FDR=stats::p.adjust(combp.v, method = "BH")))
  } else {
    return(combp.v)
  }
}