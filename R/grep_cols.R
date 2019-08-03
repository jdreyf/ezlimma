#' Find p-value or stat columns in table
#' 
#' Find p-value or stat columns in table and check their validity.
#' 
#' @inheritParams combine_pvalues
#' @details Exactly one of \code{p.cols} or \code{stat.cols} should not be \code{NULL}.
#' @return Vector of column names.

# cannot set defaults, since xor
grep_cols <- function(tab, p.cols=NULL, stat.cols=NULL){
  stopifnot(xor(is.null(p.cols), is.null(stat.cols)), !all(is.numeric(p.cols)) || all(p.cols %in% 1:ncol(tab)),  
            !all(is.numeric(stat.cols)) || all(stat.cols %in% 1:ncol(tab)), !is.null(colnames(tab)))
  
  if (!is.null(p.cols)){
    if (is.numeric(p.cols)){
      colnms <- colnames(tab)[p.cols]
    } else {
      colnms <- grep(pattern=paste0("(\\.|^)(", p.cols, ")$"), colnames(tab), ignore.case=TRUE, value = TRUE)
    }
    if (length(colnms) == 0) stop("Cannot find p cols: '", p.cols, "'.", call. = FALSE)
    if (any(tab[, colnms] < 0 | tab[, colnms] > 1, na.rm = TRUE)) stop("p-values must be 0 <= p <= 1.")
    if (any(rowMeans(as.matrix(is.na(tab[, colnms]))) == 1)) stop("Rows of p-value columns must not be all NA.")
  } else {
    if (is.numeric(stat.cols)){
      colnms <- colnames(tab)[stat.cols]
    } else {
      colnms <- grep(pattern=paste0("(\\.|^)(", stat.cols, ")$"), colnames(tab), ignore.case=TRUE, value = TRUE)
    }
    if (length(colnms) == 0) stop("Cannot find stat cols: '", stat.cols, "'.", call. = FALSE)
  }
  colnms
}