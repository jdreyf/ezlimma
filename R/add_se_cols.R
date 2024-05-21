#' Add column to toptab with standard errors
#' 
#' Add column to toptab with standard errors.
#' 
#' @param toptab Table with genes as rows and statistics as columns. Must include a column with t-statistics
#' e.g. named \code{"t"} or ends with \code{".t"} & a column of log fold changes e.g. named \code{"logFC"} or ends with \code{".logFC"}.
#' Columns with t-statistics and those with logFC's must have same prefixes.
#' @param lfc.suffix Suffix of logFC column(s).
#' @param t.suffix Suffix of t-statistic column(s).
#' @importFrom rlang :=

add_se_cols <- function(toptab, lfc.suffix="logFC", t.suffix="t"){
  lfc.colnms <- grep_cols(toptab, stat.cols = lfc.suffix)
  t.colnms <- grep_cols(toptab, stat.cols = t.suffix)
  stopifnot(length(lfc.colnms) == length(t.colnms))
  
  # prefix="" for t.colnms="t"
  prefix.lfc <- sub(pattern=paste0("(\\.|^)(", lfc.suffix, ")$"), replacement = "", lfc.colnms)
  prefix.t <- sub(pattern=paste0("(\\.|^)(", t.suffix, ")$"), replacement = "", t.colnms)
  stopifnot(sort(prefix.lfc) == sort(prefix.t))
  prefix.v <- prefix.lfc
  sep <- ifelse(prefix.lfc[1] == "", yes="", no=".")
  
  se.mat <- apply(as.matrix(prefix.v), MARGIN = 1, FUN=function(xx){
    lfc.colnm <- paste0(xx, sep, lfc.suffix)
    t.colnm <- paste0(xx, sep, t.suffix)
    toptab[, lfc.colnm]/toptab[, t.colnm]
  })
  colnames(se.mat) <- gsub(pattern = paste0(lfc.suffix, "$"), "SE", lfc.colnms)
  
  toptab2 <- toptab
  for (se.col in 1:ncol(se.mat)){
    se.colnm <- colnames(se.mat)[se.col]
    lfc.colnm <- lfc.colnms[se.col]
    toptab2 <- tibble::add_column(toptab2, rlang::`!!`(se.colnm) := se.mat[, se.colnm], .after = lfc.colnm)
  }
  toptab2
}