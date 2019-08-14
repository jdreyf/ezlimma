#' Add column to toptab with z-stats
#' 
#' Add column to toptab with z-stats.
#' 
#' @param toptab Table with genes as rows and statistics as columns. Must include a column with t-statistics that
#' is \code{"t"} or ends with \code{".t"} & a column of p-values is \code{"p"} or ends with \code{".p"}.
#' Columns with t-statistics and those with p-values must have same prefixes.
#' @param stat.suffix SUffix of numeric signed statistics column(s).
#' @param p.suffix SUffix of numeric p-value column(s).
#' @importFrom rlang :=

add_zcols <- function(toptab, stat.suffix="t", p.suffix="p"){
  stat.colnms <- grep_cols(toptab, stat.cols = stat.suffix)
  p.colnms <- grep_cols(toptab, p.cols = p.suffix)
  stopifnot(length(p.colnms) == length(stat.colnms))
  
  # prefix="" for t.colnms="t"
  prefix.stat <- sub(pattern=paste0("(\\.|^)(", stat.suffix, ")$"), replacement = "", stat.colnms)
  prefix.p <- sub(pattern=paste0("(\\.|^)(", p.suffix, ")$"), replacement = "", p.colnms)
  stopifnot(sort(prefix.p) == sort(prefix.stat))
  prefix.v <- prefix.stat
  sep <- ifelse(prefix.stat[1] == "", yes="", no=".")
  
  # moderated t-statistic df = fit$df.residual+fit$df.prior (https://support.bioconductor.org/p/6124/)
  # this is df.total
  z.mat <- apply(as.matrix(prefix.v), MARGIN = 1, FUN=function(xx){
    stat.colnm <- paste0(xx, sep, stat.suffix)
    p.colnm <- paste0(xx, sep, p.suffix)
    stats::qnorm(toptab[, p.colnm]/2, lower.tail = FALSE) * sign(toptab[, stat.colnm])
  })
  colnames(z.mat) <- gsub(pattern = paste0(stat.suffix, "$"), "z", stat.colnms)
  
  toptab2 <- toptab
  for (z.col in 1:ncol(z.mat)){
    z.colnm <- colnames(z.mat)[z.col]
    stat.colnm <- stat.colnms[z.col]
    toptab2 <- tibble::add_column(toptab2, rlang::`!!`(z.colnm) := z.mat[, z.colnm], .before = stat.colnm)
  }
  toptab2
}