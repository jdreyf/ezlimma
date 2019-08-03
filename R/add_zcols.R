#' Add column to toptab with z-stats
#' 
#' Add column to toptab with z-stats using \code{\link[limma]{zscore}}.
#' 
#' @param toptab Table with genes as rows and statistics as columns. Must include a column with t-statistics that
#' is \code{"t"}, or ends with \code{".t"}.
#' @param fit object of class \code{\link[limma]{marraylm}}.
#' @importFrom rlang :=

add_zcols <- function(toptab, fit){
  t.colnms <- grep_cols(toptab, stat.cols = "t")
  # tcol.ind <- which(colnames(toptab) %in% t.colnms)
  
  # moderated t-statistic df = fit$df.residual+fit$df.prior (https://support.bioconductor.org/p/6124/)
  # this is df.total
  z.mat <- apply(as.matrix(toptab[, t.colnms]), MARGIN = 2, FUN=function(t.v){
    limma::zscoreT(x=t.v, df=fit$df.total)
  })
  colnames(z.mat) <- gsub(pattern = "t$", "z", t.colnms)
  
  toptab2 <- toptab
  for (z.col in 1:ncol(z.mat)){
    z.colnm <- colnames(z.mat)[z.col]
    t.colnm <- t.colnms[z.col]
    toptab2 <- tibble::add_column(toptab2, rlang::`!!`(z.colnm) := z.mat[, z.colnm], .before = t.colnm)
  }
  toptab2
}