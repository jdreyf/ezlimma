#' Wrapper for limma topTable function
#' 
#' Wrapper for \pkg{limma} \code{\link[limma]{topTable}} function that subsets and changes colnames
#' 
#' @param fit Object of class \code{MArrayLM} as produced by \code{\link[limma]{lmFit}} then \code{\link[limma]{eBayes}}.
#' @param coef Column index or column name specifying which coefficient or contrast of the linear model to test. 
#' If \code{length(coef)>1}, an F-test will be performed, and \code{logFC} & \code{FC} will not be returned.
#' @inheritParams limma_contrasts
#' @inheritParams limma_cor
#' @return Data frame.
#' @details If \code{length(coef)>=2}, \code{"z"} will not be returned.
#' @seealso \code{\link[limma]{topTable}}.

# sort by p
# assume that if "logFC" in cols, then want "FC"
# limma_contrasts tests one coef at a time
eztoptab <- function(fit, cols=c("P.Value", "adj.P.Val", "logFC"), adjust.method="BH", prefix=NULL, coef=NULL){
  stopifnot(length(cols)>=1, cols %in% c("CI.L", "CI.R", "AveExpr", "z", "t", "F", "P.Value", "adj.P.Val", "B", "logFC", "SE"))
  confint.bool <- ifelse(any(c("CI.L", "CI.R") %in% cols), yes = TRUE, no = FALSE)
  
  if (!is.null(coef) && length(coef)>=2){
    # topTable tests all coefficients using F-test if at least 2 of them
    tt <- limma::topTable(fit, number=Inf, adjust.method=adjust.method, coef=coef, confint = confint.bool)
    cols <- setdiff(cols, c("logFC", "t"))
  } else {
    tt <- limma::topTable(fit, number=Inf, sort.by="P", adjust.method=adjust.method, coef=coef, confint = confint.bool)
    if ("z" %in% cols){
      tt <- add_zcols(tt, p.suffix = "P.Value")
    }
    if ("SE" %in% cols){
      tt <- add_se_cols(tt, lfc.suffix="logFC", t.suffix="t")
    }
  }
  
  # FC
  if ("logFC" %in% cols){
    tt$FC <- logfc2fc(tt$logFC)
    cols <- c(cols, "FC")
  }
  
  tt <- tt[, cols, drop=FALSE]
  colnames(tt) <- sub("P.Value", "p", colnames(tt))
  
  # p.adjust says fdr is alias for BH
  if (adjust.method %in% c("BH", "fdr")){
    colnames(tt) <- sub("adj\\.P\\.Val", "FDR", colnames(tt))
  } else {
    colnames(tt) <- sub("adj\\.P\\.Val", adjust.method, colnames(tt))
  } 
  if (!is.null(prefix)){ colnames(tt) <- paste(prefix, colnames(tt), sep=".") }
  return(tt)
}