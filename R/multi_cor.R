#' Test correlation of each row of an object to each column of pheno.tab
#' 
#' Test correlation of each row of an object to each column of \code{pheno.tab} using 
#' one of Pearson's, Kendall's, or Spearman's correlation methods, or limma 
#' regression in \code{\link{limma_cor}}.
#' 
#' @param pheno.tab Matrix-like data object with columns as sample phenotypes, with \code{nrow(pheno.tab)==ncol(object)}.
#' @param method a character string indicating which association is to be used 
#'   for the test. One of \code{"pearson"}, \code{"spearman"}, \code{"kendall"}, 
#'   from \code{\link[stats]{cor.test}} or \code{"limma"} for \code{\link{limma_cor}}.
#' @param check_names Logical indicating if \code{names(pheno.tab)=rownames(object)} should be checked.
#' @param limma.cols if \code{method="limma"}, \code{cols} from \code{\link{limma_cor}} to include.
#' @inheritParams limma_contrasts
#' @inheritParams limma_cor
#' @return Data frame with several statistical columns corresponding to each phenotype and one row per feature.
#' @details Each column of \code{pheno.tab} is tested independently.
#' @export

multi_cor <- function(object, pheno.tab, method=c('pearson', 'spearman', 'kendall', 'limma'), reorder.rows=TRUE, 
                      prefix=NULL, adjust.method='BH', check_names=TRUE,
                     limma.cols=c('AveExpr', 'P.Value', 'adj.P.Val', 'logFC')){
  method <- match.arg(method)
  
  if (is.null(dim(pheno.tab))) stop("pheno.tab needs to have rows and columns.")
  stopifnot(ncol(object)==nrow(pheno.tab))
  if (check_names){
    stopifnot(rownames(pheno.tab)==colnames(object))
  }
  
  cor.mat <- NULL
  for (i in 1:ncol(pheno.tab)){
    prefix.tmp <- ifelse(!is.null(prefix), paste(prefix, colnames(pheno.tab)[i], sep='.'), colnames(pheno.tab)[i])
    if (method=='limma'){
      cor.tmp <- data.matrix(limma_cor(object, pheno.tab[,i], reorder.rows=FALSE, prefix=prefix.tmp, cols=limma.cols))
    } else {
      cor.tmp <- ezcor(object, pheno.tab[,i], method=method, reorder.rows=FALSE, prefix=prefix.tmp, adjust.method=adjust.method)
    }
    if (!is.null(cor.mat)){ stopifnot(rownames(cor.mat)==rownames(cor.tmp)) }
    cor.mat <- cbind(cor.mat, cor.tmp)
  }
  if (reorder.rows){ cor.mat <- cor.mat[order(combine_pvalues(cor.mat)),] }
  return(cor.mat)
}