# Test association of gene sets to multiple phenotypes using rotation testing with output
#' to Excel
#'
#' Test association of gene sets to multiple phenotypes while adjusting for \code{covariates} using rotation testing with
#' \code{\link[limma]{roast}} using \pkg{limma} functions \code{mroast} or \code{fry}.
#' It returns a data frame with statistics per gene set across all phenotypes, and writes this to an
#' Excel file. The Excel file links to CSV files, which contain statistics per
#' gene set across all phenotypes.
#'
#' @param covariates Numeric vector or matrix of covariates to include in constructed \code{design} matrix.
#' @inheritParams limma_contrasts
#' @inheritParams limma_cor
#' @inheritParams roast_contrasts
#' @inheritParams multi_cor
#' @return Data frame of gene set statistics.
#' @inherit roast_contrasts details
#' @seealso \code{\link[ezlimma]{roast_cor}}.
#' @export

roast_multi_cor <- function(object, G, pheno.tab, feat.tab, name=NA, covariates=NULL, 
                    fun=c("fry", "mroast"), set.statistic = "mean", weights = NA, gene.weights=NULL, 
                    trend = FALSE, block = NULL, correlation = NULL, prefix=NULL, adjust.method = "BH", 
                    min.nfeats=3, max.nfeats=1000, alternative=c("two.sided", "less", "greater"), 
                    nrot=999, check.names=TRUE, pwy.nchar=199, seed=0){
  
  stopifnot(!is.null(dim(object)), !is.null(rownames(object)), !is.null(colnames(object)), ncol(object) > 1,
            !is.null(colnames(pheno.tab)), nrow(pheno.tab) == ncol(object), limma::isNumeric(pheno.tab), colSums(!is.na(pheno.tab)) > 2,
            length(weights)!=1 || is.na(weights), length(weights)<=1 || 
              (is.numeric(weights) && all(weights>=0) && !all(is.na(weights))), 
            length(weights)<=1 || all(dim(weights)==dim(object)) || 
              length(weights)==nrow(object) || length(weights)==ncol(object),
            is.null(gene.weights) || length(gene.weights)==nrow(object),
            is.na(name) || rownames(object) %in% rownames(feat.tab))
  if (!is.null(covariates)) stopifnot(limma::isNumeric(covariates), !is.na(covariates))
  if (check.names) stopifnot(rownames(pheno.tab)==colnames(object))
  
  # get G index
  index <- g_index(G=G, object=object, min.nfeats=min.nfeats, max.nfeats=max.nfeats)
  
  # would be faster to pre-specify size, but I don't think it matters enough here
  rc.mat <- NULL
  weights.tmp <- weights
  for (ind in 1:ncol(pheno.tab)){
    prefix.tmp <- ifelse(!is.null(prefix), paste(prefix, colnames(pheno.tab)[ind], sep="."), colnames(pheno.tab)[ind])
    # need to account for weights & block with NAs in pheno
    # na.omit() missing phenotypes in pheno.tab
    ph.idx <- 1:nrow(pheno.tab)
    if (any(is.na(pheno.tab[,ind]))){
      n.na <- sum(is.na(pheno.tab[,ind]))
      message(n.na, " NAs in ", colnames(pheno.tab)[ind], " removed")
      ph.idx <- which(!is.na(pheno.tab[,ind]))
      # block is on samples, so should be modified in case of NAs, but NULL[ph.idx] is NULL
      if (!is.null(block)) block <- block[ph.idx]
      # len > 1 excludes NA and NULL
      if (length(weights.tmp) > 1){
        if (length(weights.tmp)==ncol(object)){ 
          weights.tmp <- weights.tmp[ph.idx]
        } else {
          weights.tmp <- weights.tmp[, ph.idx, drop=FALSE]
        }
      }
    }
    
    if (is.null(covariates)){
      des.tmp <- stats::model.matrix.lm(~1+pheno.tab[, ind], na.action = stats::na.omit)
    } else {
      # model.matrix.lm, but not model.matrix, respects na.action
      # https://stackoverflow.com/questions/5616210/model-matrix-with-na-action-null
      des.tmp <- stats::model.matrix.lm(~1+pheno.tab[, ind]+covariates, na.action = stats::na.omit)
    }
    
    rc.tmp <- roast_cor(object=object[, ph.idx], G=G, feat.tab=NULL, name=NA, phenotype = NULL, design = des.tmp, prefix=prefix.tmp, weights = weights.tmp,
                        fun=fun, set.statistic = set.statistic, gene.weights=gene.weights, 
                        trend = trend, block = block, correlation = correlation, adjust.method = adjust.method, 
                        min.nfeats=min.nfeats, max.nfeats=max.nfeats, alternative=alternative, 
                        nrot=nrot, check.names=check.names, pwy.nchar=pwy.nchar, seed=seed)
    
    if (is.null(rc.mat)){ 
      rc.mat <- rc.tmp
      rc.rnms <- rownames(rc.mat)
    } else {
      stopifnot(sort(rownames(rc.mat)) == sort(rownames(rc.tmp))) 
      # don't want NGenes from additional rc.tmp
      rc.mat <- cbind(rc.mat, rc.tmp[rc.rnms, -1])
    }
  } # end for ind
  
  rc.mat <- rc.mat[order(combine_pvalues(rc.mat)), ]
  
  res.xl <- df_signif(rc.mat, digits = 8)
  # write xlsx file with links
  if (!is.na(name)){
    nm <- paste(name, fun, sep="_")
    write_linked_xl(pwy.tab=res.xl, feat.lst=index, feat.tab=feat.tab, name=nm, pwy.nchar=pwy.nchar)
  }
  return(rc.mat)
}
