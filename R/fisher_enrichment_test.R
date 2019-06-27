#' Fisher Exact Test of gene set enrichment with output to Excel
#' 
#' Test gene set enrichment using \code{\link[stats]{fisher.test}}. It returns a data frame with statistics per gene set,
#' and writes this to an Excel file. The Excel file links to CSV files, which contain statistics per gene set.
#' 
#' @param gsets Named list of gene sets whose element is a vector of (differential) gene IDs corresponding to the rownames of \code{feat.tab}
#' @inheritParams roast_contrasts
#' @return Data frame of gene set statistics.
#' @details Pathway (i.e. gene set) names are altered to be valid filenames in Windows and Linux. Numeric columns are
#' rounded to 3 significant figures.
#' @export

fisher_enrichment_test <- function(gsets, G, feat.tab, name=NA, adjust.method="BH", min.nfeats=3, max.nfeats=1000){
  
  stopifnot(!is.null(names(gsets)))
  for(gset in gsets){
    stopifnot(gset %in% rownames(feat.tab))
  }
  
  # get G index
  index <- g_index(G=G, object=feat.tab, min.nfeats=min.nfeats, max.nfeats=max.nfeats)
  
  # make each index set as vector of two-level factors
  index.b <- lapply(index, FUN=function(x) factor(-as.numeric(rownames(feat.tab) %in% x), levels=c(-1,0)))
  ngenes <- vapply(index, FUN=length, FUN.VALUE=numeric(1))
  
  # run fisher.test for each test set
  for(i in seq_along(gsets)){
    gset.b <- factor(-as.numeric(rownames(feat.tab) %in% gsets[[i]]), levels=c(-1,0))
    res.tmp <- as.data.frame(t(vapply(index.b, FUN=function(x.b){
      tb <- table(x.b, gset.b)
      c(num=tb[1,1], p= stats::fisher.test(tb, alternative="greater")$p.value)
    }, FUN.VALUE = numeric(2))))
    
    res.tmp$FDR <- stats::p.adjust(res.tmp$p, method=adjust.method)
    colnames(res.tmp) <- paste(names(gsets)[i], colnames(res.tmp), sep = '.')
    
    if(i == 1) res <- cbind(NGenes=ngenes, res.tmp) else res <- cbind(res, res.tmp[rownames(res), ])
  }
  
  # order rows by combined p-values
  res <- res[order(combine_pvalues(res)), ]
  rownames(res) <- clean_filenames(rownames(res))
  
  # change FDR to appropriate adjustment name if user doesn"t use FDR
  if (!(adjust.method %in% c("BH", "fdr"))){
    colnames(res) <- gsub("FDR$", adjust.method, colnames(res))
  }
  
  res.xl <- df_signif(res, digits = 3)
  # write xlsx file with links
  if (!is.na(name)){
    nm <- paste(name, "fisher_test", sep="_")
    write_linked_xl(pwy.tab=res.xl, feat.lst=index, feat.tab=feat.tab, name=nm)
  }
  
  return(res)
}