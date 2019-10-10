#' Fisher Exact Test of gene set enrichment with output to Excel
#' 
#' Test enrichment of one or more vectors of significant genes against the universe of genes in \code{feat.tab} per 
#' gene set using \code{\link[stats]{fisher.test}}. It returns a data frame with statistics per gene set, and can 
#' write this to Excel. The Excel file links to CSV files, which contain statistics per genes in a set.
#' 
#' @param sig.sets Named list whose elements are a vector of significant gene IDs matching \code{rownames(feat.tab)}.
#' @inheritParams roast_contrasts
#' @return Data frame of gene set statistics.
#' @details Pathway (i.e. gene set) names are altered to be valid filenames in Windows and Linux. Numeric columns are
#' rounded to 3 significant figures.
#' 
#' \code{feat.tab} is needed if \code{!is.na(name)}.
#' @export

fisher_enrichment <- function(sig.sets, G, feat.tab=NULL, name=NA, adjust.method="BH", min.nfeats=3, max.nfeats=1000){
  stopifnot(!is.null(names(sig.sets)), is.na(name) || !is.null(feat.tab))
  for (gset in sig.sets){ stopifnot(gset %in% rownames(feat.tab)) }
  
  # get G index
  index <- g_index(G=G, object=feat.tab, min.nfeats=min.nfeats, max.nfeats=max.nfeats)
  
  # make each index set as vector of two-level factors
  index.b <- lapply(index, FUN=function(x) factor(-as.numeric(rownames(feat.tab) %in% x), levels=c(-1, 0)))
  ngenes <- vapply(index, FUN=length, FUN.VALUE=numeric(1))
  
  # run fisher.test for each test set
  for(ind in seq_along(sig.sets)){
    gset.b <- factor(-as.numeric(rownames(feat.tab) %in% sig.sets[[ind]]), levels=c(-1,0))
    res.va <- vapply(index.b, FUN=function(x.b){
      tb <- table(x.b, gset.b)
      c(num=tb[1,1], p=stats::fisher.test(tb, alternative="greater")$p.value)
    }, FUN.VALUE = numeric(2))
    res.tmp <- as.data.frame(t(res.va))
    res.tmp$FDR <- stats::p.adjust(res.tmp$p, method=adjust.method)
    colnames(res.tmp) <- paste(names(sig.sets)[ind], colnames(res.tmp), sep = '.')
    
    if(ind == 1) res <- cbind(NGenes=ngenes, res.tmp) else res <- cbind(res, res.tmp[rownames(res), ])
  }
  
  # order rows by combined p-values
  res <- res[order(combine_pvalues(res)), ]
 
  # change FDR to appropriate adjustment name if user doesn't use FDR
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