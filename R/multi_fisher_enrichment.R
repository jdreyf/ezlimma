#' Fisher Exact Test of gene set enrichment with multiple vectors of significant genes with output to Excel
#' 
#' Test enrichment of multiple vectors of significant genes against the universe of genes in \code{feat.tab} per 
#' gene set using \code{\link[stats]{fisher.test}}. It returns a data frame with statistics per gene set, and can 
#' write this to Excel. The Excel file links to CSV files, which contain statistics per genes in a set. If you only
#' have one vector of significant genes, use \code{\link[ezlimma]{fisher_enrichment}}.
#' 
#' @param sig.sets Named list whose elements are vectors of significant gene IDs matching 
#' \code{rownames(feat.tab)}.
#' @inheritParams roast_contrasts
#' @return List with two elements: \code{pwy.stats}, a table of pathway statistics; and \code{feat.tab}, a table that 
#' appends a binary matrix of which genes are in which \code{sig.set} and the input \code{feat.tab}.
#' @details Pathway (i.e. gene set) names are altered to be valid filenames in Windows and Linux. Numeric columns are
#' rounded to 3 significant figures.
#' @seealso \code{\link[ezlimma]{fisher_enrichment}}
#' @export

# name is for multi _fe output
multi_fisher_enrichment <- function(sig.sets, G, feat.tab, name=NA, adjust.method="BH", min.nfeats=3, max.nfeats=1000, 
                                    pwy.nchar=199){
  stopifnot(!duplicated(names(sig.sets)), length(sig.sets) > 1)
  
  # get G index
  index <- g_index(G=G, object=feat.tab, min.nfeats=min.nfeats, max.nfeats=max.nfeats)
   
  # want gene membership matrix, which will be subset by pathway for CSVs
  # fe = fisher enrichment
  feats.all <- unique(unlist(sig.sets))
  stopifnot(feats.all %in% rownames(feat.tab))
  fe.tab <- data.frame(matrix(0, nrow=length(feats.all), ncol=length(sig.sets), dimnames=list(feats.all, names(sig.sets))))
  
  for (ind in 1:length(sig.sets)){
    nm.ss <- names(sig.sets)[[ind]]
    sig.set <- list(sig.sets[[ind]])
    names(sig.set) <- nm.ss
    fe.tmp <- fisher_enrichment(sig.set=sig.set, G=G, feat.tab=feat.tab, name=NA, adjust.method=adjust.method, 
                      min.nfeats=min.nfeats, max.nfeats=max.nfeats, pwy.nchar=pwy.nchar)
    
    if (ind == 1){
      pwy.mat <- fe.tmp
    } else {
      pwy.mat <- cbind(pwy.mat, fe.tmp[rownames(pwy.mat), -1, drop=FALSE])
    }
    
    fe.tab[, nm.ss] <- as.numeric(rownames(fe.tab) %in% sig.set[[1]])
  }
  fe.df <- data.frame(fe.tab[rownames(feat.tab),, drop=FALSE], feat.tab)

  # writexl
  res.xl <- signif(pwy.mat, digits = 3)
  # write xlsx file with links
  if (!is.na(name)){
    nm <- paste(name, "fisher_test", sep="_")
    write_linked_xl(pwy.tab=res.xl, feat.lst=index, feat.tab=fe.df, name=nm, pwy.nchar=pwy.nchar)
  }
  
  return(list(pwy.stats=pwy.mat, feat.tab=fe.df))
}