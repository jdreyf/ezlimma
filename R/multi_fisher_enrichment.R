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
#' @return Data frame of gene set statistics.
#' @details Pathway (i.e. gene set) names are altered to be valid filenames in Windows and Linux. Numeric columns are
#' rounded to 3 significant figures.
#' @seealso \code{\link[ezlimma]{fisher_enrichment}}
#' @export

# name is for multi _fe output
multi_fisher_enrichment <- function(sig.sets, G, feat.tab, name=NA, adjust.method="BH", min.nfeats=3, max.nfeats=1000, 
                                    pwy.nchar=199){
  stopifnot(!duplicated(names(sig.sets)), length(sig.sets) > 1)
   
  fe.lst <- list()
  for (ind in 1:length(sig.sets)){
    nm <- names(sig.sets)[[ind]]
    sig.set <- list(sig.sets[[ind]])
    names(sig.set) <- nm
    fe.lst[[nm]] <- fisher_enrichment(sig.set=sig.set, G=G, feat.tab=feat.tab, name=NA, adjust.method=adjust.method, 
                      min.nfeats=min.nfeats, max.nfeats=max.nfeats, return.lst=TRUE, pwy.nchar=pwy.nchar)
    
    if (ind == 1){
      pwy.mat <- fe.lst[[nm]]$pwy.stats
    } else {
      pwy.mat <- cbind(pwy.mat, fe.lst[[nm]]$pwy.stats[rownames(pwy.mat), -1])
    }
    
    # assume diff sig.sets return same pwys? yes, since subset in fisher_enrich by overlap of rownames(feat.tab) & G
    # assume diff sig.sets return same genes? yes, should all return rownames(feat.tab)
    if (ind == 1){
      gene.mem.lst <- list()
      gene.mem.tmp <- fe.lst[[nm]]$gene.membership
      pwy.nms <- colnames(gene.mem.tmp)
      for (pwy.tmp in pwy.nms){
        gene.mem.lst[[pwy.tmp]] <- gene.mem.tmp[, pwy.tmp, drop=FALSE]
        colnames(gene.mem.lst[[pwy.tmp]]) <- nm
      }
    } else {
      gene.mem.tmp <- fe.lst[[nm]]$gene.membership
      pwy.nms <- colnames(gene.mem.tmp)
      for (pwy.tmp in pwy.nms){
        gene.mem.lst[[pwy.tmp]] <- cbind(gene.mem.lst[[pwy.tmp]], 
                                         gene.mem.tmp[rownames(gene.mem.lst[[pwy.tmp]]), pwy.tmp, drop=FALSE])
        colnames(gene.mem.lst[[pwy.tmp]])[ncol(gene.mem.lst[[pwy.tmp]])] <- nm
      }
    }
  } # end for ind
  # writexl
  return(list(pwy.stats=pwy.mat, gene.membership=gene.mem.lst))
}