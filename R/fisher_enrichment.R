#' Fisher Exact Test of gene set enrichment with output to Excel
#' 
#' Test enrichment of one or more vectors of significant genes against the universe of genes in \code{feat.tab} per 
#' gene set using \code{\link[stats]{fisher.test}}. It returns a data frame with statistics per gene set, and can 
#' write this to Excel. The Excel file links to CSV files, which contain statistics per genes in a set.
#' 
#' @param sig.set Named list of length one whose sole element is a vector of significant gene IDs matching 
#' \code{rownames(feat.tab)}.
#' @param return.lst Logical; if \code{FALSE} returns table of pathway statistics, whereas if \code{TRUE} returns
#' list whose 1st element is table of pathway statistics and 2nd element is matrix whose column \code{i} is
#' \code{sig.set[[i]] G[[i]]$genes)}. Meant for internal use.  
#' @inheritParams roast_contrasts
#' @return Data frame of gene set statistics.
#' @details Pathway (i.e. gene set) names are altered to be valid filenames in Windows and Linux. Numeric columns are
#' rounded to 3 significant figures.
#' @examples
#'  G = list(s1=list(name = "s1", description=NULL, genes=letters[3:10]), 
#'    s2=list(name = "s2", description=NULL, genes=letters[2:99]))
#'  feat.tab <- matrix(rnorm(18), ncol=3, dimnames=list(letters[1:6], paste0("col", 1:3)))
#'  fisher_enrichment(sig.set = list(A=letters[1:3]), G = G, feat.tab = feat.tab)
#' @export

fisher_enrichment <- function(sig.set, G, feat.tab, name=NA, adjust.method="BH", min.nfeats=3, max.nfeats=1000, 
                              return.lst=FALSE, pwy.nchar=199){
  stopifnot(!is.null(names(sig.set)), !is.null(feat.tab), sig.set[[1]] %in% rownames(feat.tab), length(sig.set) == 1)
  
  # get G index
  index <- g_index(G=G, object=feat.tab, min.nfeats=min.nfeats, max.nfeats=max.nfeats)
  
  # make each index set as vector of two-level factors
  index.b <- lapply(index, FUN=function(x) factor(-as.numeric(rownames(feat.tab) %in% x), levels=c(-1, 0)))
  ngenes <- vapply(index, FUN=length, FUN.VALUE=numeric(1))
  
  # run fisher.test for each test set
  ind <- 1
  gset.b <- factor(-as.numeric(rownames(feat.tab) %in% sig.set[[ind]]), levels=c(-1,0))
  res.va <- vapply(index.b, FUN=function(x.b){
    tb <- table(x.b, gset.b)
    c(num=tb[1,1], p=stats::fisher.test(tb, alternative="greater")$p.value)
  }, FUN.VALUE = numeric(2))
  res.tmp <- as.data.frame(t(res.va))
  res.tmp$FDR <- stats::p.adjust(res.tmp$p, method=adjust.method)
  colnames(res.tmp) <- paste(names(sig.set)[ind], colnames(res.tmp), sep = '.')
  res <- cbind(NGenes=ngenes, res.tmp)
  
  if (return.lst){
    gene.membership <- vapply(index.b, FUN = function(xx){
      # -1 indicates membership
      xx != 0 & gset.b != 0
    }, FUN.VALUE = numeric(nrow(feat.tab)))
    rownames(gene.membership) <- rownames(feat.tab)
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
    write_linked_xl(pwy.tab=res.xl, feat.lst=index, feat.tab=feat.tab, name=nm, pwy.nchar=pwy.nchar)
  }
  
  if (!return.lst){
    return(res)
  } else {
    res.lst <- list(pwy.stats = res, gene.membership = gene.membership)
    return(res.lst)
  }
}