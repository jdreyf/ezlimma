#' A wrapper function for \code{limma::cameraPR} with output to Excel
#' 
#' 
#' Test whether a set of genes is highly ranked relative to other genes in terms of differential expression, 
#' accounting for inter-gene correlation \code{\link[limma]{camera}}.
#' To find pathways whose genes have large magnitude changes, independent of direction of their change, 
#' provide \code{abs} of gene stats, and set \code{alternative="Up"}.
#' It returns a data frame with statistics per gene set, and writes this to an Excel file. 
#' The Excel file links to CSV files, which contain statistics per gene set. 
#' 
#' @param gstats A nmatrix-like data object with rownames & one column of genewise statistic (e.g. z-scores, t-statistics) 
#' by which genes can be ranked. 
#' The names should be the same as the rownames of  *feat.tab*
#' @param alternative Alternative hypothesis; must be one of\code{"two.sided"}; \code{"greater"} or \code{"less"},
#' or their synonyms  \code{"Up"} or \code{"Down"}.
#' @inheritParams roast_contrasts
#' @inheritParams limma::geneSetTest
#' @return Data frame of gene set statistics.
#' @details Pathway (i.e. gene set) names are altered to be valid filenames in Windows and Linux. Numeric columns are
#' rounded to 3 significant figures.
#' @export

ezcamerapr <- function(gstats, G, feat.tab, name=NA, adjust.method ="BH", alternative=c("two.sided", "greater", "less", "Up", "Down"),
                      min.nfeats=3, max.nfeats=1000){
  alternative <- match.arg(alternative)
  
  # in case it's a dataframe, which is a problem in sapply
  gstats <- data.matrix(gstats)
  stopifnot(!is.null(rownames(gstats)), rownames(gstats) %in% rownames(feat.tab))
  
  # gstats must be matrix
  index <- g_index(G=G, object=gstats, min.nfeats=min.nfeats, max.nfeats=max.nfeats)

  gstats.v <- stats::setNames(gstats, nm=rownames(gstats))
  tab <- t(vapply(index, FUN=function(xx){
    # gstats must be vector
    tmp <- limma::cameraPR(statistic=gstats.v, index=xx)
    tmp$Direction <- ifelse(tmp$Direction == "Up", 1, -1)
    data.matrix(tmp)
  }, FUN.VALUE = stats::setNames(numeric(3), nm=c("NGenes", "Direction", "p"))))
  tab <- as.data.frame(tab)
  
  if (alternative!="two.sided"){
    tab$p <- two2one_tailed(tab=tab, alternative = alternative)
  }
  
  tab$FDR <- p.adjust(tab$p, method=adjust.method)
  # change FDR to appropriate adjustment name if user doesn't use FDR
  if (!(adjust.method %in% c("BH", "fdr"))){
    colnames(tab) <- gsub("FDR$", adjust.method, colnames(tab))
  }
  # order rows by p-values
  tab <- tab[order(tab$p), ]
  
  res.xl <- df_signif(as.data.frame(tab), digits=3)
  # write xlsx file with links
  if (!is.na(name)){
    nm <- paste(name, "cameraPR", sep="_")
    write_linked_xl(pwy.tab=res.xl, feat.lst=index, feat.tab=feat.tab, name=nm)
  }
  return(tab)
}