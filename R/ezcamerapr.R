#' A wrapper function for \code{cameraPR} with output to Excel
#' 
#' Test whether a set of genes is highly ranked relative to other genes in terms of differential expression, 
#' accounting for inter-gene correlation with \code{\link[limma:camera]{cameraPR}}.
#' To find pathways whose genes have large magnitude changes, independent of direction of their change, 
#' provide \code{abs} of gene stats, and set \code{alternative="Up"}.
#' It returns a data frame with statistics per gene set, and writes this to an Excel file. 
#' The Excel file links to CSV files, which contain statistics per gene set. 
#' 
#' @param stats.tab A matrix-like data object with gene row names & named columns of numeric gene-wise statistics 
#' (e.g. z-scores, t-statistics) by which genes can be ranked. 
#' The row names should be the same as the row names of \code{feat.tab}.
#' All values must be \code{\link[base:is.finite]{finite}}.
#' @param alternative Alternative hypothesis; must be one of\code{"two.sided"}; \code{"greater"} or \code{"less"},
#' or their synonyms  \code{"Up"} or \code{"Down"}.
#' @param inter.gene.cor Numeric inter-gene correlation within tested sets. Can be estimated with 
#' \code{\link[limma:camera]{interGeneCorrelation}}.
#' @inheritParams limma::camera
#' @inheritParams roast_contrasts
#' @return Data frame of gene set statistics.
#' @details Pathway (i.e. gene set) names are altered to be valid filenames in Windows and Linux. Numeric columns are
#' rounded to 8 significant figures.
#'
#' For performance, gene set membership is resolved to integer indices into \code{rownames(stats.tab)} once
#' up front, and \code{\link[limma:camera]{cameraPR}} is called once per column of \code{stats.tab} for all
#' gene sets together, rather than once per gene set. \code{cameraPR} re-derives integer indices from character
#' IDs (\code{which(ID \%in\% iset)}) and recomputes the mean and variance of the full statistic vector on every
#' call it is given, so calling it once per gene set repeats both of those over the full, often very large,
#' feature universe for every gene set and every column. Passing pre-resolved integer indices and batching all
#' gene sets into a single call per column avoids that redundant work, which otherwise scales with
#' \code{nrow(stats.tab) * length(G) * ncol(stats.tab)} and can become prohibitive for large arrays
#' (e.g. hundreds of thousands of probes) tested against thousands of gene sets across many comparisons.
#' @export

ezcamerapr <- function(stats.tab, G, feat.tab, name=NA, adjust.method ="BH", alternative=c("two.sided", "greater", "less", "Up", "Down"),
                      min.nfeats=3, max.nfeats=1000, inter.gene.cor=0.01, pwy.nchar=199){
  alternative <- match.arg(alternative)
  if (is.data.frame(stats.tab)){ stats.tab <- as.matrix(stats.tab) }
  stopifnot(!is.null(rownames(stats.tab)), !is.null(colnames(stats.tab)), rownames(stats.tab) %in% rownames(feat.tab),
            is.finite(stats.tab), is.numeric(stats.tab))

  # stats.tab must be matrix
  index <- g_index(G=G, object=stats.tab, min.nfeats=min.nfeats, max.nfeats=max.nfeats)
  # resolve gene set membership to integer indices into rownames(stats.tab) once, since this mapping is the
  # same for every column of stats.tab -- see performance note in @details above
  index.int <- lapply(index, function(nm) match(nm, rownames(stats.tab)))

  for (col.ind in 1:ncol(stats.tab)){
    stats.tab.v <- stats.tab[, col.ind]
    # test all gene sets in a single cameraPR() call, rather than one call per gene set
    tmp <- limma::cameraPR(statistic=stats.tab.v, index=index.int, inter.gene.cor=inter.gene.cor, sort=FALSE)
    tmp$Direction <- ifelse(tmp$Direction == "Up", yes = 1, no = -1)
    # drop cameraPR's own FDR column (if present); it's recomputed below using adjust.method
    tab.tmp <- tmp[, c("NGenes", "Direction", "PValue")]
    colnames(tab.tmp) <- c("NGenes", "Direction", "p")
    if (alternative!="two.sided"){
      tab.tmp$p <- two2one_tailed(tab=tab.tmp, alternative = alternative)[,1]
    }
    tab.tmp$FDR <- stats::p.adjust(tab.tmp$p, method=adjust.method)
    # change FDR to appropriate adjustment name if user doesn't use FDR
    if (!(adjust.method %in% c("BH", "fdr"))){
      colnames(tab.tmp) <- gsub("FDR$", adjust.method, colnames(tab.tmp))
    }
    # don't name ngenes
    colnames(tab.tmp)[-1] <- paste(colnames(stats.tab)[col.ind], colnames(tab.tmp)[-1], sep=".")
    # NGenes must be conserved, b/c all stats must be finite
    if (col.ind == 1) tab <- tab.tmp else tab <- data.frame(tab, tab.tmp[rownames(tab), setdiff(colnames(tab.tmp), "NGenes")])
  }
  # order rows by p-values
  tab <- tab[order(combine_pvalues(tab)), ]
  
  res.xl <- df_signif(as.data.frame(tab), digits=8)
  # write xlsx file with links
  if (!is.na(name)){
    nm <- paste(name, "cameraPR", sep="_")
    write_linked_xl(pwy.tab=res.xl, feat.lst=index, feat.tab=feat.tab, name=nm, pwy.nchar=pwy.nchar)
  }
  return(tab)
}