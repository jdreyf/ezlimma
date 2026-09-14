#' Fisher Exact Test of gene set enrichment with output to Excel
#' 
#' Test enrichment of one or more vectors of significant genes against the universe of genes in \code{feat.tab} per 
#' gene set using \code{\link[stats]{fisher.test}}. It returns a data frame with statistics per gene set, and can 
#' write this to Excel. The Excel file links to CSV files, which contain statistics per genes in a set.
#' 
#' @param sig.set Named list of length one whose sole element is a vector of significant gene IDs matching 
#' \code{rownames(feat.tab)}.
#' @inheritParams roast_contrasts
#' @return Table of pathway statistics with the number of genes from \code{feat.tab} in the pathway, the number of these genes that are 
#' in \code{sig.set}, the p-value, and the adjusted p-value from the one-sided Fisher exact test.
#' @details Pathway (i.e. gene set) names are altered to be valid filenames in Windows and Linux. Numeric columns are
#' rounded to 8 significant figures.
#' @examples
#'  G = list(s1=list(name = "s1", description=NULL, genes=letters[3:10]), 
#'    s2=list(name = "s2", description=NULL, genes=letters[2:99]))
#'  feat.tab <- matrix(rnorm(18), ncol=3, dimnames=list(letters[1:6], paste0("col", 1:3)))
#'  fisher_enrichment(sig.set = list(A=letters[1:3]), G = G, feat.tab = feat.tab)
#' @export

fisher_enrichment <- function(sig.set, G, feat.tab, name=NA, adjust.method="BH", min.nfeats=3, max.nfeats=1000,
                              pwy.nchar=199){
  stopifnot(!is.null(names(sig.set)), !is.null(feat.tab), sig.set[[1]] %in% rownames(feat.tab), length(sig.set) == 1)

  # get G index
  index <- g_index(G=G, object=feat.tab, min.nfeats=min.nfeats, max.nfeats=max.nfeats)
  .fisher_enrichment_index(sig.set=sig.set, index=index, feat.tab=feat.tab, name=name, adjust.method=adjust.method,
                           pwy.nchar=pwy.nchar)
}

#' Core of \code{fisher_enrichment}, taking a pre-computed gene set \code{index}
#'
#' Not exported. Split out of \code{\link{fisher_enrichment}} so that \code{\link{multi_fisher_enrichment}} can
#' resolve gene set membership once (via \code{\link{g_index}}) and reuse it across all \code{sig.sets}, rather
#' than redoing that \code{nrow(feat.tab) * length(G)}-scale work for every one.
#'
#' Also avoids materializing a length-\code{nrow(feat.tab)} indicator vector per gene set (as the original
#' implementation did, via \code{table()} on two such vectors): each 2x2 contingency table only needs the
#' set size, the \code{sig.set} size, and their intersection size, so intersection counts for all sets are
#' computed with a single vectorized group-by over the (much smaller) total gene set membership instead.
#'
#' @inheritParams fisher_enrichment
#' @param index Gene set list as returned by \code{\link{g_index}}.
#' @noRd

.fisher_enrichment_index <- function(sig.set, index, feat.tab, name=NA, adjust.method="BH", pwy.nchar=199){
  n.universe <- nrow(feat.tab)
  ngenes <- lengths(index)

  ind <- 1
  sig.v <- sig.set[[ind]]
  n.sig <- sum(rownames(feat.tab) %in% sig.v)

  # intersection size per set, via one group-by over all set memberships, instead of one O(nrow(feat.tab))
  # indicator vector per set
  nsets <- length(index)
  probe.v <- unlist(index, use.names = FALSE)
  set.v <- rep.int(seq_len(nsets), lengths(index))
  in.sig <- probe.v %in% sig.v
  a <- as.vector(tapply(in.sig, factor(set.v, levels=seq_len(nsets)), sum))
  a[is.na(a)] <- 0

  # 2x2 contingency table per set: a = in-set & in-sig, b = in-set & not-in-sig,
  # c = not-in-set & in-sig, d = not-in-set & not-in-sig
  b <- ngenes - a
  cc <- n.sig - a
  d <- n.universe - ngenes - n.sig + a

  pv <- vapply(seq_len(nsets), FUN=function(i){
    tb <- matrix(c(a[i], cc[i], b[i], d[i]), nrow=2)
    stats::fisher.test(tb, alternative="greater")$p.value
  }, FUN.VALUE = numeric(1))

  res.tmp <- data.frame(N.DE=a, p=pv, row.names=names(index))
  res.tmp$FDR <- stats::p.adjust(res.tmp$p, method=adjust.method)
  colnames(res.tmp) <- paste(names(sig.set)[ind], colnames(res.tmp), sep = '.')
  res <- cbind(NGenes=ngenes, res.tmp)

  # order rows by combined p-values
  res <- res[order(combine_pvalues(res)), ]

  # change FDR to appropriate adjustment name if user doesn't use FDR
  if (!(adjust.method %in% c("BH", "fdr"))){
    colnames(res) <- gsub("FDR$", adjust.method, colnames(res))
  }

  res.xl <- df_signif(res, digits = 8)
  # write xlsx file with links
  if (!is.na(name)){
    nm <- paste(name, "fisher_test", sep="_")
    write_linked_xl(pwy.tab=res.xl, feat.lst=index, feat.tab=feat.tab, name=nm, pwy.nchar=pwy.nchar)
  }

  return(res)
}