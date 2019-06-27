#' A wrapper function for \code{\link[limma]{geneSetTest}} with output to Excel
#' 
#' Test whether a set of genes is highly ranked. It returns a data frame with statistics per gene set, and writes this to an Excel file. 
#' The Excel file links to CSV files, which contain statistics per gene set. 
#' 
#' @param gsats A named vectoror, matrix or data.frame of genewise statistic (e.g. z-scores, t-statistics) by which genes can be ranked. 
#' The names or rownames should be the same as the rownames of  *feat.tab*
#' @inheritParams roast_contrasts
#' @inheritParams limma::geneSetTest
#' @return Data frame of gene set statistics.
#' @details Pathway (i.e. gene set) names are altered to be valid filenames in Windows and Linux. Numeric columns are
#' rounded to 3 significant figures.
#' @export

multi.gst <- function(gstats, G, feat.tab, name=NA, adjust.method="BH", alternative="mixed", type="auto", adjust.method ="BH", 
                      min.nfeats=3, max.nfeats=1000, ranks.only=TRUE, nsim=9999){
  
  if (is.vector(gstats)){ 
    gstats <- as.matrix(gstats)
    colnames(gstats) <- "gst"
    }
  stopifnot(!is.null(rownames(gstats)), rownames(gstats) %in% rownames(feat.tab))
  
  # get G index
  index <- g_index(G=G, object=gstats, min.nfeats=min.nfeats, max.nfeats=max.nfeats)
  
  # run geneSetTest
  for(i in 1:ncol(gstats)){
    gstat <- gstats[, i, drop=TRUE]
    names(gstat) <- rownames(gstats)
    
    pval <- vapply(index, function(x){
      limma::geneSetTest(x, statistics=gstat, alternative=alternative, type=type, ranks.only=ranks.only, nsim=nsim)
    })
    
    # some pval may > 1
    # the problems arises because the two-sided test method is to double the smaller of the two tail probabilities.  If both of the tail
    # probabilities are greater than 0.5, then the final p-value ends up being > 1.
    pval[pval > 1] <- 1
    
    fdr <- stats::p.adjust(pval, method=adjust.method)
    
    if(alternative == "mixed"){
      res.tmp <- data.frame(p=pval, FDR=fdr)
    } else {
      avg <- vapply(index, FUN=function(x) mean(gstat[x]), FUN.VALUE = numeric(1))
      direction <- ifelse(avg > 0, "Up", "Down")
      res.tmp <- data.frame(Direction=direction, p=pval, FDR=fdr)
    }
    
    colnames(res.tmp) <- paste(colnames(gstats)[i], colnames(res.tmp), sep=".")
    
    ngenes <- Vapply(index, FUN=length, FUN.VALUE=numeric(1))
    if(i == 1) res <- cbind(NGenes=ngenes, res.tmp)
    else res <- cbind(res, res.tmp[rownames(res), ])
  }
  
  # order rows by combined p-values
  res <- res[order(combine_pvalues(res)), ]
  rownames(res) <- clean_filenames(rownames(res))
  
  # change FDR to appropriate adjustment name if user doesn"t use FDR
  if (!(adjust.method %in% c("BH", "fdr"))){
    colnames(res) <- gsub("FDR$", adjust.method, colnames(res))
  }
  
  res.xl <- df_signif(res, digits=3)
  # write xlsx file with links
  if (!is.na(name)){
    nm <- paste(name, "gst", sep="_")
    write_linked_xl(pwy.tab=res.xl, feat.lst=index, feat.tab=feat.tab, name=nm)
  }
  
  return(res)
}
