#' Map pathway list gene IDs to another ID type
#' 
#' Map pathway list with old gene IDs (e.g. entrez IDs) to another ID type (e.g. gene symbols) using a table with both
#' annotation types. The original IDs can be separated within a column of annot with a string, e.g. " /// ".
#' The other IDs should be the rownames of annot.
#' 
#' @param annot annotation for the features that has a column of the same type as in gene set list *G*.
#' @param sep.str strings that separates symbols if there are multiple symbols for a feature.   
#' @param symbol.col column name or index for the symbol column in *annot*.  
#' @inheritParams roast_contrasts
#' @return Gene set list *G* with the symbols replaced by the mapped IDs.
#' @details Gene annotations are transformed to upper-case to avoid missing matches of different cases.
#' @export

map_glist <- function(G, annot, sep.str=" /// ", symbol.col="Gene.Symbol"){
  
  stopifnot(symbol.col %in% colnames(annot) || symbol.col %in% 1:ncol(annot), !is.null(rownames(annot)))
  
  sym.v <- toupper(annot[, symbol.col, drop = TRUE])
  sym.lst <- strsplit(sym.v, split = sep.str)
  row.nms <- rep(rownames(annot), times = vapply(sym.lst, FUN=length, FUN.VALUE = numeric(1)))
  map <- data.frame(Rowname = row.nms, Symbol = unlist(sym.lst), stringsAsFactors = FALSE)
  for(i in seq_along(G)){
    G[[i]][[3]] <- unique(map$Rowname[ map$Symbol %in% toupper(G[[i]][[3]]) ])
  }
  return(G)
}
