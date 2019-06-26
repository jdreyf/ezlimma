#' Gets the mapped IDs (rownames) of features (e.g. genes) for a gene set list of a particular type (e.g. gene symbols)
#' 
#' @param G Gene set list as returned from \code{\link{read_gmt}}.
#' @param annot annotation for the features that has a column of the same type as in gene set list *G*.
#' @param sep.str strings that separates symbols if there are multiple symbols for a feature.   
#' @param symbol.col column name or index for the symbol column in *annot*.  
#' @return Gene set list *G* with the symbols replaced by the mapped IDs.
#' @export

map_glist <- function(G, annot, sep.str=" /// ", symbol.col="Gene.Symbol"){
  
  stopifnot(symbol.col %in% colnames(annot) || symbol.col %in% 1:ncol(annot))
  
  sym.v <- toupper(annot[, symbol.col, drop = TRUE])
  sym.lst <- strsplit(sym.v, split = sep.str)
  row.nms <-rep(rownames(annot), times = vapply(sym.lst, FUN=length, FUN.VALUE = numeric(1)))
  map <- data.frame(Rowname = row.nms, Symbol = unlist(sym.lst))
  for(i in seq_along(G)){
    G[[i]][[3]] <- unique(map$Rowname[ map$Symbol %in% toupper(G[[i]][[3]]) ])
  }
  return(G)
}
