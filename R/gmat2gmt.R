#' Transform analyte-by-pathway matrix into list compatible with Gene Matrix Transposed (GMT) file format
#' 
#' Transform analyte-by-pathway matrix into list compatible with Gene Matrix Transposed (GMT) file format, described at the \href{https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29}{Broad Institute}.
#'
#' @param gmat A numeric Matrix of features (proteins and/or metabolites) as rows and pathways as columns, indicating 
#' pathway membership with positive elements.
#' @return A list of analyte sets. Each element has a set \code{name}, \code{description}, and \code{genes}.
#' @export

gmat2gmt <- function(gmat){
  stopifnot(ncol(gmat) > 0, nrow(gmat) > 0, !is.null(colnames(gmat)), !is.null(rownames(gmat)))
  
  gmt.lst <- list()
  for (g.col in colnames(gmat)){
    gmt.lst[[g.col]]$name <- g.col
    #repeat name in description
    gmt.lst[[g.col]]$description <- g.col
    gmt.lst[[g.col]]$genes <- rownames(gmat)[gmat[,g.col] > 0]
  }
  return(gmt.lst)
}
