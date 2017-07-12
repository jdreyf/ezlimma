#' Transform a log2 fold-change to a fold-change
#'
#' Transform a possibly negative log2 fold-change to a fold-change.
#' 
#' @param logFC A numeric log2 fold-change.
#' @return The fold-change.
#' @details This function assumes that log fold-changes less than 1 were converted 
#' to their negative multiplicative reciprocal, e.g. lfc = 0.5 -> lfc = -2, as is
#' commonly done in biology.
#' @examples
#' logfc2fc(2) #4
#' logfc2fc(-2) #-4

#note: function doesn't need to be exported
logfc2fc <- function(logFC){
  #sign is -1 if logFC<0; 1 if logFC>=0
  sgn <- (-1)^(1+as.numeric(logFC>=0))
  fc <- sgn*2^abs(logFC)
  return(fc)
}