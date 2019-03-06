#' Transform a log2 fold-change to a fold-change
#'
#' Transform a log2 fold-change to a (possibly negative) fold-change.
#' 
#' @param logFC Numeric vector of log2 fold-changes.
#' @return Fold-change.
#' @details This function converts log fold-changes < 0 to their negative multiplicative reciprocal, e.g. 
#' fc = 0.5 -> fc = -2, as is commonly done in biology. This function is not meant to be called directly by the user.
#' @examples 
#' \dontrun{
#' # not run
#' logfc2fc(2) #4
#' logfc2fc(-2) #-4
#' }

logfc2fc <- function(logFC){
  # sign is -1 if logFC<0; 1 if logFC>=0
  sgn <- (-1)^(1+as.numeric(logFC>=0))
  fc <- sgn*2^abs(logFC)
  return(fc)
}