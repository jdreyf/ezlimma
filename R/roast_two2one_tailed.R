#' Convert roast two-tailed results to one-tailed
#' 
#' Convert \code{mroast} or \code{fry} two-tailed results to one-tailed and remove 
#' non-directional, termed "Mixed", statistics.
#' 
#' @param roast.res Result from \code{mroast} or \code{fry}.
#' @inheritParams limma_contrasts
#' @inheritParams roast_contrasts
#' @return Modified result from \code{mroast} or \code{fry}.

roast_two2one_tailed <- function(roast.res, fun, alternative, nrot, adjust.method){
  if (fun=="fry"){
    # two2one_tailed returns a matrix
    roast.res[,"PValue"] <- two2one_tailed(tab = roast.res, p.cols="PValue", stat.cols="Direction", 
                                           alternative = alternative)[,1]
  } else {
    roast.res[,"PValue"] <- two2one_tailed(tab = roast.res, p.cols="PValue", stat.cols="Direction", 
                                           alternative = alternative, nperm = nrot)[,1]
  }
  roast.res[,"FDR"] <- stats::p.adjust(roast.res[,"PValue"], method = adjust.method)
  mixed.cols <- grep("Mixed", colnames(roast.res))
  if (length(mixed.cols) > 0){ roast.res <- roast.res[, -mixed.cols] }
  return(roast.res)
}