#' Convert roast two-tailed results to one-tailed
#' 
#' Convert \code{mroast} or \code{fry} two-tailed results to one-tailed and remove 
#' non-directional, termed "Mixed", statistics.
#' 
#' @param roast.res Result from \code{mroast} or \code{fry}.
#' @inheritParams limma_contrasts
#' @inheritParams roast_contrasts
#' @return Modified result from \code{mroast} or \code{fry}.

roast_one_tailed <- function(roast.res, fun, alternative, nrot, adjust.method){
  direction <- sub("greater", "Up", sub("less", "Down", alternative))
  if (fun=="fry"){
    roast.res[,"PValue"] <- fry_two2one_tailed(tab = roast.res, direction = direction)
  } else {
    roast.res[,"PValue"] <- mroast_two2one_tailed(tab = roast.res, direction = direction, nrot = nrot)
  }
  roast.res[,"FDR"] <- stats::p.adjust(roast.res[,"PValue"], method = adjust.method)
  mixed.cols <- grep("Mixed", colnames(roast.res))
  if (length(mixed.cols)>0){ roast.res <- roast.res[,-mixed.cols] }
  return(roast.res)
}