#' Convert roast two-tailed results to one-tailed
#' 
#' Convert \code{mroast} or \code{fry} two-tailed results to one-tailed and remove 
#' non-directional, termed "Mixed", statistics.
#' 
#' @param roast.res Result from \code{mroast} or \code{fry}. 
#' @param fun function to use, either \code{fry} or \code{mroast}.
#' @param alternative indicates the alternative hypothesis and must be one of
#'  \code{"two.sided"}, \code{"greater"} or \code{"less"}. \code{"greater"}
#'  corresponds to positive association, \code{"less"} to negative association.
#' @param nrot number of rotations used to estimate the p-values for \code{mroast}.
#' @param adjust.method method used to adjust the p-values for multiple testing.
#' Only for \code{mroast}.
#' @return Modified result from \code{mroast} or \code{fry}.
#' @details This function is not meant to be called directly by the user.

roast_one_tailed <- function(roast.res, fun, alternative, nrot, adjust.method){
  direction <- sub("greater", "Up", sub("less", "Down", alternative))
  if (fun=="fry"){
    roast.res[,'PValue'] <- fry_two2one_tailed(tab = roast.res, direction = direction)
  } else {
    roast.res[,'PValue'] <- mroast_two2one_tailed(tab = roast.res, direction = direction, nrot = nrot)
  }
  roast.res[,'FDR'] <- stats::p.adjust(roast.res[,'PValue'], method = adjust.method)
  mixed.cols <- grep('Mixed', colnames(roast.res))
  if (length(mixed.cols)>0){ roast.res <- roast.res[,-mixed.cols] }
  return(roast.res)
}