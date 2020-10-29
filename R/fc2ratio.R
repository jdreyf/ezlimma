#' Transform fold-change (FC) into a ratio
#'
#' Transform fold-change (FC) into a ratio. It's assumed \pkg{ezlimma} transformed a ratio less than 1 into
#' a negative FC, e.g. a ratio of 0.5 gave a FC of -2. This function reverses that operation.
#'
#' @param FC Numeric vector of fold-changes.
#' @return Vector of corresponding ratios.
#' @export

fc2ratio <- function(FC){
  ret <- rep(NA, length(FC))
  fc.pos <- which(FC >= 0)
  if (length(fc.pos) > 0) ret[fc.pos] <- FC[fc.pos]-1
  fc.neg <- which(FC < 0)
  if (length(fc.neg) > 0) ret[fc.neg] <- -1*(1 - 1/abs(FC[fc.neg]))
  ret
}
