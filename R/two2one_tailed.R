#' Convert 2-tailed into 1-tailed p-values
#' 
#' Convert 2-tailed into 1-tailed p-values from a matrix-like object.
#' 
#' @param stat.col Index or \code{\link{regexp}} with column name  or column name suffix with numeric signed statistics
#' or with \code{"Up", "Down"} values.
#' @param p.col Index or \code{\link{regexp}} with column name suffix or full column names of p-value columns.
#' @param alternative Vector of direction of change: either \code{"greater"} or \code{"less"}. Also permitted are  
#' \code{"Up"} or \code{"Down"}, in which case the \code{stat.col} must include elements \code{"Up"} or \code{"Down"}.
#' @inheritParams combine_pvalues
#' @inheritParams prop_changed
#' @return p-values.

# alternative can't be 2-sided
two2one_tailed <- function(tab, stat.col="logFC|slope|r|Direction", p.col="p|PValue", 
                           alternative=c("greater", "less", "Up", "Down")){
  stopifnot(ncol(tab) >= 1, nrow(tab) >= 1, !is.null(colnames(tab)))
  alternative <- match.arg(alternative)
  p.colnm <- grep_cols(tab, p.cols=p.col)
  stat.colnm <- grep_cols(tab, stat.cols=stat.col)
  
  if (alternative %in% c("Up", "Down")){
    if (!all(tab[, stat.colnm] %in% c("Up", "Down"))){
      stop("alternative is ", alternative, "so stats elements must be 'Up' or 'Down'.")
    }
    ifelse(tab[, stat.colnm] == "Up", 1, -1)
  }
  alt.sign <- ifelse(alternative %in% c("Up", "greater"), 1, -1)
  sign.v <- sign(tab[, stat.colnm])
  
  new_pv <- stats::setNames(tab[, p.col]/2, nm=rownames(tab))
  if (any(sign.v != alt.sign)){
    opp.colnm <- which(sign.v != alt.sign)
    new_pv[opp.colnm] <- 1 - new_pv[opp.colnm]
  }
  return(new_pv)
}