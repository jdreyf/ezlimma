#' Convert 2-tailed into 1-tailed p-values
#' 
#' Convert 2-tailed into 1-tailed p-values, possibly using signs of associated statistics, from a matrix-like object.
#' 
#' @param alternative Direction of change: either \code{"greater"} or \code{"less"}, or their synonyms  
#' \code{"Up"} or \code{"Down"}.
#' @param nperm Number of rotations used to estimate the p-values from \code{\link[limma]{roast}}, or number of 
#' permutations, where p-value calculated according to Phipson & Smyth (2010). Ignored if \code{NULL}.
#' @inheritParams combine_pvalues
#' @inheritParams prop_changed
#' @return Matrix of p-values.
#' @references Phipson & Smyth (2010). Permutation P-values Should Never Be Zero: Calculating Exact P-values When 
#' Permutations Are Randomly Drawn. Statistical Applications in Genetics and Molecular Biology, 9(1).

# alternative can't be 2-sided
# nperm:
# pv=(b+1)/(nrot+1) where b = number of rotations giving a more extreme stat
# if dir.col==direction, then by symmetry half of the previously extreme rotations 
# are still extreme, so pv'=(b/2 + 1)/(nrot+1)
# if dir.col!=direction, then nrot-b/2 of the previously extreme rotations are extreme,
# which is the equivalent of 1-p/2, so pv' = (nrot - b/2 + 1) / (nrot+1)
# this is also consistent with 0.5<=pv'<=1

two2one_tailed <- function(tab, p.cols="p|PValue", stat.cols="logFC|slope|cor|rho|Direction", 
                           alternative=c("greater", "less", "Up", "Down"), nperm=NULL){
  stopifnot(ncol(tab) >= 1, nrow(tab) >= 1, !is.null(colnames(tab)))
  alternative <- match.arg(alternative)
  p.colnms <- grep_cols(tab, p.cols=p.cols)
  stat.colnms <- grep_cols(tab, stat.cols=stat.cols)
  # columns either numeric or up/down
  # if matrix, must include numeric b/c of p-values
  if (is.data.frame(tab)){
    stat.char.colnms <- stat.colnms[!sapply(tab[, stat.colnms, drop=FALSE], is.numeric)]
    if (length(stat.char.colnms) > 0){
      tab[, stat.char.colnms] <- apply(tab[, stat.char.colnms, drop=FALSE], MARGIN=2, FUN=function(v){
        if (all(v %in% c("Up", "Down"))){
          ifelse(v == "Up", 1, -1)
        } else {
          if (any(is.na(suppressWarnings(as.numeric(v))))){
            stop("All stat columns must be numeric or have elements 'Up' or 'Down'.", call. = FALSE)
          } else {
            as.numeric(v)
          }
        }
      })
    } # stat.char.colnms
    tab <- data.matrix(tab[, colnames(tab) %in% c(stat.colnms, p.colnms)])
  } # end df
  alt.sign <- ifelse(alternative %in% c("Up", "greater"), 1, -1)
  
  sign.mat <- sign(data.matrix(tab[, stat.colnms, drop=FALSE]))
  
  # sign=0 associated p-values are not transformed to 1-new_pv 
  if (!is.null(nperm)){
    stopifnot(nperm > 0)
    b.mat <- tab[, p.colnms, drop=FALSE]*(nperm+1)-1
    # initialize as if dir.col==direction
    new_pv <- matrix((b.mat/2+1) / (nperm+1), nrow=nrow(tab), dimnames=list(rownames(tab), p.colnms))
    if (any(sign.mat != alt.sign)){
      opp.ind <- which(sign.mat != alt.sign)
      new_pv[opp.ind] <- (nperm - b.mat[opp.ind]/2 + 1)/(nperm+1)
    }
  } else {
    new_pv <- matrix(tab[, p.colnms, drop=FALSE]/2, nrow=nrow(tab), dimnames=list(rownames(tab), p.colnms))
    if (any(sign.mat != alt.sign)){
      opp.ind <- which(sign.mat != alt.sign)
      new_pv[opp.ind] <- 1 - new_pv[opp.ind]
    }
  }
  new_pv
}