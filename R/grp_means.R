#' Calculate group means per feature
#' 
#' Calculate group means per feature (e.g. per gene), accounting for possible intra-sample duplication.
#' 
#' @param args.lst List of arguments to be used in \code{do.call(limma::lmFit, args.list)}. It is only checked
#' whether \code{is.null(args.lst$ndups)} and \code{is.null(args.lst$spacing)}, which could be true or false.
#' @inheritParams limma_contrasts
#' @return Matrix of features (e.g. genes) by groups of means of object with column names ending in \code{.avg}.

# no need to export
# this is tested as part of limma_contrasts
grp_means <- function(object, grp, args.lst=list()){
  ndups <- spacing <- 1
  mtrx <- limma::getEAWP(object)$exprs
  stopifnot(nrow(mtrx)>1, ncol(mtrx)>1, ncol(mtrx)==length(grp), !is.null(dimnames(mtrx)), is.list(args.lst))
  printed <- ifelse(!is.matrix(object), yes=!is.null(object$printer), no=FALSE)
  
  if (!is.null(args.lst$ndups)){
    ndups <- args.lst$ndups
  } else if (printed && !is.null(object$printer$ndups)) ndups <- object$printer$ndups
  
  if (!is.null(args.lst$spacing)){
    spacing <- args.lst$spacing
  } else if (printed && !is.null(object$printer$spacing)) spacing <- object$printer$spacing
  
  stopifnot(is.numeric(ndups), ndups>=1, is.numeric(spacing), spacing>=1)
  
  if (ndups >= 2){
    rnms <- limma::uniquegenelist(rownames(mtrx), ndups = ndups, spacing = spacing)
    mtrx <- limma::unwrapdups(M=mtrx, ndups=ndups, spacing=spacing)
    rownames(mtrx) <- rnms
    grp <- rep(x=grp, each=ndups)
  }
  grp.means <- t(apply(X=mtrx, MARGIN=1, FUN=function(v) tapply(X=v, INDEX=grp, FUN=mean, na.rm=TRUE)))
  colnames(grp.means) <- paste(colnames(grp.means), "avg", sep=".")
  
  grp.means
}