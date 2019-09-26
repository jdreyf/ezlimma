#' Two one-sided test procedure with \pkg{limma} for statistical equivalence test
#' 
#' Two one-sided test procedure with \pkg{limma} for statistical equivalence test of one contrast.
#' 
#' @param tost.lfc Numeric; smallest logFC that could be of interest.
#' @inheritParams limma_contrasts
#' 
#' @details \code{length(contrast.v)} must be 1. This function is based on \code{\link{limma_contrasts}}, so see there
#' for more details.
#' @seealso \code{\link{limma_contrasts}}.
#' @export

limma_tost <- function(object, grp=NULL, contrast.v, tost.lfc, design=NULL, weights=NA, trend=FALSE, block=NULL, 
                       correlation=NULL, adjust.method="BH", add.means=!is.null(grp), check.names=TRUE){
  
  stopifnot(all(is.na(weights)) || is.null(weights) || dim(weights)==dim(object) || length(weights)==nrow(object) || 
              length(weights)==ncol(object), length(contrast.v)==1, tost.lfc>0)
  if (is.vector(object)) stop("'object' must be a matrix-like object; you can coerce it to one with 'as.matrix()'")
  if (is.null(design) || add.means){
    stopifnot(ncol(object)==length(grp))
    if (check.names){ stopifnot(colnames(object)==names(grp)) }
  }
  if (any(duplicated(rownames(object)))) stop("object cannot have duplicated rownames.")
  if (any(rownames(object)=="")) stop("object cannot have an empty rowname ''.")
  
  if (is.null(design)){
    design <- stats::model.matrix(~0+grp)
    colnames(design) <- sub("grp", "", colnames(design), fixed=TRUE)
  }
  
  if (length(weights)!=1 || !is.na(weights)){
    if (!is.matrix(object) && !is.null(object$weights)){ warning("object$weights are being ignored") }
    fit <- limma::lmFit(object, design=design, block = block, correlation = correlation, weights=weights)
  } else {
    fit <- limma::lmFit(object, design=design, block = block, correlation = correlation)
  }
  
  contr.mat <- limma::makeContrasts(contrasts=contrast.v, levels=design)
  fit2 <- limma::contrasts.fit(fit, contr.mat)
  
  fit2 <- ezebayes(fit2, moderated = TRUE, trend=trend)

  # limma ignores names of contrast.v when it's given as vector
  if (!is.null(names(contrast.v))){
    stopifnot(colnames(fit2$contrasts)==contrast.v)
    colnames(fit2$contrasts) <- names(contrast.v)
  }
  # limma loses rowname if object has only one row
  if (nrow(object) == 1) rownames(fit2$coefficients) <- rownames(object)
  
  # tost
  stopifnot(ncol(fit2$coefficients)==1)
  acoef <- abs(fit2$coefficients)
  se <- as.matrix(fit2$stdev.unscaled) * sqrt(fit2$s2.post)
  tstat.right <- (acoef - tost.lfc)/se
  # equivalent to ptost=max(p1, p2) from TOSTER::TOSTtwo.raw & identical expression from equivalence::tost
  pv <- pt(tstat.right, df=fit2$df.total, lower.tail=TRUE)
  qv <- p.adjust(pv, method=adjust.method)
  mtt <- cbind(fit2$coefficients, tstat.right, pv, qv)
  colnames(mtt) <- paste0(names(contrast.v), "_tost.", c("logFC", "t", "p", "FDR"))
  mtt <- mtt[order(mtt[,2]),, drop=FALSE]
  
  # cbind grp means
  if (add.means){
    if (nrow(object) > 1){
      grp.means <- t(apply(object, MARGIN=1, FUN=function(v) tapply(v, INDEX=grp, FUN=mean, na.rm=TRUE)))
    } else {
      grp.means <- data.matrix(t(tapply(object[1,], INDEX=grp, FUN=mean, na.rm=TRUE)))
      # rownames gets lost when only one row, so grab it from input object
      rownames(grp.means) <- rownames(object)
    }
    colnames(grp.means) <- paste(colnames(grp.means), "avg", sep=".")
    mtt <- cbind(grp.means[rownames(mtt),,drop=FALSE], mtt)
  }
  return(mtt)
}