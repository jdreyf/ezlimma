#' Apply limma's lmFit, contrasts.fit, & eBayes to one or more contrasts and return a table
#' 
#' Apply \pkg{limma}'s \code{lmFit}, \code{contrasts.fit}, & \code{eBayes} to one or more contrasts, and return
#' a table. See examples in vignette.
#' 
#' @param object Matrix-like data object containing log-ratios or log-expression values, with rows corresponding to 
#' features (e.g. genes) and columns to samples. Must have rownames that are non-duplicated and non-empty.
#' @param grp Vector of sample groups. These must be valid variable names in R and the same length as 
#' \code{ncol(object)}.
#' @param contrast.v Named vector of contrasts, passed to \code{\link[limma]{makeContrasts}}.
#' @param design Design matrix of the experiment, with rows corresponding to samples and columns to coefficients to be 
#' estimated.
#' @param weights Non-negative observation weights. Can be a numeric matrix of individual weights of same size as the 
#' \code{object}, or a numeric vector of sample weights with length \code{ncol(object)}, or a numeric vector of gene 
#' weights with length equal to\code{nrow(object)}. Set to \code{NULL} to ignore \code{object$weights}. 
#' \code{weights=NA} (with length one) doesn't pass weights to \code{limma}.
#' @param trend Logical; should an intensity-trend be allowed for the prior variance? Default is that the prior variance 
#' is constant.
#' @param block Vector specifying a blocking variable on the samples. Has length = \code{ncol(object)}.
#' @param correlation Inter-duplicate or inter-technical replicate correlation.
#' @param adjust.method Method used to adjust the p-values for multiple testing. Options, in increasing conservatism, 
#' include \code{"none"}, \code{"BH"}, \code{"BY"}, and \code{"holm"}. See \code{\link[stats]{p.adjust}} for the complete
#' list of options. A \code{NULL} value will result in the default adjustment method, which is \code{"BH"}.
#' @param add.means Logical indicating if (unweighted) group means per row should be added to the output.
#' @param treat.lfc Vector of logFC passed to \code{\link[limma:lmFit]{treat}} \code{lfc}. It is recycled as needed to match
#' rows of \code{object}. If given, \code{length(contrast.v)} must be 1. McCarthy & Smyth suggest a 10% fold-change,
#' which is \code{treat.lfc=log2(1.1)}.
#' @param moderated Logical; should \code{\link[limma]{eBayes}} be used? Otherwise an unmoderated version for 
#' \pkg{limma} to produce ordinary least squares statistics is used.
#' @param check.names Logical; should \code{names(grp)==rownames(object)} be checked? Ignored if \code{is.null(design)}
#' and \code{add.means} is \code{FALSE}.
#' @param cols Columns of \code{topTable} output to include. Possibilities include 
#' \code{"logFC", "AveExpr", "z", "t", "P.Value", "adj.P.Val", "B"}. Some of these column names are then changed here. 
#' If \code{logFC} is specified, \code{FC} will automatically also be given.
#' @return Data frame.
#' @details If \code{design} is \code{NULL} and \code{grp} is given, design will be calculated as 
#' \code{model.matrix(~0+grp)}. However, \code{grp} isn't needed if \code{design} is provided & \code{add.means} 
#' is \code{FALSE}.
#' 
#' When \code{moderated} is FALSE, an error is generated if \code{!is.null(treat.lfc)} or \code{trend} is TRUE.  
#' @references McCarthy DJ & Smyth GK (2009). Testing significance relative to a fold-change threshold is a TREAT. 
#' Bioinformatics 25, 765-771.
#' @seealso \code{\link[limma]{lmFit}}; \code{\link[limma]{eBayes}}; \code{\link{limma_cor}}.
#' @export

# don't include parameters for robust fitting, since ppl unlikely to use
limma_contrasts <- function(object, grp=NULL, contrast.v, design=NULL, weights=NA, trend=FALSE, block=NULL, 
                            correlation=NULL, adjust.method="BH", add.means=!is.null(grp), treat.lfc=NULL, 
                            moderated=TRUE, check.names=TRUE, cols=c("P.Value", "adj.P.Val", "logFC")){
  
  stopifnot(all(is.na(weights)) || is.null(weights) || dim(weights)==dim(object) || length(weights)==nrow(object) || 
            length(weights)==ncol(object), is.null(treat.lfc) || length(contrast.v)==1, moderated || !trend)
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
  
  # can't set weights=NULL in lmFit when using voom, since lmFit only assigns
  # weights "if (missing(weights) && !is.null(y$weights))"
  # can't make this into separate function, since then !missing(weights)
  # length(NULL)=0; other weights should have length > 1
  if (length(weights)!=1 || !is.na(weights)){
    if (!is.matrix(object) && !is.null(object$weights)){ warning("object$weights are being ignored") }
    fit <- limma::lmFit(object, design=design, block = block, correlation = correlation, weights=weights)
  } else {
    fit <- limma::lmFit(object, design=design, block = block, correlation = correlation)
  }
  
  contr.mat <- limma::makeContrasts(contrasts=contrast.v, levels=design)
  fit2 <- limma::contrasts.fit(fit, contr.mat)
  
  if (is.null(treat.lfc)){
    if (trend && !moderated) stop("'trend' must be FALSE when 'moderated' is FALSE.")
    fit2 <- ezebayes(fit2, moderated = moderated, trend=trend)
  } else {
    if (!moderated) stop("'treat.lfc' must be NULL when 'moderated' is FALSE.")
    fit2 <- limma::treat(fit2, lfc=treat.lfc, trend=trend)
  }
  
  # limma ignores names of contrast.v when it's given as vector
  if (!is.null(names(contrast.v))){
    stopifnot(colnames(fit2$contrasts)==contrast.v)
    colnames(fit2$contrasts) <- names(contrast.v)
  }
  
  mtt <- multiTopTab(fit2, cols=cols, adjust.method=adjust.method)
  
  # cbind grp means
  if (add.means){
    if (nrow(object) > 1){
      grp.means <- t(apply(object, MARGIN=1, FUN=function(v) tapply(v, INDEX=grp, FUN=mean, na.rm=TRUE)))
    } else {
      grp.means <- data.matrix(t(tapply(object[1,], INDEX=grp, FUN=mean, na.rm=TRUE)))
      rownames(grp.means) <- rownames(mtt)
    }
    colnames(grp.means) <- paste(colnames(grp.means), "avg", sep=".")
    mtt <- cbind(grp.means[rownames(mtt),,drop=FALSE], mtt)
  }
  return(mtt)
}