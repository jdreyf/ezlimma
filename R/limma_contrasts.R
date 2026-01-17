#' Apply limma's lmFit, contrasts.fit, & eBayes to one or more contrasts and return a table
#' 
#' Apply \pkg{limma}'s \code{lmFit}, \code{contrasts.fit}, & \code{eBayes} to one or more contrasts, and return
#' a table. See examples in vignette.
#' 
#' @param object Matrix-like data object containing log-ratios or log-expression values, with rows corresponding to 
#' features (e.g. genes) and columns to samples. Must have row names that are non-duplicated and non-empty.
#' @param grp Vector of sample groups. These must be valid variable names in R and the same length as 
#' \code{ncol(object)}.
#' @param contrast.v Named vector of contrasts, passed to \code{\link[limma]{makeContrasts}}.
#' @param design Design matrix of the experiment, with rows corresponding to samples and columns to coefficients to be 
#' estimated.
#' @param weights Non-negative observation weights. Can be a numeric matrix of individual weights of same size as the 
#' \code{object}, or a numeric vector of sample weights with length \code{ncol(object)}, or a numeric vector of gene 
#' weights with length equal to \code{nrow(object)}. Set to \code{NULL} to ignore \code{object$weights}. 
#' \code{weights=NA} (with length one) doesn't pass weights to \code{limma}.
#' @param trend Logical; should an intensity-trend be allowed for the prior variance? Default is that the prior variance 
#' is constant.
#' @param ndups Positive integer giving the number of times each distinct probe is measured in each sample. See Details.
#' @param spacing Positive integer giving the spacing between duplicate occurrences of the same probe, with 
#' \code{spacing=1} for consecutive rows.  See Details.
#' @param block Vector specifying a blocking variable on the samples. Has length = \code{ncol(object)}. 
#' Must be \code{NULL} if \code{ndups > 1}.
#' @param correlation Inter-duplicate or inter-technical replicate correlation. Must be given if 
#' \code{ndups>1} or \code{!is.null(block)}.
#' @param adjust.method Method used to adjust the p-values for multiple testing. Options, in increasing conservatism, 
#' include \code{"none"}, \code{"BH"}, \code{"BY"}, and \code{"holm"}. See \code{\link[stats]{p.adjust}} for the complete
#' list of options. A \code{NULL} value will result in the default adjustment method, which is \code{"BH"}.
#' @param add.means Logical indicating if (unweighted) group means per row should be added to the output.
#' @param treat.lfc Vector of logFC passed to \code{\link[limma:lmFit]{treat}} \code{lfc}. It is recycled as needed to match
#' rows of \code{object}. If given, \code{length(contrast.v)} must be 1. McCarthy & Smyth suggest a 10 percent fold-change,
#' which is \code{treat.lfc=log2(1.1)}. Limma Treat tests for logFC's outside the interval [-treat.lfc, treat.lfc].
#' @param moderated Logical; should \code{\link[limma]{eBayes}} be used? Otherwise an unmoderated version for 
#' \pkg{limma} to produce ordinary least squares statistics is used.
#' @param check.names Logical; should \code{names(grp)==rownames(object)} be checked? Ignored if \code{is.null(design)}
#' and \code{add.means} is \code{FALSE}.
#' @param cols Columns of \code{topTable} output to include. Possibilities include 
#' \code{"logFC", "SE", "AveExpr", "z", "t", "P.Value", "adj.P.Val", "CI.L", "CI.R", "B"}. Some of these column names are then changed here. 
#' If \code{logFC} is specified, \code{FC} will automatically also be given.
#' @return Data frame.
#' @details If \code{design} is \code{NULL} and \code{grp} is given, design will be calculated as 
#' \code{model.matrix(~0+grp)}. However, \code{grp} isn't needed if \code{design} is provided & \code{add.means} 
#' is \code{FALSE}.
#' 
#' The defaults of arguments \code{ndups} and \code{spacing} are set to \code{NULL}, 
#' which allows these arguments to be overridden by the elements \code{object$printer$ndups} and
#' \code{object$printer$spacing}, respectively, if these exist. Whereas, if an element does not exist,
#' the corresponding argument is treated as being its default in \code{\link[limma]{lmFit}}, 
#' i.e. \code{ndups=1} or \code{spacing=1}. 
#' If either of these arguments are specified, they would override the 
#' respective element of \code{object$printer}, if the element existed.
#' 
#' When \code{moderated} is FALSE, an error is generated if \code{!is.null(treat.lfc)} or \code{trend} is TRUE.
#' 
#' Cols \code{CI.L} and \code{CI.R} give the 95% confidence interval of the \code{logFC}. \code{SE} gives the empirical Bayes standard error.
#'   
#' @references McCarthy DJ & Smyth GK (2009). Testing significance relative to a fold-change threshold is a TREAT. 
#' Bioinformatics 25, 765-771.
#' @seealso \code{\link[limma]{lmFit}}; \code{\link[limma]{eBayes}}; \code{\link{limma_cor}}.
#' @export

# don't include parameters for robust fitting, since ppl unlikely to use
# need to set ndups=NULL, st users can set ndups=1 to override object$printer$ndups
limma_contrasts <- function(object, grp=NULL, contrast.v, design=NULL, weights=NA, trend=FALSE, ndups=NULL, spacing=NULL,
                            block=NULL, correlation=NULL, adjust.method="BH", add.means=!is.null(grp), treat.lfc=NULL, 
                            moderated=TRUE, check.names=TRUE, cols=c("P.Value", "adj.P.Val", "logFC")){
  
  # spacing could be given if ndups=1, since lmFit might use ndups=object$printer$ndups
  stopifnot(!is.null(dim(object)), !is.null(rownames(object)), !is.null(colnames(object)), ncol(object)>1,
            length(weights)!=1 || is.na(weights), length(weights)<=1 || 
              (is.numeric(weights) && all(weights>=0) && !all(is.na(weights))), 
            length(weights)<=1 || all(dim(weights)==dim(object)) || 
              length(weights)==nrow(object) || length(weights)==ncol(object), 
            is.null(ndups) || (is.numeric(ndups) && ndups >= 1), 
            is.null(spacing) || (is.numeric(spacing) && spacing>=1), 
            is.null(correlation) || (is.numeric(correlation) && abs(correlation)<=1),
            is.null(ndups) || ndups < 2 || is.null(block),
            is.null(treat.lfc) || length(contrast.v)==1, moderated || !trend)
  
  if ((!is.null(block) || (!is.null(ndups) && ndups > 1)) && is.null(correlation))
    stop("!is.null(block) or ndups>1, so correlation must not be NULL.")
  
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
  
  # lmFit tests if is.null(block), rather than missing(block), so can safely provide block=NULL
  # lmFit requires !missing(correlation) if !(ndups < 2 && is.null(block))
  args.lst <- list(object=object, design=design, block = block, correlation = correlation)
  
  # 22 Nov 2020: limma 3.47.1
  #- Explicitly setting `weights=NULL` in a call to lmFit() no longer over-rides the `weights` value found in `object`.
  # - Default settings in lmFit() changed from `ndups=` and `spacing=1` to `ndups=NULL` and `spacing=NULL`. 
  # No change to function behavior from a user point of view.
  
  # handles weights=NULL or weights as vector or matrix
  if (length(weights)!=1 || !is.na(weights)) args.lst <- c(args.lst, list(weights=weights))
  if (!is.null(ndups)) args.lst <- c(args.lst, ndups=ndups)
  if (!is.null(spacing)) args.lst <- c(args.lst, spacing=spacing)
  
  # use do.call so that can specify args as list()
  # fit <- limma::lmFit(object, design=design, block = block, correlation = correlation)
  fit <- do.call(limma::lmFit, args = args.lst)
  
  # if ndups>1, lmFit returns y$genes=uniquegenelist(y$probes,ndups=ndups,spacing=spacing), but does not fix row names
  # assign row names before multitoptab, since it reorders rows 
  if (!is.null(ndups) && ndups >= 2){
    if (is.null(spacing)) spacing <- 1
    rownames(fit$coefficients) <- limma::uniquegenelist(rownames(object), ndups=ndups, spacing=spacing)
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
  
  # limma until 3.66.0 ignores names of contrast.v when it's given as vector
  if (!is.null(names(contrast.v)) && any(colnames(fit2$contrasts) != names(contrast.v))){
    stopifnot(colnames(fit2$contrasts)==contrast.v)
    colnames(fit2$contrasts) <- names(contrast.v)
  }
  
  mtt <- multiTopTab(fit2, cols=cols, adjust.method=adjust.method)
  
  # cbind grp means
  if (add.means){
    # if only one row, need to fix rownames(mtt), also
    if (nrow(object)==1){
      grp.means <- as.matrix(t(tapply(object[1,], INDEX=grp, FUN=mean, na.rm=TRUE)))
      rownames(grp.means) <- rownames(mtt) <- rownames(object)[1]
    } else{
      grp.means <- grp_means(object=object, grp=grp, args.lst = args.lst)
    }
    mtt <- cbind(grp.means[rownames(mtt),,drop=FALSE], mtt)
  }
  return(mtt)
}