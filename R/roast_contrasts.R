#' Test contrasts of gene sets between groups using rotation testing with output to Excel
#'
#' Test contrasts of gene sets using \code{\link[limma]{roast}} with functions \code{mroast} or \code{fry}. It 
#' returns a data frame with statistics per gene set, and writes this to an Excel file. The Excel file links to 
#' CSV files, which contain statistics per gene set. Some of the arguments only apply to \code{mroast}. See 
#' examples in vignette.
#' 
#' @param G Gene set list as returned from \code{\link{read_gmt}}.
#' @param feat.tab Table of feature (e.g. gene) statistics that the Excel table can link to.
#' @param fun Function to use, either \code{fry} or \code{mroast}.
#' @param set.statistic Summary set statistic. Possibilities are \code{"mean"}, \code{"floormean"}, \code{"mean50"}, or 
#' \code{"msq"}. Only for \code{mroast}.
#' @param name Name for the folder and Excel file that get written. Set to \code{NA} to avoid writing output.
#' @param gene.weights Numeric vector of directional (positive or negative) genewise weights. These represent each 
#' gene's contribution to pathways. They are not for precision weights, from \code{weights}. This vector must have 
#' length equal to \code{nrow(object)}. Only for \code{mroast}.
#' @param adjust.method Method used to adjust the p-values for multiple testing. Only for \code{mroast}.
#' @param min.nfeats Minimum number of features (e.g. genes) needed in a gene set for testing.
#' @param max.nfeats Maximum number of features (e.g. genes) needed in a gene set for testing.
#' @param nrot Number of rotations used to estimate the p-values for \code{mroast}.
#' @param alternative Alternative hypothesis; must be one of \code{"two.sided"}, \code{"greater"} or \code{"less"}. 
#' \code{"greater"} corresponds to positive association, \code{"less"} to negative association.
#' @param check.names Logical; should \code{names(grp)==rownames(object)} be checked? Ignored if \code{is.null(design)}.
#' @param pwy.nchar Numeric maximum number of characters allowed in pathway name.
#' @param seed Integer seed to set for reproducility if \code{fun="mroast"}, since \code{mroast} uses random 
#' simulations. Ignored if \code{fun="fry"}.
#' @inheritParams limma_contrasts
#' @return Data frame of gene set statistics.
#' @details Pathway (i.e. gene set) names are altered to be valid filenames in Windows and Linux. Numeric columns are
#' rounded to 3 significant figures.
#' 
#' This function does not accept \code{ndups} nor \code{spacing}, because we found that these arguments do not impact
#' \code{limma::mroast} nor \code{limma::fry}.
#' @export

# limma 3.34.6 fixed array weights bug, but I don't require this version of limma, since don't have it on server
roast_contrasts <- function(object, G, feat.tab, grp=NULL, contrast.v, design=NULL, fun=c("fry", "mroast"), 
                            set.statistic = "mean", name=NA, weights = NA, gene.weights = NULL, trend = FALSE, 
                            block = NULL, correlation = NULL, adjust.method = "BH", min.nfeats=3, max.nfeats=1000, nrot=999,
                            alternative=c("two.sided", "less", "greater"), check.names=TRUE, pwy.nchar=199, seed=0){
  
  stopifnot(!is.null(dim(object)), !is.null(rownames(object)), !is.null(colnames(object)), ncol(object) > 1,
    rownames(object) %in% rownames(feat.tab), !is.null(design) || !is.null(grp),
            length(weights)!=1 || is.na(weights), length(weights)<=1 || 
              (is.numeric(weights) && all(weights>=0) && !all(is.na(weights))), 
            length(weights)<=1 || dim(weights)==dim(object) || 
              length(weights)==nrow(object) || length(weights)==ncol(object),
            is.null(gene.weights) || length(gene.weights)==nrow(object))
  
  if (!is.null(block) && is.null(correlation))
    stop("!is.null(block), so correlation must not be NULL.")
  
  # only mroast takes some arguments
  if (fun=="fry" && (!is.null(gene.weights)||adjust.method!="BH")){
    warning("fry method does not take arguments: gene.weights or adjust.method. These arguments will be ignored.")
  }
  fun <- match.arg(fun)
  alternative <- match.arg(alternative)
  
  if (fun=="mroast") set.seed(seed=seed)
  
  # get G index
  index <- g_index(G=G, object=object, min.nfeats=min.nfeats, max.nfeats=max.nfeats)
  
  if (is.null(design)){
    stopifnot(ncol(object) == length(grp))
    design <- stats::model.matrix(~0+grp)
    colnames(design) <- sub("grp", "", colnames(design), fixed=TRUE)
    if (check.names){ stopifnot(colnames(object) == names(grp)) }
  }
  
  contr.mat <- limma::makeContrasts(contrasts = contrast.v, levels = design)
  args.lst <- list(y=object, index=index, design=design, block = block, correlation = correlation, trend=trend)
  
  # deal with weights
  if (!is.matrix(object)){
    if (!is.null(object$weights)){
      if (length(weights) != 1){
        warning("object$weights are being ignored")
        args.lst <- c(args.lst, list(weights=weights))
      } else {
        if (is.null(dimnames(object$weights))) dimnames(object$weights) <- dimnames(object)
        args.lst <- c(args.lst, list(weights=object$weights))
      }
    }
  } else if (length(weights) != 1) args.lst <- c(args.lst, list(weights=weights))
  
  # run fry or mroast for each contrast
  # block & correlation from lmFit, trend from eBayes
  for (i in seq_along(contrast.v)){
    args.tmp <- c(args.lst, list(contrast = contr.mat[, i]))
    if (fun=="fry"){
      res.tmp <- do.call(limma::fry, args = args.tmp)
    } else {
      args.tmp <- c(args.tmp, list(gene.weights = gene.weights, set.statistic = set.statistic, nrot=nrot))
      res.tmp <- do.call(limma::mroast, args = args.tmp)

      # PropUp & PropDown don't use feat.tab & threshold at z=sqrt(2) -> p=0.1
      res.tmp <- res.tmp[, setdiff(colnames(res.tmp), c("PropDown", "PropUp"))]
    }# end roast
    # need to coerce "direction" from factor to char
    res.tmp$Direction <- as.character(res.tmp$Direction)
    
    # add propup & propdown from feat.tab if same contr
    contr.nm <- names(contrast.v)[i]
    contr.cols <- paste(contr.nm, c("logFC", "p"), sep=".")
    if (all(contr.cols %in% colnames(feat.tab))){
      prop.diff <- signif(prop_changed(feat.tab=feat.tab[, contr.cols, drop=FALSE], feat.lst=index), 2)
      res.tmp <- data.frame(res.tmp[, 1:2, drop=FALSE], prop.diff[rownames(res.tmp),, drop=FALSE], 
                            res.tmp[, -(1:2), drop=FALSE])
    }
    
    # if want one-sided test, change p-values, calc new FDRs, then remove Mixed columns
    if (alternative!="two.sided"){
      res.tmp <- roast_two2one_tailed(roast.res=res.tmp, fun=fun, alternative=alternative, 
                                      nrot=nrot, adjust.method=adjust.method)
    }# end if one.tailed
    colnames(res.tmp) <- gsub("PValue", "p", 
                              gsub("FDR.Mixed", "Mixed.FDR", 
                                   gsub("PValue.Mixed", "Mixed.p", colnames(res.tmp))))
    # add contrast names to each column except 1st, which is NGenes
    colnames(res.tmp)[-1] <- paste(names(contrast.v[i]), colnames(res.tmp)[-1], sep = ".")
    if (i == 1){
      res <- res.tmp
    } else {
      res <- cbind(res, res.tmp[rownames(res), -1])
    }
  }# end for i
  
  # let combine_pvalues find pvalue columns
  if (nrow(res) > 1) res <- res[order(combine_pvalues(res)), ]
  
  # change FDR to appropriate adjustment name if user doesn"t use FDR
  if (!(adjust.method %in% c("BH", "fdr"))){
    colnames(res) <- gsub("FDR$", adjust.method, colnames(res))
  }
  
  res.xl <- df_signif(res, digits = 3)
  # write xlsx file with links
  if (!is.na(name)){
    nm <- paste(name, fun, sep="_")
    write_linked_xl(pwy.tab=res.xl, feat.lst=index, feat.tab=feat.tab, name=nm, pwy.nchar=pwy.nchar)
  }# end if !is.na(name)
  return(res)
}