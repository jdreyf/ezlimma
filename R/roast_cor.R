# Test association of gene sets to phenotype using rotation testing with output
#' to Excel
#'
#' Test association of gene sets to phenotype using rotation testing with
#' \code{\link[limma]{roast}} using \pkg{limma} functions \code{mroast} or \code{fry}.
#' It returns a data frame with statistics per gene set, and writes this to an
#' Excel file. The Excel file links to CSV files, which contain statistics per
#' gene set. See example in vignette of \code{\link[ezlimma]{roast_contrasts}}.
#'
#' @inheritParams limma_contrasts
#' @inheritParams limma_cor
#' @inheritParams roast_contrasts
#' @return Data frame of gene set statistics.
#' @inherit roast_contrasts details
#' @seealso \code{\link[ezlimma]{roast_contrasts}}.
#' @export

roast_cor <- function(object, G, feat.tab=NULL, name=NA, phenotype = NULL, design = NULL,
                    fun=c("fry", "mroast"), set.statistic = "mean", weights = NA, gene.weights=NULL,
                    trend = FALSE, block = NULL, correlation = NULL, prefix=NULL, adjust.method = "BH",
                    min.nfeats=3, max.nfeats=1000, alternative=c("two.sided", "less", "greater"),
                    nrot=999, check.names=TRUE, pwy.nchar=199, seed=0){

  fun <- match.arg(fun)
  alternative <- match.arg(alternative)
  # get G index
  index <- g_index(G=G, object=object, min.nfeats=min.nfeats, max.nfeats=max.nfeats)
  .roast_cor_index(object=object, index=index, feat.tab=feat.tab, name=name, phenotype=phenotype, design=design,
                   fun=fun, set.statistic=set.statistic, weights=weights, gene.weights=gene.weights,
                   trend=trend, block=block, correlation=correlation, prefix=prefix, adjust.method=adjust.method,
                   alternative=alternative, nrot=nrot, check.names=check.names, pwy.nchar=pwy.nchar, seed=seed)
}

#' Core of \code{roast_cor}, taking a pre-computed gene set \code{index}
#'
#' Not exported. Split out of \code{\link{roast_cor}} so that \code{\link{roast_multi_cor}} can resolve gene set
#' membership once (via \code{\link{g_index}}) and reuse it, and pre-resolved integer indices (via
#' \code{index.int}), across all phenotype columns, rather than redoing that \code{nrow(object) * length(G)}-scale
#' work on every column.
#'
#' @param index Gene set list as returned by \code{\link{g_index}}, i.e. character-vector probe/gene IDs per set.
#' @param index.int As \code{index}, but with sets already resolved to integer positions into
#' \code{rownames(object)}. If \code{NULL} (the default), this is derived from \code{index}. Passing it in lets a
#' caller that loops over many calls with the same \code{object} rows (e.g. \code{roast_multi_cor}) do that
#' resolution once instead of on every call.
#' @inheritParams roast_cor
#' @noRd

.roast_cor_index <- function(object, index, index.int=NULL, feat.tab=NULL, name=NA, phenotype = NULL, design = NULL,
                    fun="fry", set.statistic = "mean", weights = NA, gene.weights=NULL,
                    trend = FALSE, block = NULL, correlation = NULL, prefix=NULL, adjust.method = "BH",
                    alternative="two.sided", nrot=999, check.names=TRUE, pwy.nchar=199, seed=0){

  stopifnot(!is.null(dim(object)), !is.null(rownames(object)), !is.null(colnames(object)), ncol(object) > 1,
            !is.null(design)|!is.null(phenotype),
            length(weights)!=1 || is.na(weights), length(weights)<=1 ||
              (is.numeric(weights) && all(weights>=0) && !all(is.na(weights))),
            length(weights)<=1 || all(dim(weights)==dim(object)) ||
              length(weights)==nrow(object) || length(weights)==ncol(object),
            is.null(gene.weights) || length(gene.weights)==nrow(object),
            is.na(name) || all(rownames(object) %in% rownames(feat.tab)))

  if (!is.null(block) && is.null(correlation))
    stop("!is.null(block), so correlation must not be NULL.")

  # only mroast takes some arguments
  if (fun=="fry" && (!is.null(gene.weights) || set.statistic!="mean")){
    warning("fry method does not take the argument gene.weights or set.statistic, so these will be ignored.")
  }
  if (fun=="fry" && alternative == "two.sided" && adjust.method!="BH"){
    warning("When alternative is 'two.sided', fry method does not take the argument adjust.method, so it will be ignored.")
  }

  if (!is.null(phenotype)){
    stopifnot(length(phenotype)==ncol(object), limma::isNumeric(phenotype))
    if (check.names){
      stopifnot(names(phenotype)==colnames(object))
    }
  }

  if (fun=="mroast") set.seed(seed=seed)

  # resolve gene set membership to integer indices into rownames(object), unless already done by the caller.
  # limma::fry()/mroast() re-derive integer indices from character IDs (which(geneid %in% iset)) on every call,
  # so reusing pre-resolved integer indices across many calls (e.g. once per phenotype in roast_multi_cor)
  # avoids repeating that O(nrow(object)) matching, per gene set, every time.
  if (is.null(index.int)) index.int <- lapply(index, function(x) if (is.character(x)) match(x, rownames(object)) else x)

  if (is.null(design)){
      # model.matrix clips NAs in phenotype, so need to also remove from object, weights
      n.na <- sum(is.na(phenotype))
      if (n.na > 0){
          message(n.na, " NAs in phenotype removed")
          pheno.nona <- phenotype[!is.na(phenotype)]
          # len > 1 excludes NA and NULL
          if(length(weights) > 1){
            if (length(weights)==ncol(object)){ 
              weights <- weights[!is.na(phenotype)]
            } else {
              weights <- weights[,!is.na(phenotype), drop=FALSE]
            }
          } # end if !is.null weights
          object <- object[,!is.na(phenotype), drop=FALSE]
      } else {
          pheno.nona <- phenotype
      } # end if/else n.na > 0
      design <- stats::model.matrix(~1+pheno.nona) 
      colnames(design) <- gsub("\\(|\\)", "", gsub("pheno.nona", "", colnames(design), fixed=TRUE))
  } else {
    if (!is.null(phenotype)) warning("Phenotype is ignored, since design is given.")
  }
  stopifnot(is.numeric(design[,2]))
  
  # deal with weights
  if (!is.matrix(object)){
      if (!is.null(object$weights)){
          # NULL has len = 0, given weights should have len > 1 -> only weights=NA has len 1
          if (length(weights) != 1){
            warning("object$weights are being ignored")
          } else {
              # weights = NA
              if (is.null(dimnames(object$weights))) dimnames(object$weights) <- dimnames(object)
              weights <- object$weights
          }
      }
  }# end if(!is.matrix(object))
  
  # run for contrast = 2
  # mroast vs fry x is.na(weights) vs not, so 4 conditions
  if (fun=="fry"){
    if (length(weights) == 1 && is.na(weights)){
      res <- limma::fry(y = object, index = index.int, design = design, contrast = 2, trend = trend, block = block,
                        correlation = correlation)
    } else {
      res <- limma::fry(y = object, index = index.int, design = design, contrast = 2,
                        weights = weights, trend = trend, block = block, correlation = correlation)
    }
  } else {
    if (length(weights) == 1 && is.na(weights)){
      res <- limma::mroast(y = object, index = index.int, design = design, contrast = 2,
                           set.statistic = set.statistic, gene.weights = gene.weights,
                           trend = trend, block = block, correlation = correlation,
                           adjust.method = adjust.method, nrot = nrot)
    } else {
      res <- limma::mroast(y = object, index = index.int, design = design, contrast = 2,
                           set.statistic = set.statistic, weights = weights, gene.weights = gene.weights,
                           trend = trend, block = block, correlation = correlation,
                           adjust.method = adjust.method, nrot = nrot)
    }
  }
  # need to coerce Direction from factor to char
  res$Direction <- as.character(res$Direction)
  
  # if want one-sided test, change p-values, calc new FDRs, then remove Mixed columns
  if (alternative!="two.sided"){
    res <- roast_two2one_tailed(roast.res=res, fun=fun, alternative=alternative, nrot=nrot, adjust.method=adjust.method)
  }
  colnames(res) <- gsub("PValue", "p", 
                        gsub("FDR.Mixed", "Mixed.FDR", 
                             gsub("PValue.Mixed", "Mixed.p", colnames(res))))
  # add prefix to each column except 1st, which is NGenes
  if (!is.null(prefix)) colnames(res)[-1] <- paste(prefix, colnames(res)[-1], sep = ".")
  # let combine_pvalues find pvalue columns
  res <- res[order(combine_pvalues(res)), ]
  
  # change FDR to appropriate adjustment name if user doesn"t use FDR
  if (!(adjust.method %in% c("BH", "fdr"))){
    colnames(res) <- gsub("FDR$", adjust.method, colnames(res))
  }
  
  res.xl <- df_signif(res, digits = 8)
  # write xlsx file with links
  if (!is.na(name)){
    nm <- paste(name, fun, sep="_")
    write_linked_xl(pwy.tab=res.xl, feat.lst=index, feat.tab=feat.tab, name=nm, pwy.nchar=pwy.nchar)
  }
  return(res)
}
