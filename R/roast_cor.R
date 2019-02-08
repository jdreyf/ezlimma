# Test association of gene sets to phenotype using rotation testing with output
#' to Excel
#'
#' Test association of gene sets to phenotype using rotation testing with
#' \code{\link[limma]{roast}} using \code{limma} functions \code{mroast} or \code{fry}.
#' It returns a data frame with statistics per gene set, and writes this to an
#' Excel file. The Excel file links to CSV files, which contain statistics per
#' gene set.
#'
#' @inheritParams limma_contrasts
#' @inheritParams limma_cor
#' @inheritParams roast_contrasts
#' @return Data frame of gene set statistics.
#' @details Pathway (i.e. gene set) names are altered to be valid filenames in Windows and Linux. Numeric columns 
#' are rounded to 3 significant figures.
#' @export

roast_cor <- function(object, G, feat.tab, name=NA, phenotype = NULL, design = NULL, 
                    fun=c("fry", "mroast"), set.statistic = "mean",
                    weights = NA, gene.weights=NULL, trend = FALSE, block = NULL, 
                    correlation = NULL, prefix=NULL, adjust.method = "BH", min.nfeats=3, max.nfeats=1000, 
                    alternative=c("two.sided", "less", "greater"), nrot=999, check_names=TRUE, seed=0){
  
  stopifnot(rownames(object) %in% rownames(feat.tab), !is.null(design)|!is.null(phenotype),
            is.null(gene.weights)|length(gene.weights)==nrow(object), ncol(object) > 1)
  #only mroast takes some arguments
  if (fun=="fry" && (set.statistic!="mean" || !is.null(gene.weights) || adjust.method!="BH")){
    warning("fry method does not take arguments: set.statistic, gene.weights, or adjust.method. These arguments will be ignored.")
  }
  
  if (!is.null(phenotype)){
    stopifnot(length(phenotype)==ncol(object), limma::isNumeric(phenotype))
    if (check_names){
      stopifnot(names(phenotype)==colnames(object))
    }
  }
  fun <- match.arg(fun)
  alternative <- match.arg(alternative)
  
  if (fun=="mroast") set.seed(seed=seed)
  
  ##get G index
  index <- g_index(G=G, object=object, min.nfeats=min.nfeats, max.nfeats=max.nfeats)
  
  if (is.null(design)){
      #model.matrix clips NAs in phenotype, so need to also remove from object, weights
      n.na <- sum(is.na(phenotype))
      if (n.na > 0){
          warning(n.na, " NAs in phenotype removed")
          pheno.nona <- phenotype[!is.na(phenotype)]
          #len > 1 excludes NA and NULL
          if(length(weights) > 1){
            if (length(weights)==ncol(object)){ 
              weights <- weights[!is.na(phenotype)]
            } else {
              weights <- weights[,!is.na(phenotype)]
            }
          } #end if !is.null weights
          object <- object[,!is.na(phenotype), drop=FALSE]
      } else {
          pheno.nona <- phenotype
      } #end if/else n.na > 0
      design <- stats::model.matrix(~1+pheno.nona) 
      colnames(design) <- gsub("\\(|\\)", "", gsub("pheno.nona", "", colnames(design), fixed=TRUE))
  } else {
    if (!is.null(phenotype)) warning("Phenotype is ignored, since design is given.")
  }
  stopifnot(is.numeric(design[,2]))
  
  #deal with weights
  if (!is.matrix(object)){
      if (!is.null(object$weights)){
          #NULL has len = 0, given weights should have len > 1 -> only weights=NA has len 1
          if (length(weights) != 1){
            warning("object$weights are being ignored")
          } else {
              #weights = NA
              if (is.null(dimnames(object$weights))) dimnames(object$weights) <- dimnames(object)
              weights <- object$weights
          }
      }
  }#end if(!is.matrix(object))
  
  ##run for contrast = 2
  #mroast vs fry x is.na(weights) vs not, so 4 conditions
  if (fun=="fry"){
    if (length(weights) == 1 && is.na(weights)){
      res <- limma::fry(y = object, index = index, design = design, contrast = 2, trend = trend, block = block, correlation = correlation)
    } else {
      res <- limma::fry(y = object, index = index, design = design, contrast = 2,
                        weights = weights, trend = trend, block = block, correlation = correlation)
    }
  } else {
    if (length(weights) == 1 && is.na(weights)){
      res <- limma::mroast(y = object, index = index, design = design, contrast = 2,
                           set.statistic = set.statistic, gene.weights = gene.weights,  
                           trend = trend, block = block, correlation = correlation,
                           adjust.method = adjust.method, nrot = nrot)
    } else {
      res <- limma::mroast(y = object, index = index, design = design, contrast = 2,
                           set.statistic = set.statistic, weights = weights, gene.weights = gene.weights,  
                           trend = trend, block = block, correlation = correlation,
                           adjust.method = adjust.method, nrot = nrot)
    }
  }
  #need to coerce Direction from factor to char
  res$Direction <- as.character(res$Direction)
  
  #if want one-sided test, change p-values, calc new FDRs, then remove Mixed columns
  if (alternative!="two.sided"){
    res <- roast_one_tailed(roast.res=res, fun=fun, alternative=alternative, nrot=nrot, adjust.method=adjust.method)
  }
  colnames(res) <- gsub("PValue", "p", gsub("FDR.Mixed", "Mixed.FDR", gsub("PValue.Mixed", "Mixed.p", colnames(res))))
  #add prefix to each column except 1st, which is NGenes
  if (!is.null(prefix)) colnames(res)[-1] <- paste(prefix, colnames(res)[-1], sep = ".")
  #let combine_pvalues find pvalue columns
  res <- res[order(combine_pvalues(res)), ]
  
  #change FDR to appropriate adjustment name if user doesn"t use FDR
  if (!(adjust.method %in% c("BH", "fdr"))){
    colnames(res) <- gsub("FDR$", adjust.method, colnames(res))
  }
  
  res.xl <- df_signif(res, digits = 3)
  #write xlsx file with links
  if (!is.na(name)){
    nm <- paste(name, fun, sep="_")
    write_top_xl(pwy.tab=res.xl, feat.lst=index, feat.tab=feat.tab, name=nm)
  }
  return(res)
}#end fcn