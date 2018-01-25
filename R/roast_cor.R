#'Test association of gene sets to phenotype using rotation testing with output
#'to Excel
#'
#'Test association of gene sets to phenotype using rotation testing with
#'\code{\link[limma]{roast}} using \code{limma} functions \code{mroast} or \code{fry}.
#'It returns a data frame with statistics per gene set, and writes this to an
#'Excel file. The Excel file links to CSV files, which contain statistics per
#'gene set.
#'
#'@param object A matrix-like data object containing log-ratios or 
#'  log-expression values for a series of arrays, with rows corresponding to 
#'  genes and columns to samples.
#'@param G a gene set list as returned from \link{read_gmt}.
#'@param stats.tab a table of feature (e.g. gene) statistics that the Excel 
#'  table can link to.
#'@param name a name for the folder and Excel file that get written. Set to \code{NA} 
#'  to avoid writing output.
#'@param phenotype Vector of phenotypes of the samples. Should be same length as
#'  \code{ncol(object)}. If the vector is named, names should match 
#'  \code{colnames(object)}.
#'@param design the design matrix of the experiment, with rows corresponding to 
#'  arrays and columns to coefficients to be estimated. Can be used to provide 
#'  covariates.
#'@param fun function to use, either \code{fry} or \code{mroast}.
#'@param set.statistic summary set statistic, if using \code{mroast}. 
#'  Possibilities are \code{"mean"},\code{"floormean"}, \code{"mean50"}, or 
#'  \code{"msq"}.
#'@param weights non-negative observation weights passed to \code{lmFit}. Can be
#'  a numeric matrix of individual weights, of same size as the object 
#'  expression matrix, or a numeric vector of array weights with length equal to
#'  \code{ncol} of the expression matrix, or a numeric vector of gene weights 
#'  with length equal to \code{nrow} of the expression matrix.
#'@param gene.weights numeric vector of directional (positive or negative) 
#'  genewise weights, passed to \code{mroast} or \code{fry}. This vector must 
#'  have length equal to \code{nrow(object)}.
#'@param trend logical, should an intensity-trend be allowed for the prior 
#'  variance? Default is that the prior variance is constant.
#'@param block vector or factor specifying a blocking variable on the arrays. 
#'  Has length equal to the number of arrays.
#'@param correlation the inter-duplicate or inter-technical replicate 
#'  correlation.
#'@param adjust.method method used to adjust the p-values for multiple testing.
#'@param min.ngenes minimum number of genes needed in a gene set for testing.
#'@param max.ngenes maximum number of genes needed in a gene set for testing.
#'@param alternative indicates the alternative hypothesis and must be one of 
#'  \code{"two.sided"}, \code{"greater"} or \code{"less"}. \code{"greater"} 
#'  corresponds to positive association, \code{"less"} to negative association.
#'@param n.toptabs number of gene set toptables to write to CSV and link to from
#'  Excel
#'@param nrot number of rotations used to compute the p-values in \code{mroast}.
#'@return data frame of gene set statistics.
#'@export

roast_cor <- function(object, G, stats.tab, name=NA, phenotype = NULL, design = NULL, 
                    fun=c("fry", "mroast"), set.statistic = 'mean',
                    weights = NULL, gene.weights=NULL, trend = FALSE, block = NULL, 
                    correlation = NULL, adjust.method = 'BH', min.ngenes=3, max.ngenes=1000, 
                    alternative=c("two.sided", "less", "greater"), n.toptabs = Inf, nrot=999){
  stopifnot(rownames(object) %in% rownames(stats.tab), !is.null(design)|!is.null(phenotype),
            is.null(gene.weights)|length(gene.weights)==nrow(object))
  if (!is.null(phenotype)){
    stopifnot(length(phenotype)==ncol(object), is.numeric(phenotype), 
              names(phenotype)==colnames(object))
  }
  fun <- match.arg(fun)
  alternative <- match.arg(alternative)
  
  ##get G index
  index <- g_index(G=G, object=object, min.ngenes=min.ngenes, max.ngenes=max.ngenes)
  
  if (is.null(design)){
      #model.matrix clips NAs in phenotype, so need to also remove from object, weights
      n.na <- sum(is.na(phenotype))
      if (n.na > 0){
          warning(n.na, ' NAs in phenotype removed')
          pheno.nona <- phenotype[!is.na(phenotype)]
          object <- object[,!is.na(phenotype), drop=FALSE]
          if(!is.null(weights)){
            if (length(weights)==ncol(object)){ 
              weights <- weights[!is.na(phenotype)]
            } else {
              weights <- weights[,!is.na(phenotype)]
            }
          }#end if !is.null weights
      } else {
          pheno.nona <- phenotype
      }
      design <- model.matrix(~1+pheno.nona) 
      colnames(design) <- gsub('\\(|\\)', '', gsub('pheno.nona', '', colnames(design), fixed=TRUE))
  }
  stopifnot(colnames(design)[1] == 'Intercept' & is.numeric(design[,2]))
    
  if (!is.matrix(object)){
      if (!is.null(object$weights)){
          if (!is.null(weights)){
            warning('object$weights are being ignored') 
          } else {
              if (is.null(dimnames(object$weights))) dimnames(object$weights) <- dimnames(object)
              weights <- object$weights
          }
      }
  }#end if(!is.matrix(object))
  
  ##run for contrast = 2
  if (fun=="fry"){
    res <- fry(y = object, index = index, design = design, contrast = 2,
               weights = weights, gene.weights = gene.weights, trend.var = trend, 
               block = block, correlation = correlation, adjust.method = adjust.method)
  } else {
    res <- mroast(y = object, index = index, design = design, contrast = 2,
                  set.statistic = set.statistic, weights = weights, gene.weights = gene.weights,  
                  trend.var = trend, block = block, correlation = correlation,
                  adjust.method = adjust.method, nrot = nrot)
  }
  
  ##if want one-sided test, change p-values, calc new FDRs, then remove Mixed columns
  if (alternative!="two.sided"){
    res <- roast_one_tailed(roast.res=res, fun=fun, alternative=alternative, nrot=nrot, adjust.method=adjust.method)
  }
  colnames(res) <- gsub("PValue", "p", gsub("FDR.Mixed", "Mixed.FDR", gsub("PValue.Mixed", "Mixed.p", colnames(res))))
  #add prefix to each column except 1st, which is NGenes
  #but colnames(design)[2] is ""
  #colnames(res)[-1] <- paste(colnames(design)[2], colnames(res)[-1], sep = '.')
  res <- res[order(combine_pvalues(res, grep('\\.p$', colnames(res)))), ]
  
  #change FDR to appropriate adjustment name if user doesn't use FDR
  if (!(adjust.method %in% c("BH", "fdr"))){
    colnames(res) <- gsub("\\.FDR$", adjust.method, colnames(res))
  }
  
  ##write xlsx file with links
  if (!is.na(name)){
    write_linked_xlsx(name=name, fun=fun, res=res, index=index, stats.tab=stats.tab, n.toptabs=n.toptabs)
  }
  return(res)
}#end fcn