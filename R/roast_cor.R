#'Test association of gene sets to phenotype using rotation testing with output
#'to Excel
#'
#'Test association of gene sets to phenotype using rotation testing with
#'\code{limma} function \code{\link[limma]{fry}} or \code{\link[limma]{mroast}}.
#'It returns a dataframe with statistics per gene set, and writes this to an
#'Excel file. The Excel file links to CSV files, which contain statistics per
#'gene set.
#'
#'@param object A matrix-like data object containing log-ratios or 
#'  log-expression values for a series of arrays, with rows corresponding to 
#'  genes and columns to samples.
#'@param G a gene set list returned from \code{\link{read.gmt}}.
#'@param stats.tab a table of feature (e.g. gene) statistics, that the Excel 
#'  table can link to.
#'@param name a name for the folder and Excel file that get written.
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
#'@return dataframe of gene set statistics.
#'@export

roast_cor <- function(object, G, stats.tab, name, phenotype = NULL, design = NULL, 
                    fun=c("fry", "mroast"), set.statistic = 'mean',
                    weights = NULL, gene.weights=NULL, trend = FALSE, block = NULL, 
                    correlation = NULL, adjust.method = 'BH', min.ngenes=3, max.ngenes=1000, 
                    alternative=c("two.sided", "less", "greater"), n.toptabs = Inf){
  stopifnot(rownames(object) %in% rownames(stats.tab), length(gene.weights)==nrow(object),
            !is.null(design)|!is.null(phenotype))
  if (!is.null(phenotype)){
    stopifnot(length(phenotype)==ncol(object), is.numeric(phenotype), 
              names(phenotype)==colnames(object))
  }
  fun <- match.arg(fun)
  alternative <- match.arg(alternative)
  
  ##define objects
  index <- lapply(G, function(g) rownames(object)[rownames(object) %in% g$genes])
  names(index) <- sapply(G, function(g) g$name)
  ##remove gene sets of the wrong size
  index <- index[sapply(index, function(x) length(x) >= min.ngenes & length(x) <= max.ngenes)]
  
  if (is.null(design)){
      #model.matrix clips NAs in phenotype, so need to also remove from object, weights
      n.na <- sum(is.na(phenotype))
      if (n.na > 0){
          warning(n.na, 'NAs in phenotype removed')
          pheno.nona <- phenotype[!is.na(phenotype)]
          object <- object[,!is.na(phenotype), drop=FALSE]
          if(!is.null(weights)){
            if (length(weights)==ncol(object)){ 
              weights <- weights[!is.na(pheno)]
            } else {
              weights <- weights[,!is.na(pheno)]
            }  
          }#end if !is.null weights
      } else {
          pheno.nona <- pheno
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
  #if want one-sided test, change p-values, calc new FDRs, then remove Mixed columns
  if (alternative!="two.sided"){ 
    direction <- sub("greater", "Up", sub("less", "Down", alternative))
    if (fun=="fry"){
      res[,'PValue'] <- fry_two2one_tailed(res, direction = direction)
    } else {
      res[,'PValue'] <- mroast_two2one_tailed(res, direction = one.tailed, nrot = nrot)
    }
    res[,'FDR'] <- p.adjust(res[,'PValue'], method = adjust.method)
    mixed.cols <- grep('Mixed', colnames(res))
    if (length(mixed.cols) > 0){ res <- res[,-mixed.cols] }
  }#end if one.tailed
  
  colnames(res) <- gsub("PValue", "p", gsub("FDR.Mixed", "Mixed.FDR", gsub("PValue.Mixed", "Mixed.p", colnames(res))))
  #add prefix to each column except 1st, which is NGenes
  colnames(res)[-1] <- paste(colnames(design)[2], colnames(res)[-1], sep = '.')
  res <- res[order(combine.pvalues(res, grep('\\.p$', colnames(res)))), ]
  
  #change FDR to appropriate adjustment name if user doesn't use FDR
  if (!(adjust.method %in% c("BH", "fdr"))){
    colnames(res) <- gsub("\\.FDR$", adjust.method, colnames(res))
  }
  
  ##write xlsx file with links
  name <- paste(name, fun, sep='_')
  dir.create(name)
  dir.create(paste0(name, '/pathways'))
  
  if (n.top > nrow(res)) n.top <- nrow(res)
  pwys <- rownames(res)[1:n.top]
  for(pwy in pwys){
      stat <- stats.tab[index[[pwy]], ]
      stat <- stat[order(combine.pvalues(stat)), ]
      write.csv(stat, paste0(name, '/pathways/', substr(pwy, 1, 150), '.csv'))
  }
  
  wb <- xlsx::createWorkbook()
  sheet <- xlsx::createSheet(wb, sheetName = name)
  xlsx::addDataFrame(df.signif(res, 3), sheet)
  rows  <- xlsx::getRows(sheet)
  cells <- xlsx::getCells(rows, 1)
  
  for(i in seq_along(pwys)){
    xlsx::addHyperlink(cells[[paste0(i+1, '.1')]], paste0('pathways/', substr(pwys[i], 1, 150), '.csv'), 'FILE')
  }
  xlsx::saveWorkbook(wb, paste0(name, '/', name, '.xlsx'))
  return(res)
}#end fcn