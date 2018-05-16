#' Write Excel XLSX file with links to CSVs
#' 
#' Write Excel XLSX file using package \code{xlsx} with links to CSVs
#' 
#' @param name a name for the folder and Excel file that get written.
#' @param fun function to use, either \code{fry} or \code{mroast}.
#' @param res result of \code{mroast} or \code{fry}.
#' @param index index of pathways in G
#' @param stats.tab a table of feature (e.g. gene) statistics that the Excel file can link to
#' @param n.toptabs number of gene set toptables to write to CSV and link to from Excel.
#' @details This function requires package \code{xlsx}. However, loading it may not work automatically, as explained in
#' this Stack Overflow \href{https://stackoverflow.com/questions/43738366/r-importing-xlsx-package-to-my-own-packag-doesnt-work}{thread}, 
#' in which case you will get an error instructing you to call \code{library(xlsx)}. This function is not meant to be 
#' called directly by the user.

# don't @import xlsx, since don't want it to be installed with ezlimma
write_linked_xlsx <- function(name, fun, res, index, stats.tab, n.toptabs){
  #https://stackoverflow.com/questions/43738366/r-importing-xlsx-package-to-my-own-package-doesnt-work 
  if (!requireNamespace("xlsx", quietly = TRUE)){
    stop("Package 'xlsx' needed for this function to work. Please install it.", call. = FALSE)
  }
  
  name <- paste(name, fun, sep="_")
  dir.create(name)
  dir.create(paste0(name, '/pathways'))
  
  if (n.toptabs > nrow(res)) n.toptabs <- nrow(res)
  
  pwys <- rownames(res)[1:n.toptabs]
  for(pwy in pwys){
    stat <- stats.tab[index[[pwy]], ]
    stat <- stat[order(combine_pvalues(stat)), ]
    utils::write.csv(stat, paste0(name, '/pathways/', substr(pwy, 1, 150), '.csv'))
  }
  
  wb <- xlsx::createWorkbook()
  sheet <- xlsx::createSheet(wb, sheetName = name)
  #below throws "Error in .jfindClass(as.character(class)) : class not found" if xlsx not loaded & attached
  try.xlsx <- try(xlsx::addDataFrame(x=df_signif(res, 3), sheet=sheet))
  if (class(try.xlsx)=="try-error"){
    stop("You need to load and attach 'xlsx' via 'library(xlsx)'.")
  }
  rows  <- xlsx::getRows(sheet)
  cells <- xlsx::getCells(rows, 1)
  
  for(i in seq_along(pwys)){
    xlsx::addHyperlink(cells[[paste0(i+1, '.1')]], paste0('pathways/', substr(pwys[i], 1, 150), '.csv'), 'FILE')
  }
  xlsx::saveWorkbook(wb, paste0(name, '/', name, '.xlsx'))
  return(TRUE)
}