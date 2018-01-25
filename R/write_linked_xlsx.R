#' Write Excel xlsx file with links to CSVs
#' 
#' Write Excel xlsx file with links to CSVs
#' 
#' @param name a name for the folder and Excel file that get written.
#' @param fun function to use, either \code{fry} or \code{mroast}.
#' @param res result of \code{mroast} or \code{fry}.
#' @param index index of pathways in G
#' @param stats.tab a table of feature (e.g. gene) statistics that the Excel file can link to
#' @param n.toptabs number of gene set toptables to write to CSV and link to from Excel.

write_linked_xlsx <- function(name, fun, res, index, stats.tab, n.toptabs){
  name <- paste(name, fun, sep="_")
  dir.create(name)
  dir.create(paste0(name, '/pathways'))
  
  if (n.toptabs > nrow(res)) n.toptabs <- nrow(res)
  
  pwys <- rownames(res)[1:n.toptabs]
  for(pwy in pwys){
    stat <- stats.tab[index[[pwy]], ]
    stat <- stat[order(combine_pvalues(stat)), ]
    write.csv(stat, paste0(name, '/pathways/', substr(pwy, 1, 150), '.csv'))
  }
  
  wb <- xlsx::createWorkbook()
  sheet <- xlsx::createSheet(wb, sheetName = name)
  xlsx::addDataFrame(df_signif(res, 3), sheet)
  rows  <- xlsx::getRows(sheet)
  cells <- xlsx::getCells(rows, 1)
  
  for(i in seq_along(pwys)){
    xlsx::addHyperlink(cells[[paste0(i+1, '.1')]], paste0('pathways/', substr(pwys[i], 1, 150), '.csv'), 'FILE')
  }
  xlsx::saveWorkbook(wb, paste0(name, '/', name, '.xlsx'))
  return(TRUE)
}