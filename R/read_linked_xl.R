#' Read Excel XLSX file with links to CSVs
#' 
#' Read Excel XLSX file with links to CSVs written by \code{write_linked_xl}.
#' 
#' @param file.xlsx Path of XLSX file.
#' @return Data frame.
#' @export

read_linked_xl <- function(file.xlsx){
  if (!requireNamespace("readxl", quietly = TRUE)){
    stop("Install 'readxl' package.", call. = FALSE)
  }
  if (!requireNamespace("tidyxl", quietly = TRUE)){
    stop("Install 'tidyxl' package.", call. = FALSE)
  }
  
  df1 <- suppressMessages(readxl::read_xlsx(path=file.xlsx))[,-1]
  df1 <- as.data.frame(df1)
  
  df2 <- tidyxl::xlsx_cells(path=file.xlsx)
  df2 <- df2[df2$col == 1,]
  pwy.nms <- gsub(".+,|\\)", "", gsub('\"', "", df2$formula))
  rownames(df1) <- pwy.nms
  return(df1)
}