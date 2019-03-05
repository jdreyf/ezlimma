#' Read Excel XLSX file with links to CSVs
#' 
#' Read Excel XLSX file with links to CSVs written by \code{write_linked_xl}.
#' 
#' @param file.path Path to the XLSX file.
#' @return Data frame.
#' @export

read_linked_xl <- function(file.path){
  if (!requireNamespace("readxl", quietly = TRUE)){
    stop("Install 'readxl' package.", call. = FALSE)
  }
  if (!requireNamespace("tidyxl", quietly = TRUE)){
    stop("Install 'tidyxl' package.", call. = FALSE)
  }
  
  df1 <- suppressMessages(readxl::read_xlsx(path=file.path))[,-1]
  df1 <- as.data.frame(df1)
  
  df2 <- tidyxl::xlsx_cells(path=file.path)
  df2 <- df2[df2$col == 1,]
  pwy.nms <- gsub(".+,|\\)", "", gsub('\"', "", df2$formula))
  rownames(df1) <- pwy.nms
  return(df1)
}