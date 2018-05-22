#' Alter names to valid filenames
#'
#' Alter names to be valid Uniform Resource Identifier (URIs) for rJava, valid filenames in both Linux and Windows, 
#' and ensure they are unique. Improper characters are replaced with "_".
#'
#' @param nm Character string of names.

# https://stackoverflow.com/questions/1976007/what-characters-are-forbidden-in-windows-and-linux-directory-names.
rescue_filenames <- function(nm){
  #safe characters
  nm <- gsub("[^[:alnum:]\\.]", "_", nm)
  #can't end in " " or "."
  nm <- gsub("\\.$", "_", nm)
  unsafe <- c('CON', 'PRN', 'AUX', 'NUL',
              'COM1', 'COM2', 'COM3', 'COM4', 'COM5', 'COM6', 'COM7', 'COM8', 'COM9',
              'LPT1', 'LPT2', 'LPT3', 'LPT4', 'LPT5', 'LPT6', 'LPT7', 'LPT8', 'LPT9')
  unsafe.ind <- which(toupper(nm) %in% unsafe)
  if (length(unsafe.ind) > 0){
    nm[unsafe.ind] <- paste0(nm[unsafe.ind], "_")
  }
  
  while (any(duplicated(nm))){
    dup.ind <- which(duplicated(nm))
    nm[dup.ind] <- paste0(nm[dup.ind], "_")
  }
  
  return(nm)
}