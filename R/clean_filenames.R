#' Alter names to valid filenames
#'
#' Alter names to be valid filenames in both Linux and Windows, and ensure they are unique. Improper characters 
#' are replaced with "_". Filenames are limited to 199 characters for Windows.
#'
#' @param nm Character string of names.
#' @export

clean_filenames <- function(nm){
  # https://stackoverflow.com/questions/1976007/what-characters-are-forbidden-in-windows-and-linux-directory-names
  # safe characters
  nm <- gsub("[^[:alnum:]\\.]", "_", nm)
  # can't end in " " or "."
  nm <- gsub("\\.$", "_", nm)
  unsafe <- c("CON", "PRN", "AUX", "NUL",
              "COM1", "COM2", "COM3", "COM4", "COM5", "COM6", "COM7", "COM8", "COM9",
              "LPT1", "LPT2", "LPT3", "LPT4", "LPT5", "LPT6", "LPT7", "LPT8", "LPT9")
  unsafe.ind <- which(toupper(nm) %in% unsafe)
  if (length(unsafe.ind) > 0){
    nm[unsafe.ind] <- paste0(nm[unsafe.ind], "_")
  }
  
  while (any(duplicated(nm))){
    dup.ind <- which(duplicated(nm))
    nm[dup.ind] <- paste0(nm[dup.ind], "_")
  }
  
  # limit filenames to 199 characters
  # https://stackoverflow.com/questions/265769/maximum-filename-length-in-ntfs-windows-xp-and-windows-vista
  nm <- substr(x=nm, start=1, stop=199)
  
  return(nm)
}