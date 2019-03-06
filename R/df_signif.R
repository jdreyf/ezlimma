#' Apply \code{signif} function to numeric columns of data frames
#' 
#' Apply \code{signif} function to numeric columns of data frames.
#' 
#' @param tab A data frame
#' @param digits number of digits passed to \code{signif}.
#' @export

# to improve can look at
# https://stackoverflow.com/questions/9063889/how-to-round-a-data-frame-in-r-that-contains-some-character-variables
df_signif <- function(tab, digits=3){
  stopifnot(is.data.frame(tab))
  cols.num <- sapply(tab, is.numeric)
  tab.num <- tab[,cols.num, drop=FALSE]
  tab.num <- signif(tab.num, digits=digits)
  tab.notnum <- tab[,!cols.num, drop=FALSE]
  tab2 <- data.frame(tab.num, tab.notnum, check.names=FALSE)
  # revert column names to original
  tab2 <- tab2[,colnames(tab)]
  return(tab2)
}