#'Apply \code{signif} function to numeric columns of data frames
#'
#'Apply \code{signif} function to numeric columns of data frames
#'
#'@param df A data frame
#'@param digits number of digits passed to \code{signif}.
#'@export

#to improve can look at
#https://stackoverflow.com/questions/9063889/how-to-round-a-data-frame-in-r-that-contains-some-character-variables
df_signif <- function(df, digits=3){
  stopifnot(is.data.frame(df))
  cols.num <- sapply(df, is.numeric)
  df.num <- df[,cols.num, drop=FALSE]
  df.num <- signif(df.num, digits=digits)
  df.notnum <- df[,!cols.num, drop=FALSE]
  df2 <- data.frame(df.num, df.notnum, check.names=FALSE)
  #revert column names to original
  df2 <- df2[,colnames(df)]
  return(df2)
}