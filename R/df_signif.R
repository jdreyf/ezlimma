#'Apply \code{signif} function to numeric columns of data frames
#'
#'Apply \code{signif} function to numeric columns of data frames
#'
#'@param df A dataframe
#'@param digits number of digits to take significance to

#don't export
#to improve can look at
#https://stackoverflow.com/questions/9063889/how-to-round-a-data-frame-in-r-that-contains-some-character-variables
df_signif <- function(df, digits=3){
  cols.num <- sapply(df, is.numeric)
  mat <- df[,cols.num]
  mat <- signif(mat, digits=digits)
  nonmat <- data.frame(df[,!cols.num], check.names=FALSE)
  if (sum(!cols.num)==1){ colnames(nonmat) <- colnames(df)[!cols.num] }
  df2 <- data.frame(mat, nonmat, check.names=FALSE)
  #revert column names to original
  df2 <- df2[,colnames(df)]
  return(df2)
}