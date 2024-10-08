% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_cor.R
\name{multi_cor}
\alias{multi_cor}
\title{Test correlation of each row of an object to each column of pheno.tab}
\usage{
multi_cor(
  object,
  pheno.tab,
  method = c("pearson", "spearman", "kendall", "limma"),
  reorder.rows = TRUE,
  prefix = NULL,
  block = NULL,
  correlation = NULL,
  adjust.method = "BH",
  covariates = NULL,
  check.names = TRUE,
  limma.cols = c("AveExpr", "P.Value", "adj.P.Val", "logFC")
)
}
\arguments{
\item{object}{Matrix-like data object containing log-ratios or log-expression values, with rows corresponding to 
features (e.g. genes) and columns to samples. Must have row names that are non-duplicated and non-empty.}

\item{pheno.tab}{Matrix-like data object with columns as sample phenotypes, with \code{nrow(pheno.tab)==ncol(object)}.}

\item{method}{Character string indicating which association is to be used 
for the test. One of \code{"pearson"}, \code{"spearman"}, \code{"kendall"}, 
from \code{\link[stats]{cor.test}} or \code{"limma"} for \code{\link{limma_cor}}.}

\item{reorder.rows}{Logical, should rows be reordered by p-value?}

\item{prefix}{Character string to add to beginning of column names. \code{NULL} does not add a prefix.}

\item{block}{Vector specifying a blocking variable on the samples. Has length = \code{ncol(object)}. 
Must be \code{NULL} if \code{ndups > 1}.}

\item{correlation}{Numeric vector of inter-duplicate or inter-technical replicate correlations. Must be given if 
\code{!is.null(block)}. Its length should be the same as the number of columns of \code{pheno.tab}.}

\item{adjust.method}{Method used to adjust the p-values for multiple testing. Options, in increasing conservatism, 
include \code{"none"}, \code{"BH"}, \code{"BY"}, and \code{"holm"}. See \code{\link[stats]{p.adjust}} for the complete
list of options. A \code{NULL} value will result in the default adjustment method, which is \code{"BH"}.}

\item{covariates}{If \code{method="limma"}, numeric vector or matrix of covariates to include in 
\code{\link{limma_cor}} \code{design} matrix.}

\item{check.names}{Logical; should \code{rownames(pheno.tab)=colnames(object)} be checked?}

\item{limma.cols}{If \code{method="limma"}, \code{cols} from \code{\link{limma_cor}} to include.}
}
\value{
Data frame with several statistical columns corresponding to each phenotype and one row per feature.
}
\description{
Test correlation of each row of an object to each column of \code{pheno.tab} using 
one of Pearson's, Kendall's, or Spearman's correlation methods, or limma 
regression in \code{\link{limma_cor}}. See examples in vignette.
}
\details{
Each column of \code{pheno.tab} is tested independently. Arguments \code{covariates}, \code{block}, and
\code{correlation} only apply if \code{method="limma"}. When each individual \code{pheno.tab} column is tested, 
if some samples have \code{NA}s for that column, those samples are omitted for that column only.
}
