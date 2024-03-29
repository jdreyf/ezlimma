% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_fisher_enrichment.R
\name{multi_fisher_enrichment}
\alias{multi_fisher_enrichment}
\title{Fisher Exact Test of gene set enrichment with multiple vectors of significant genes with output to Excel}
\usage{
multi_fisher_enrichment(
  sig.sets,
  G,
  feat.tab,
  name = NA,
  adjust.method = "BH",
  min.nfeats = 3,
  max.nfeats = 1000,
  pwy.nchar = 199
)
}
\arguments{
\item{sig.sets}{Named list whose elements are vectors of significant gene IDs matching 
\code{rownames(feat.tab)}.}

\item{G}{Gene set list as returned from \code{\link{read_gmt}}.}

\item{feat.tab}{Table of feature (e.g. gene) statistics that the Excel table can link to.}

\item{name}{Name for the folder and Excel file that get written. Set to \code{NA} to avoid writing output.}

\item{adjust.method}{Method used to adjust the p-values for multiple testing. Only for \code{mroast}.}

\item{min.nfeats}{Minimum number of features (e.g. genes) needed in a gene set for testing.}

\item{max.nfeats}{Maximum number of features (e.g. genes) needed in a gene set for testing.}

\item{pwy.nchar}{Numeric maximum number of characters allowed in pathway name.}
}
\value{
List with two elements: \code{pwy.stats}, a table of pathway statistics; and \code{feat.tab}, a table that 
appends a binary matrix of which genes are in which \code{sig.set} and the input \code{feat.tab}.
}
\description{
Test enrichment of multiple vectors of significant genes against the universe of genes in \code{feat.tab} per 
gene set using \code{\link[stats]{fisher.test}}. It returns a data frame with statistics per gene set, and can 
write this to Excel. The Excel file links to CSV files, which contain statistics per genes in a set. If you only
have one vector of significant genes, use \code{\link[ezlimma]{fisher_enrichment}}.
}
\details{
Pathway (i.e. gene set) names are altered to be valid filenames in Windows and Linux. Numeric columns are
rounded to 8 significant figures.
}
\seealso{
\code{\link[ezlimma]{fisher_enrichment}}
}
