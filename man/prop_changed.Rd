% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prop_changed.R
\name{prop_changed}
\alias{prop_changed}
\title{Calculate proportion of changed features per pathway}
\usage{
prop_changed(feat.tab, feat.lst, stat.col = 1, p.col = 2)
}
\arguments{
\item{feat.tab}{Table of feature (e.g. gene) statistics.}

\item{feat.lst}{Named list, with each element containing names or indices of rows of \code{feat.tab}.}

\item{stat.col}{Column name or column index with signed statistics.}

\item{p.col}{Column name or column index with p-values.}
}
\description{
Calculate proportion of features up & down at p=0.05 per pathway.
}
