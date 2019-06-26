# ezlimma
R package that streamlines & extends the popular R/Bioconductor package `limma`.

[![Build Status](https://travis-ci.org/jdreyf/ezlimma.svg?branch=master)](https://travis-ci.org/jdreyf/ezlimma)
[![Coverage Status](https://img.shields.io/codecov/c/github/jdreyf/ezlimma/master.svg)](https://codecov.io/github/jdreyf/ezlimma?branch=master)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

## Install
On Windows, you should have [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

Install `ezlimma` from GitHub using `remotes` within R. You must install `remotes`, e.g. with `install.packages("remotes", repos="https://cloud.r-project.org")`, if you haven't before. `ezlimma` depends on `limma` so you must also install this using instruction below if you haven't before.
```
#if haven't already installed limma
install.packages("BiocManager", repos="https://cloud.r-project.org") #if haven't already installed BiocManager
library(BiocManager)
BiocManager::install("limma")

library(remotes)
remotes::install_github(repo="jdreyf/ezlimma", build_opts = c("--no-resave-data", "--no-manual"))
```

## Usage
The vignette presents a tutorial. To see the vignette:
```
library(limma)
library(ezlimma)
library(rmarkdown)
browseVignettes(package="ezlimma")
```
and click on HTML.
