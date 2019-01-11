# ezlimma
R package that streamlines & extends the popular R/Bioconductor package `limma`.

[![Build Status](https://travis-ci.org/jdreyf/ezlimma.svg?branch=master)](https://travis-ci.org/jdreyf/ezlimma)
[![Coverage Status](https://img.shields.io/codecov/c/github/jdreyf/ezlimma/master.svg)](https://codecov.io/github/jdreyf/ezlimma?branch=master)

## Install
Install `ezlimma` from GitHub using `remotes` within R. You must install `remotes` if you haven't before. `ezlimma` depends on `limma` so you must also install this if you haven't before.
```
#if haven't already installed limma
source("http://bioconductor.org/biocLite.R")
biocLite("limma")

install.packages("remotes") #if haven't already installed remotes
library(remotes)
remotes::install_github(repo="jdreyf/ezlimma", build_opts = c("--no-resave-data", "--no-manual"))
```

## Usage
The vignette presents a tutorial. To see the vignette:
```
library(limma)
library(ezlimma)
browseVignettes(package="ezlimma")
```
and click on HTML.