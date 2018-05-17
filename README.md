# ezlimma
`ezlimma` is an R/Bioconductor package that streamlines & extends the popular R/Bioconductor package `limma`.

[![Build Status](https://travis-ci.org/jdreyf/ezlimma.svg?branch=master)](https://travis-ci.org/jdreyf/ezlimma)
[![Coverage Status](https://img.shields.io/codecov/c/github/jdreyf/ezlimma/master.svg)](https://codecov.io/github/jdreyf/ezlimma?branch=master)

## Install
Install `ezlimma` from GitHub using `devtools` within R. You must install `devtools` if you haven't before. `ezlimma` depends on `limma` so you must also install this if you haven't before.
```
source("http://bioconductor.org/biocLite.R")
biocLite("limma") #if haven't already installed limma
install.packages("devtools") #if haven't already installed devtools
library("devtools")
install_github("jdreyf/ezlimma")
```

## Usage
See the vignette for a tutorial.

## Acknowledgement
Thanks to Grace Daher (Joslin Diabetes Center) for testing the package and making suggestions.
