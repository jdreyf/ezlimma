#' Wrapper for limma's eBayes function that also allows for computing OLS statistics
#' 
#' Wrapper for \pkg{limma}'s \code{\link[limma]{eBayes}}, which computes empirical Bayes statistics, that allows for 
#' using the rest of the \pkg{limma} workflow in \pkg{ezlimma} to compute ordinary least squares statistics.
#' 
#' @param fit Object of class \code{MArrayLM} as produced by \code{\link[limma]{lmFit}} alone (without 
#' \code{\link[limma]{eBayes}}).
#' @param moderated Logical; should \code{\link[limma]{eBayes}} be used (otherwise an unmoderated version for \pkg{limma} to 
#' produce ordinary least squares statistics is used)?
#' @inheritParams limma_contrasts
#' @return Object of \code{\link[limma]{MArrayLM-class}}.
#' @details \code{trend} is only applicable if \code{moderated} is \code{TRUE}.

ezebayes <- function(fit, moderated=TRUE, trend=TRUE){
  stopifnot(is.logical(moderated), is.logical(trend))
  if (moderated){
    fit <- limma::eBayes(fit, trend=trend)
  } else {
    # coefficients <- fit$coefficients
    # stdev.unscaled <- fit$stdev.unscaled
    # sigma <- fit$sigma
    # df.residual <- fit$df.residual
    if (is.null(fit$coefficients) || is.null(fit$stdev.unscaled) || is.null(fit$sigma) || is.null(fit$df.residual)){
      stop("No data, or argument is not a valid lmFit object")
    }
    if (all(fit$df.residual == 0)) stop("No residual degrees of freedom in linear model fits")
    if (all(!is.finite(fit$sigma))) stop("No finite residual standard deviations")
    
    # limma user guide, pg 85: ordinary.t <- fit$coef / fit$stdev.unscaled / fit$sigma
    # https://support.bioconductor.org/p/113833/ with Gordon Smyth: 
    # fit2$t <- fit2$coef/fit2$stdev.unscaled/fit2$sigma
    # fit2$p.value <- 2 * pt(-abs(fit2$t), df = fit2$df.residual)
    fit$t <- fit$coefficients/fit$stdev.unscaled/fit$sigma
    fit$p.value <- 2 * pt(-abs(fit$t), df = fit$df.residual)

    if (!is.null(fit$design) && is.fullrank(fit$design)){
      F.stat <- limma::classifyTestsF(fit, fstat.only = TRUE)
      fit$F <- as.vector(F.stat)
      df1 <- attr(F.stat, "df1")
      df2 <- attr(F.stat, "df2")
      if (df2[1] > 1e+06){
        fit$F.p.value <- stats::pchisq(df1 * fit$F, df1, lower.tail = FALSE)
      } else {
        fit$F.p.value <- stats::pf(fit$F, df1, df2, lower.tail = FALSE)
      }
    } # end !is.null(fit$design) && is.fullrank(fit$design) 
  } # end else
  fit
}