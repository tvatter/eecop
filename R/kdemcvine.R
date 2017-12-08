#' Multi-class kernel density estimator based on simplified C-vine copulas
#'
#' Implements a multi-class C-vine based kernel density estimator.
#'
#' @param x (\eqn{n x (d-1)}) data matrix.
#' @param y (\eqn{n x 1}) data vector for the central node.
#' @param w (\eqn{n x 1}) data vector for the classes.
#' @param mult_1d numeric; all bandwidhts for marginal kernel density estimation
#'   are multiplied with \code{mult_1d}. Defaults to `log(1 + d)`.
#' @param xmin numeric vector of length d; see \code{\link[kdevine]{kde1d}}.
#' @param xmax numeric vector of length d; see \code{\link[kdevine]{kde1d}}.
#' @param ... further arguments passed to \code{\link[kdevine]{kde1d}} or
#'   \code{\link[kdevine]{kdevinecop}}.
#'
#' @return An object of class \code{kdemcvine}.
#' 
#' @author Thibault Vatter
#'
#' @seealso \code{\link{dkdecvine}} \code{\link[kdevine]{kdevine}} \code{\link[kdevine]{kdevine}}
#'
#' @examples
#' TODO
#'
#' @export
kdemcvine <- function(x, y, w, mult_1d = 1, xmin = NULL, xmax = NULL, ...) {
  
  if (is.vector(x)) {
    d <- 2
  } else {
    d <- ncol(x)+1
  }

  res <- lapply(unique(w), function(d) {
    sel <- w == d
    kdecvine(x[sel,], y[sel], mult_1d, xmin, xmax, ...)
  })
  
  class(res) <- "kdemcvine"
  res
}

#' Predictions for multi-class kernel based C-vine model
#' 
#' Use a fitted multi-class C-vine model to predict the treatment effect.
#'
#' @param object a `kdemcvine` object.
#' @param newdata points where the fit shall be evaluated.
#' @param ... 
#'
#' @return A vector of predictions.
#'
#' @author Thibault Vatter
#'
#' @seealso \code{\link{kdemcvine}}
#'
#' @examples
#' TODO
#'
#' @importFrom kdevine dkde1d pkde1d dkdevinecop
#' @importFrom kdecopula dkdecop hkdecop
#' @export
predict.kdemcvine <- function(object, newdata = NULL, ...) {
  
  stopifnot(all(is.element(what, c("mean", "quantile"))))
  what <- unique(what)
  
  if (any(what == "quantile")) {
    if (!is.null(list(...)$tau)) {
      tau <- list(...)$tau
      stopifnot(is.vector(tau) && is.numeric(tau) && all(tau > 0 & tau < 1))
    } else {
      tau <- c(0.1,0.25,0.5,0.75,0.9)
    }
  }
  
  # if (any(what == "causal")) {
  #   stopifnot(!is.null(list(...)$w) && all(list(...)$w == 0 | list(...)$w == 1))
  #   w <- list(...)$w
  # }
  
  d <- length(object$marg.dens)
  if (!is.null(newdata)) {
    stopifnot((is.vector(newdata) && length(newdata) == d-1) ||
                (is.matrix(newdata) && ncol(newdata) == d-1))
    if (is.vector(newdata)) {
      newdata <- matrix(newdata, nrow = 1)
    }
  } else {
    newdata <- sapply(2:d, function(k) object$marg.dens[[k]]$x_cc)
  }
  
  y <- as.numeric(object$marg.dens[[1]]$x_cc)
  uy <- pkde1d(y, object$marg.dens[[1]])
  cxy <- apply(newdata, 1, function(x) getcxy(x, uy, object))
  
  f <- function(type) {
    switch(
      type,
      "mean"     = getMean(y, cxy),
      "quantile" = getQuantile(y, cxy, tau))
  }
  browser()
  
  if (length(what) > 1) {
    res <- lapply(what, f)
  } else {
    res <- f(what)
  }
  
  return(res)
}