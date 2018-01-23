#' EE-vine copula
#'
#' TODO
#'
#' @param x (\eqn{n x (d-1)}) data matrix.
#' @param y (\eqn{n x 1}) data vector for the central node.
#' @param ... further arguments passed to [kde1d()] or [vinecop()].
#'
#' @return An object of class `eevine`.
#' @examples
#' set.seed(0)
#' 
#' ## The predictors
#' p <- 10
#' n <- 1e3
#' X <- matrix(2 * runif(n * p) - 1, n, p)
#' colnames(X) <- paste0('X', 1:p)
#' 
#' ## mean and variance shifts
#' Y1 <- rnorm(n) + exp(X[, 1])
#' Y2 <- rnorm(n) * (1 + exp(X[, 1]))
#' 
#' ## For the evaluation of the results
#' ticks = 1e2
#' X.test = matrix(0, ticks, p)
#' X.test[, 1] = seq(-0.9, 0.9, length.out = ticks)
#' colnames(X.test) <- paste0('X', 1:p)
#' 
#' ## The truth
#' truth1 <- cbind(qnorm(0.1) + exp(X.test[, 1]),
#'                 exp(X.test[,1]),
#'                 qnorm(0.9) + exp(X.test[, 1]))
#' truth2 <- cbind(qnorm(0.1) * (1 + exp(X.test[, 1])),
#'                 0,
#'                 qnorm(0.9) * (1 + exp(X.test[, 1])))
#' 
#' ## range for the margin estimates
#' xmin <- c(-Inf,rep(-1, p))
#' xmax <- c(Inf,rep(1, p))
#' 
#' ## Can we identify the right structure?
#' rvr1 <- eevine(x = X, y = Y1, xmin = xmin, xmax = xmax)
#' rvr2 <- eevine(x = X, y = Y2, xmin = xmin, xmax = xmax)
#' 
#' ## Quantile regression results at levels 10%, 50% and 90%
#' rq1 <- predict(rvr1, X.test, 'quantile', tau = c(0.1, 0.5, 0.9))
#' rq2 <- predict(rvr2, X.test, 'quantile', tau = c(0.1, 0.5, 0.9))
#' 
#' ## Plot the results
#' par(mfrow = c(1,2))
#' plot(X[,1], Y1, xlab = 'X1', ylab = 'Y1')
#' matlines(X.test[,1], truth1, col = 'green', lty = c(2, 1, 2), lwd = 2)
#' matlines(X.test[,1], rq1, col = 'red', lty = c(2, 1, 2), lwd = 2)
#' legend('topleft', legend = c('data', 'truth', 'estimate'),
#'        col = c(1, 'green', 'red'),
#'        pch = c(1, rep(NA, 2)),
#'        lty = c(0, rep(1, 2)),
#'        cex = 0.7)
#' plot(X[,1], Y2, xlab = 'X1', ylab = 'Y2')
#' matlines(X.test[,1], truth2, col = 'green', lty = c(2, 1, 2), lwd = 2)
#' matlines(X.test[,1], rq2, col = 'red', lty = c(2, 1, 2), lwd = 2)
#'
#' @importFrom assertthat assert_that is.number is.count is.scalar
#' @importFrom kdevine kde1d pkde1d
#' @importFrom rvinecopulib bicop vinecop hbicop dbicop
#' @export
eevine <- function(x, y, ...) {
  
  ## sanity checks
  assert_that(is.vector(y) && is.numeric(y))
  assert_that((is.vector(x) || is.matrix(x)) && is.numeric(x))
  if (is.vector(x)) {
    d <- 2
    assert_that(length(y) == length(x))
  } else {
    d <- ncol(x) + 1
    assert_that(length(y) == nrow(x))
  }
  
  data <- cbind(y, x)
  n <- nrow(data)
  
  ## default arguments
  args <- list(...)
  if (is.null(args[["mult_1d"]])) 
    args[["mult_1d"]] <- 1
  if (is.null(args[["xmin"]])) 
    args[["xmin"]] <- rep(-Inf, d)
  if (is.null(args[["xmax"]])) 
    args[["xmax"]] <- rep(Inf, d)
  if (is.null(args[["family_set"]])) 
    args[["family_set"]] <- c("indep", "tll")
  if (is.null(args[["par_method"]])) 
    args[["par_method"]] <- "mle"
  if (is.null(args[["nonpar_method"]])) 
    args[["nonpar_method"]] <- "quadratic"
  if (is.null(args[["mult"]])) 
    args[["mult"]] <- rep(1, d - 1)
  if (is.null(args[["selcrit"]])) 
    args[["selcrit"]] <- "aic"
  if (is.null(args[["presel"]])) 
    args[["presel"]] <- TRUE
  if (is.null(args[["trunc_lvl"]])) {
    args[["trunc_lvl"]] <- Inf
  } else {
    assert_that(is.scalar(args[["trunc_lvl"]]) && 
                  (is.na(args[["trunc_lvl"]]) || is.count(args[["trunc_lvl"]])))
  }
  if (is.null(args[["tree_crit"]])) 
    args[["tree_crit"]] <- "hoeffd"
  if (is.null(args[["threshold"]])) 
    args[["threshold"]] <- 0
  if (is.null(args[["keep_data"]])) 
    args[["keep_data"]] <- FALSE
  if (is.null(args[["show_trace"]])) 
    args[["show_trace"]] <- FALSE
  if (is.null(args[["cores"]])) 
    args[["cores"]] <- 1
  
  ## estimation of the margins
  margin_models <- as.list(numeric(d))
  for (k in 1:d) {
    margin_models[[k]] <- kde1d(data[, k], 
                                xmin = args[["xmin"]][k], 
                                xmax = args[["xmax"]][k], 
                                bw = args[["bw"]][k], 
                                mult = args[["mult_1d"]])
  }
  res <- list(margin_models = margin_models)
  
  ## estimation of the first tree
  udata <- sapply(1:d, function(k) pkde1d(data[, k], margin_models[[k]]))
  res$first_tree <- lapply(2:d, function(k) {
    bicop(udata[, c(1, k)], 
          family_set = args[["family_set"]], 
          par_method = args[["par_method"]], 
          nonpar_method = args[["nonpar_method"]], 
          mult = args[["mult"]], 
          selcrit = args[["selcrit"]], 
          presel = args[["presel"]], 
          keep_data = args[["keep_data"]], 
          cores = args[["cores"]])
  })
  
  ## estimation of the other trees
  if (d > 2 && (is.na(args[["trunc_lvl"]]) || args[["trunc_lvl"]] > 1)) {
    udata_first_tree <- sapply(2:d, function(k) 
      hbicop(udata[, c(1, k)], cond_var = 1, res$first_tree[[k - 1]]))
    res$vine <- vinecop(udata_first_tree, 
                        family_set = args[["family_set"]], 
                        par_method = args[["par_method"]], 
                        nonpar_method = args[["nonpar_method"]], 
                        mult = args[["mult"]][k - 1], 
                        selcrit = args[["selcrit"]], 
                        presel = args[["presel"]], 
                        trunc_lvl = ifelse(is.na(args[["trunc_lvl"]]), 
                                           NA, args[["trunc_lvl"]] - 1), 
                        tree_crit = args[["tree_crit"]], 
                        threshold = args[["threshold"]], 
                        keep_data = args[["keep_data"]], 
                        show_trace = args[["show_trace"]], 
                        cores = args[["cores"]])
  } else {
    res$vine = NULL
  }
  
  class(res) <- "eevine"
  res
}

#' Working with a `eevine` object
#'
#' @aliases dkeeevine
#'
#' @param x \eqn{m x 2} matrix of evaluation points (should be in (0,1) if 
#' \code{copula_data = TRUE}).
#' @param object `eevine` object.
#' @param copula_data if `TRUE``, returns the copula density only.
#'
#' @return A numeric vector for the density.
#'
#' @seealso \code{\link{eevine}}, \code{\link[kdevine]{dkdevine}}
#'
#' @examples
#' ## TODO
#'
#' @importFrom kdevine dkde1d pkde1d
#' @importFrom rvinecopulib dbicop hbicop dvinecop
#' @export
deevine <- function(x, object, copula_data = FALSE) {
  
  assert_that(class(object) == "eevine")
  assert_that(ncol(x) == length(object$margin_models))
  d <- ncol(x)
  
  ## contributions of the marginals
  if (copula_data == FALSE) {
    marg_vals <- sapply(1:d, function(k) 
      dkde1d(x[, k], object$margin_models[[k]]))
    udata <- sapply(1:d, function(k) pkde1d(x[, k], object$margin_models[[k]]))
  } else {
    marg_vals <- rep(1, nrow(x))
    udata <- x
  }
  
  ## contributions of the first tree
  first_tree_vals <- sapply(2:d, function(k) 
    dbicop(udata[, c(1, k)], object$first_tree[[k - 1]]))
  
  ## contribution of the other trees
  if (!is.null(object$vine)) {
    udata_first_tree <- sapply(2:d, function(k) 
      hbicop(udata[, c(1, k)], cond_var = 1, object$first_tree[[k - 1]]))
    vine_vals <- dvinecop(udata_first_tree, object$vine)
  } else {
    vine_vals <- rep(1, nrow(x))
  }
  apply(cbind(marg_vals, first_tree_vals, vine_vals), 1, prod)
}

## the copula weights
getcxy <- function(x, uy, object) {
  d <- length(object$margin_models)
  n <- length(uy)
  ux <- sapply(2:d, function(k) pkde1d(x[k - 1], object$margin_models[[k]]))
  ux <- matrix(data = rep(ux, each = n), nrow = n, ncol = d - 1)
  deevine(cbind(uy, ux), object, copula_data = TRUE)
}

#' Predictions for kernel based C-vine model
#' 
#' Use a fitted C-vine model to predict the mean or quantiles of the response 
#' variable (i.e., the central node).
#'
#' @param object a `eevine` object.
#' @param newdata points where the fit shall be evaluated.
#' @param what what to predict, one or two of `'mean'`, `'sd'`, `'var'`, 
#' or `'quantile'`.
#' @param ... if `what='quantile'`, then `tau` is used to specify the quantile 
#' of interest (default `tau = c(0.1,0.25,0.5,0.75,0.9)`).
#'
#' @return Predictions, either in a vector (if `what='mean'`, `what='sd'` or 
#' `what='var'`), a matrix (if `what='quantile'`), or list (if e.g. 
#' `what=c('mean', 'quantile')`).
#'
#' @examples
#' set.seed(0)
#' 
#' # Quantiles of interest
#' alpha <- c(0.1,0.5,0.9)
#' 
#' ## The predictors
#' p <- 10
#' n <- 1e3
#' X <- matrix(2 * runif(n * p) - 1, n, p)
#' colnames(X) <- paste0('X', 1:p)
#' 
#' ## mean and variance shifts
#' Y1 <- rnorm(n) + exp(X[,1])
#' Y2 <- rnorm(n) * (1 + exp(X[,1]))
#' 
#' ## For the evaluation of the results
#' ticks = 1e2
#' X.test = matrix(0, ticks, p)
#' X.test[,1] = seq(-0.9, 0.9, length.out = ticks)
#' colnames(X.test) <- paste0('X', 1:p)
#' 
#' ## The truth
#' truth1 <- cbind(qnorm(0.1) + exp(X.test[,1]),
#'                 exp(X.test[,1]),
#'                 qnorm(0.9) + exp(X.test[,1]))
#' truth2 <- cbind(qnorm(0.1) * (1 + exp(X.test[,1])),
#'                 0,
#'                 qnorm(0.9) * (1 + exp(X.test[,1])))
#' 
#' ## range for the margin estimates
#' xmin <- c(-Inf,rep(-1, p))
#' xmax <- c(Inf,rep(1, p))
#' 
#' ## Can we identify the right structure?
#' rvr1 <- eevine(x = X, y = Y1, xmin = xmin, xmax = xmax,
#'                  test.level = 0.05, trunc.level = 1, treecrit = 'hoeffd')
#' rvr2 <- eevine(x = X, y = Y2, xmin = xmin, xmax = xmax,
#'                  test.level = 0.05, trunc.level = 1, treecrit = 'hoeffd')
#' 
#' ## Quantile regression results
#' rq1 <- predict(rvr1, X.test, 'quantile', tau = c(0.1, 0.5, 0.9))
#' rq2 <- predict(rvr2, X.test, 'quantile', tau = c(0.1, 0.5, 0.9))
#' 
#' ## Plot the results
#' par(mfrow = c(1,2))
#' plot(X[,1], Y1, xlab = 'X1', ylab = 'Y1')
#' matlines(X.test[,1], truth1, col = 'green', lty = c(2,1,2), lwd = 2)
#' matlines(X.test[,1], rq1, col = 'red', lty = c(2,1,2), lwd = 2)
#' legend('topleft', legend = c('data', 'truth', 'estimate'),
#'        col = c(1, 'green', 'red'),
#'        pch = c(1, rep(NA,2)),
#'        lty = c(0, rep(1,2)),
#'        cex = 0.7)
#' plot(X[,1], Y2, xlab = 'X1', ylab = 'Y2')
#' matlines(X.test[,1], truth2, col = 'green', lty = c(2,1,2), lwd = 2)
#' matlines(X.test[,1], rq2, col = 'red', lty = c(2,1,2), lwd = 2)
#'
#' @importFrom kdevine dkde1d pkde1d
#' @importFrom rvinecopulib dbicop hbicop dvinecop
#' @export
predict.eevine <- function(object, newdata = NULL, what = "mean", ...) {
  
  assert_that(all(is.element(what, c("mean", "sd", "var", "quantile", "cxy"))))
  what <- unique(what)
  
  if (any(what == "quantile")) {
    if (!is.null(list(...)[["tau"]])) {
      tau <- list(...)[["tau"]]
      assert_that(is.vector(tau) && is.numeric(tau) && all(tau > 0 & tau < 1))
    } else {
      tau <- c(0.1, 0.25, 0.5, 0.75, 0.9)
    }
  }
  
  d <- length(object$margin_models)
  if (!is.null(newdata)) {
    assert_that((is.vector(newdata) && d == 2) || 
                  (is.matrix(newdata) && ncol(newdata) == d - 1))
    if (is.vector(newdata)) {
      newdata <- matrix(newdata, ncol = 1)
    }
  } else {
    newdata <- sapply(2:d, function(k) object$margin_models[[k]]$x_cc)
  }
  
  y <- as.numeric(object$margin_models[[1]]$x_cc)
  uy <- pkde1d(y, object$margin_models[[1]])
  cxy <- apply(newdata, 1, function(x) getcxy(x, uy, object))
  
  f <- function(type) {
    switch(type, 
           mean = conditional_mean(y, cxy), 
           sd = sqrt(conditional_variance(y, cxy)), 
           var = conditional_variance(y, cxy), 
           quantile = conditional_quantile(y, cxy, tau),
           cxy = cxy)
  }
  
  if (length(what) > 1) {
    res <- lapply(what, f)
  } else {
    res <- f(what)
  }
  
  return(res)
}

conditional_mean <- function(y, cxy) {
  as.numeric(y %*% cxy)/apply(cxy, 2, sum)
}

conditional_variance <- function(y, cxy) {
  z <- apply(cxy, 2, sum)
  m1 <- as.numeric(y %*% cxy)/z
  m2 <- as.numeric(y^2 %*% cxy)/z
  m2 - m1^2
}

conditional_quantile <- function(y, cxy, tau) {
  t(apply(cxy, 2, function(w) wquantile(y, w, tau)))
}

## weighted quantiles (inspired from bigvis)
wquantile <- function(x, w, probs = seq(0, 1, 0.25)) {
  
  # Ensure x and w in ascending order of x
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  
  # Find closest x just below and above index
  n <- sum(w)
  index <- 1 + (n - 1) * probs
  j <- floor(index)
  
  wts <- cumsum(w)
  lo <- x[sapply(j, function(i) which(wts - i > 0)[1])]  # X_j
  hi <- x[sapply(j + 1, function(i) which(wts - i > 0)[1])]
  
  g <- index - j
  ifelse(lo == hi, lo, (1 - g) * lo + g * hi)
}