#' Kernel density estimatior based on simplified C-vine copulas
#'
#' Implements the C-vine based estimator of Nagler and Czado (2016). We refer
#' to \code{\link[kdevine]{kdevine}} for more details.
#'
#' @param x (\eqn{n x (d-1)}) data matrix.
#' @param y (\eqn{n x 1}) data vector for the central node.
#' @param mult_1d numeric; all bandwidhts for marginal kernel density estimation
#'   are multiplied with \code{mult_1d}. Defaults to `log(1 + d)`.
#' @param xmin numeric vector of length d; see \code{\link[kdevine]{kde1d}}.
#' @param xmax numeric vector of length d; see \code{\link[kdevine]{kde1d}}.
#' @param ... further arguments passed to \code{\link[kdevine]{kde1d}} or
#'   \code{\link[kdevine]{kdevinecop}}.
#'
#' @return An object of class \code{kdecvine}.
#' 
#' @author Thibault Vatter
#' 
#' @references Nagler, T., Czado, C. (2016) *Evading the curse of
#'   dimensionality in nonparametric density estimation with simplified vine
#'   copulas.* Journal of Multivariate Analysis 151, 69-89
#'   (doi:10.1016/j.jmva.2016.07.003)
#'
#' @seealso \code{\link{dkdecvine}} \code{\link[kdevine]{kdevine}} \code{\link[kdevine]{kdevine}}
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
#' colnames(X) <- paste0("X", 1:p)
#' 
#' ## mean and variance shifts
#' Y1 <- rnorm(n) + exp(X[,1])
#' Y2 <- rnorm(n) * (1 + exp(X[,1]))
#' 
#' ## For the evaluation of the results
#' ticks = 1e2
#' X.test = matrix(0, ticks, p)
#' X.test[,1] = seq(-0.9, 0.9, length.out = ticks)
#' colnames(X.test) <- paste0("X", 1:p)
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
#' rvr1 <- kdecvine(x = X, y = Y1, xmin = xmin, xmax = xmax,
#'                  test.level = 0.05, trunc.level = 1, treecrit = "hoeffd")
#' rvr2 <- kdecvine(x = X, y = Y2, xmin = xmin, xmax = xmax,
#'                  test.level = 0.05, trunc.level = 1, treecrit = "hoeffd")
#' 
#' ## Quantile regression results
#' rq1 <- predict(rvr1, X.test, "quantile", tau = c(0.1, 0.5, 0.9))
#' rq2 <- predict(rvr2, X.test, "quantile", tau = c(0.1, 0.5, 0.9))
#' 
#' ## Plot the results
#' par(mfrow = c(1,2))
#' plot(X[,1], Y1, xlab = "X1", ylab = "Y1")
#' matlines(X.test[,1], truth1, col = "green", lty = c(2,1,2), lwd = 2)
#' matlines(X.test[,1], rq1, col = "red", lty = c(2,1,2), lwd = 2)
#' legend("topleft", legend = c("data", "truth", "estimate"),
#'        col = c(1, "green", "red"),
#'        pch = c(1, rep(NA,2)),
#'        lty = c(0, rep(1,2)),
#'        cex = 0.7)
#' plot(X[,1], Y2, xlab = "X1", ylab = "Y2")
#' matlines(X.test[,1], truth2, col = "green", lty = c(2,1,2), lwd = 2)
#' matlines(X.test[,1], rq2, col = "red", lty = c(2,1,2), lwd = 2)
#'
#' @importFrom kdevine kde1d pkde1d kdevinecop
#' @importFrom kdecopula kdecop hkdecop
#' @export
kdecvine <- function(x, y, mult_1d = 1, xmin = NULL, xmax = NULL, ...) {
  
  if (is.vector(x)) {
    d <- 2
  } else {
    d <- ncol(x)+1
  }
  data <- cbind(y,x)
  n <- nrow(data)

  args <- list(...)
  if (is.null(mult_1d)) 
    mult_1d <- log(d+1)
  if (is.null(args$info)) 
    args$info <- FALSE
  if (is.null(args$method)) 
    args$method <- "TLL2"
  if (is.null(args$mult)) 
    args$mult <- 1
  if (is.null(args$test.level)) 
    args$test.level <- 1
  if (is.na(args$test.level)) 
    args$test.level <- 1
  if (is.null(args$trunc.level)) 
    args$trunc.level <- d
  if (is.na(args$trunc.level)) 
    args$trunc.level <- d
  if (is.null(args$treecrit)) 
    args$treecrit <- "hoeffd"
  if (is.na(args$treecrit)) 
    args$treecrit <- "hoeffd"
  if (is.null(args$cores)) 
    args$cores <- 1
  
  ## estimation of the margins
  marg.dens <- as.list(numeric(d))
  for (k in 1:d) {
    marg.dens[[k]] <- kde1d(data[, k],
                            xmin = args$xmin[k],
                            xmax = args$xmax[k],
                            bw   = args$bw[k],
                            mult = mult_1d)
  }
  res <- list(marg.dens = marg.dens)
  
  ## estimation of the first tree
  udata <- sapply(1:d, function(k) pkde1d(data[, k], marg.dens[[k]]))
  res$first_tree <- fit_first_tree(udata, args$method, args$mult, 
                                   args$info, args$test.level)
  
  ## estimation of the other trees
  if (d > 2 && args$trunc.level > 1) {
    udata_first_tree <- sapply(2:d, function(k) 
      hkdecop(udata[,c(1,k)], res$first_tree[[k-1]], cond.var = 1))
    res$vine <- kdevinecop(udata_first_tree, 
                           method = args$method, 
                           mult = args$mult, 
                           info = args$info, 
                           test.level = args$test.level, 
                           trunc.level = args$trunc.level - 1, 
                           treecrit = args$treecrit, 
                           cores = args$cores)
  } else {
    res$vine = NULL
  }
  
  class(res) <- "kdecvine"
  res
}

fit_first_tree <- function(udata, method, mult, info, test.level) {
  d <- ncol(udata)
  indepinfo <- list(effp = 0,
                    likvalues = rep(1, nrow(udata)),
                    loglik = 0,
                    effp = 0,
                    AIC = 0,
                    cAIC = 0,
                    BIC = 0)
  lapply(2:d, function(k) {
    indep <- ifelse(test.level < 1,
                    hoeffd(udata[,c(1,k)])$p.value >= test.level,
                    FALSE)
    if (indep) {
      pcfit <- list()
      if (info)
        pcfit$info <- indepinfo
      class(pcfit) <- c("kdecopula", "indep.copula")
    } else {
      pcfit <- kdecop(udata[,c(1,k)], method = method,
                      mult = mult, info = info)
    }
    return(pcfit)
  })
}

#' Working with a \code{kdecvine} object
#'
#' A C-vine density estimate (stored in a \code{kdecvine} object)
#' can be evaluated on arbitrary points with \code{dkecvine}.
#'
#' @aliases dkecvine
#'
#' @param x \eqn{m x 2} matrix of evaluation points (should be in (0,1) if 
#' \code{copula_data = TRUE}).
#' @param object \code{kdecvine} object.
#' @param copula_data if \code{TRUE}, returns the copula density only.
#'
#' @return A numeric vector of the density.
#'
#' @author Thibault Vatter
#'
#' @seealso \code{\link{kdecvine}}, \code{\link[kdevine]{dkdevine}}
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
#' colnames(X) <- paste0("X", 1:p)
#' 
#' ## mean and variance shifts
#' Y1 <- rnorm(n) + exp(X[,1])
#' Y2 <- rnorm(n) * (1 + exp(X[,1]))
#' 
#' ## For the evaluation of the results
#' ticks = 1e2
#' X.test = matrix(0, ticks, p)
#' X.test[,1] = seq(-0.9, 0.9, length.out = ticks)
#' colnames(X.test) <- paste0("X", 1:p)
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
#' rvr1 <- kdecvine(x = X, y = Y1, xmin = xmin, xmax = xmax,
#'                  test.level = 0.05, trunc.level = 1, treecrit = "hoeffd")
#' rvr2 <- kdecvine(x = X, y = Y2, xmin = xmin, xmax = xmax,
#'                  test.level = 0.05, trunc.level = 1, treecrit = "hoeffd")
#' 
#' ## Quantile regression results
#' rq1 <- predict(rvr1, X.test, "quantile", tau = c(0.1, 0.5, 0.9))
#' rq2 <- predict(rvr2, X.test, "quantile", tau = c(0.1, 0.5, 0.9))
#' 
#' ## Plot the results
#' par(mfrow = c(1,2))
#' plot(X[,1], Y1, xlab = "X1", ylab = "Y1")
#' matlines(X.test[,1], truth1, col = "green", lty = c(2,1,2), lwd = 2)
#' matlines(X.test[,1], rq1, col = "red", lty = c(2,1,2), lwd = 2)
#' legend("topleft", legend = c("data", "truth", "estimate"),
#'        col = c(1, "green", "red"),
#'        pch = c(1, rep(NA,2)),
#'        lty = c(0, rep(1,2)),
#'        cex = 0.7)
#' plot(X[,1], Y2, xlab = "X1", ylab = "Y2")
#' matlines(X.test[,1], truth2, col = "green", lty = c(2,1,2), lwd = 2)
#' matlines(X.test[,1], rq2, col = "red", lty = c(2,1,2), lwd = 2)
#'
#' @importFrom kdevine dkde1d pkde1d dkdevinecop
#' @importFrom kdecopula dkdecop hkdecop
#' @export
dkdecvine <- function(x, object, copula_data = FALSE) {
  
  stopifnot(class(object) == "kdecvine")
  stopifnot(ncol(x) == length(object$marg.dens))
  d <- ncol(x)
  
  ## contributions of the marginals
  if (copula_data == FALSE) {
    marg_vals <- sapply(1:d, function(k) dkde1d(x[, k], object$marg.dens[[k]]))
    udata <- sapply(1:d, function(k) pkde1d(x[, k], object$marg.dens[[k]]))
  } else {
    marg_vals <- rep(1, nrow(x))
    udata <- x
  }
  
  ## contributions of the first tree
  first_tree_vals <- sapply(2:d, function(k) dkdecop(udata[,c(1,k)], 
                                                     object$first_tree[[k-1]]))
  
  ## contribution of the other trees
  if (!is.null(object$vine)) {
    udata_first_tree <- sapply(2:d, function(k) 
      hkdecop(udata[,c(1,k)], object$first_tree[[k-1]], cond.var = 1))
    vine_vals <- dkdevinecop(udata_first_tree, object$vine, stable = TRUE)
  }
  else {
    vine_vals <- rep(1, nrow(x))
  }
  apply(cbind(marg_vals, first_tree_vals, vine_vals), 1, prod)
}

## the copula weights
getcxy <- function(x, uy, object) {
  d <- length(object$marg.dens)
  n <- length(uy)
  ux <- sapply(2:d, function(k) pkde1d(x[k-1], object$marg.dens[[k]]))
  ux <- matrix(data = rep(ux, each = n), nrow = n, ncol = d-1)
  dkdecvine(cbind(uy, ux), object, copula_data = TRUE)
}

#' Predictions for kernel based C-vine model
#' 
#' Use a fitted C-vine model to predict the mean or quantiles of the response 
#' variable (i.e., the central node).
#'
#' @param object a `kdecvine` object.
#' @param newdata points where the fit shall be evaluated.
#' @param what what to predict, one or two of `"mean"`, `"sd"`, `"var"`, 
#' or `"quantile"`.
#' @param ... if `what="quantile"`, then `tau` is used to specify the quantile 
#' of interest (default `tau = c(0.1,0.25,0.5,0.75,0.9)`).
#'
#' @return Either a vector (if `what="mean"`), a matrix (if `what="quantile"`), 
#' or list (if `what=c("mean", "quantile")`) of predictions.
#'
#' @author Thibault Vatter
#'
#' @seealso \code{\link{kdecvine}}, \code{\link{dkdecvine}}
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
#' colnames(X) <- paste0("X", 1:p)
#' 
#' ## mean and variance shifts
#' Y1 <- rnorm(n) + exp(X[,1])
#' Y2 <- rnorm(n) * (1 + exp(X[,1]))
#' 
#' ## For the evaluation of the results
#' ticks = 1e2
#' X.test = matrix(0, ticks, p)
#' X.test[,1] = seq(-0.9, 0.9, length.out = ticks)
#' colnames(X.test) <- paste0("X", 1:p)
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
#' rvr1 <- kdecvine(x = X, y = Y1, xmin = xmin, xmax = xmax,
#'                  test.level = 0.05, trunc.level = 1, treecrit = "hoeffd")
#' rvr2 <- kdecvine(x = X, y = Y2, xmin = xmin, xmax = xmax,
#'                  test.level = 0.05, trunc.level = 1, treecrit = "hoeffd")
#' 
#' ## Quantile regression results
#' rq1 <- predict(rvr1, X.test, "quantile", tau = c(0.1, 0.5, 0.9))
#' rq2 <- predict(rvr2, X.test, "quantile", tau = c(0.1, 0.5, 0.9))
#' 
#' ## Plot the results
#' par(mfrow = c(1,2))
#' plot(X[,1], Y1, xlab = "X1", ylab = "Y1")
#' matlines(X.test[,1], truth1, col = "green", lty = c(2,1,2), lwd = 2)
#' matlines(X.test[,1], rq1, col = "red", lty = c(2,1,2), lwd = 2)
#' legend("topleft", legend = c("data", "truth", "estimate"),
#'        col = c(1, "green", "red"),
#'        pch = c(1, rep(NA,2)),
#'        lty = c(0, rep(1,2)),
#'        cex = 0.7)
#' plot(X[,1], Y2, xlab = "X1", ylab = "Y2")
#' matlines(X.test[,1], truth2, col = "green", lty = c(2,1,2), lwd = 2)
#' matlines(X.test[,1], rq2, col = "red", lty = c(2,1,2), lwd = 2)
#'
#' @importFrom kdevine dkde1d pkde1d dkdevinecop
#' @importFrom kdecopula dkdecop hkdecop
#' @export
predict.kdecvine <- function(object, newdata = NULL, what = "mean", ...) {

  stopifnot(all(is.element(what, c("mean", "sd", "var", "quantile"))))
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
      "mean" = getMean(y, cxy),
      "sd" = sqrt(getVar(y, cxy)),
      "var" = getVar(y, cxy),
      "quantile" = getQuantile(y, cxy, tau))
  }
  
  if (length(what) > 1) {
    res <- lapply(what, f)
  } else {
    res <- f(what)
  }
  
  return(res)
}

getMean <- function(y, cxy) {
  as.numeric(y %*% cxy)/apply(cxy,2,sum)
  # g <- function(par, x, w) {
  #   return(matrix(w*(par - x), ncol = 1))
  # }
  # my <- mean(y)
  # yy <- range(y)
  # apply(cxy, 2, function(w) gmm(function(par, x) g(par, x, w), y, my, 
  #                               optfct="optimize", interval=yy)$coefficients)
}

getVar <- function(y, cxy) {
  z <- apply(cxy,2,sum)
  m1 <- as.numeric(y %*% cxy)/z
  m2 <- as.numeric(y^2 %*% cxy)/z
  m2-m1^2
}

getQuantile <- function(y, cxy, tau) {
  t(apply(cxy, 2, function(w) wquantile(y, w, tau)))
  # res <- t(apply(cxy, 2, function(w) wquantile(y, w, tau)))
  # g <- function(theta, x, w, tau) {
  #   ind <- t(sapply(theta, function(t) ifelse(x > t, 1, 0)))
  #   return(w*t(ind * tau - (1-tau)*(1-ind)))
  # }
  # yy <- range(y)
  # res2 <- t(apply(cxy, 2, function(w) sapply(tau, function(t)
  #   uniroot(function(theta) sum(g(theta, y, w, t)), interval = yy)$root)))
  # my <- mean(y)
  # res3 <- t(apply(cxy, 2, function(w) sapply(tau, function(t)
  #   gmm(function(theta, x) g(theta, x, w, t), y, my, 
  #       optfct="optimize", interval=yy)$coefficients)))
  # my <- quantile(y, tau)
  # res4 <- t(apply(cxy, 2, function(w)  
  #   gmm(function(theta, x) g(theta, x, w, tau), y, my)$coefficients))
  # matplot(res, type = "l", lty = 1, col = 1)
  # matlines(res2, col = 2)
  # matlines(res3, col = 3)
  # matlines(res4, col = 4)
}


## inspired from bigvis
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
  lo <- x[sapply(j, function(i) which(wts - i > 0)[1])]     # X_j
  hi <- x[sapply(j+1, function(i) which(wts - i > 0)[1])]
  
  g <- index - j
  ifelse(lo == hi, lo, (1 - g) * lo + g * hi)
}

hoeffd <- function(x) {
  n <- nrow(x)
  R <- apply(x,2,rank)-1
  Q <- sapply(1:n, function(i) sum(x[,1] < x[i,1] & x[,2] < x[i,2]))
  A <- sum(apply(R*(R-1),1,prod))
  B <- sum(apply((R-1),1,prod)*Q)
  C <- sum(Q*(Q-1))
  D <- (A - 2*(n-2)*B + (n-2)*(n-3)*C)/(n*(n-1)*(n-2)*(n-3)*(n-4))
  return(list(D = D, p.value = phoeffd(D, n)))
}

## inspired from Hmisc
#' @importFrom stats approx
phoeffd <- function(d, n)
{
  b <- d + 1 / 36 / n
  z <- .5 * (pi ^ 4) * n * b
  zz <- as.vector(z)
  zz[is.na(zz)] <- 1e50   # so approx won't bark
  
  tabvals <- c(5297,4918,4565,4236,3930,
               3648,3387,3146,2924,2719,2530,2355,
               2194,2045,1908,1781,1663,1554,1453,
               1359,1273,1192,1117,1047,0982,0921,
               0864,0812,0762,0716,0673,0633,0595,
               0560,0527,0496,0467,0440,0414,0390,
               0368,0347,0327,0308,0291,0274,0259,
               0244,0230,0217,0205,0194,0183,0173,
               0163,0154,0145,0137,0130,0123,0116,
               0110,0104,0098,0093,0087,0083,0078,
               0074,0070,0066,0063,0059,0056,0053,
               0050,0047,0045,0042,0025,0014,0008,
               0005,0003,0002,0001)/10000
  
  ifelse(z < 1.1 | z > 8.5,
         pmax(1e-8, pmin(1, exp(.3885037 -1.164879 * z))),
         approx(c(seq(1.1, 5, by=.05),
                  seq(5.5,8.5,by=.5)),
                tabvals, zz)$y)
}