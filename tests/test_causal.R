rm(list=ls())
library(gradient.forest)
library(gmmcop)
set.seed(0)

## The predictors
p <- 10
n <- 1e3
X <- matrix(rnorm(n*p), n, p)
colnames(X) <- paste0("X", 1:p)

## Treatment
W <- rbinom(n, 1, 0.5)
Y <- pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)

## For the evaluation of the results
ticks <- 1e2
X.test <- matrix(0, ticks, p)
X.test[,1] <- seq(-2, 2, length.out = ticks)
colnames(X.test) <- paste0("X", 1:p)

## The truth
truth <- pmax(X.test[,1], 0)

fit <- kdemcvine(X, Y, W, test.level = 0.05, trunc.level = 1, treecrit = "hoeffd")
fit2 <- kdecvine(X, Y, test.level = 0.05, trunc.level = 1, treecrit = "hoeffd")
predict(fit2, X.test, what = "mean", w = W)

tau.forest = causal.forest(X, Y, W)
tau.hat = predict(tau.forest, X.test)
plot(X.test[,1], tau.hat$predictions, ylim = range(tau.hat$predictions, 0, 2), xlab = "x", ylab = "tau", type = "l")
lines(X.test[,1], pmax(0, X.test[,1]), col = 2, lty = 2)