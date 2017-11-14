# Implementation from https://johanndejong.wordpress.com/2016/01/17/svm-with-recursive-feature-elimination/

svm_rfe <- function(X, y, elim_frac = 1 / ncol(X), ...) {
  
  library(kernlab)
  
  # keep track of the iteration during which
  # a feature was eliminated
  ii <- rep(NA, ncol(X))
  i <- 0
  while ( any(is.na(ii)) ) {
    # indices of remaining features
    not_elim_yet <- which(is.na(ii))
    # number of features to eliminate
    n_to_elim <- ceiling ( elim_frac * length(not_elim_yet) )
    # train the classifier on the remaining features
    fit <- ksvm(X[,not_elim_yet], y, scaled = FALSE, ...)
    # compute the primal problem coefficients from the dual
    # problem coefficients
    sv_i <- alphaindex(fit)[[1]]
    w <- t( coef(fit)[[1]] ) %*% X[ sv_i, not_elim_yet ]
    # eliminate the features with the smallest squared weights
    to_elim <- not_elim_yet[ head(order( w * w ), n_to_elim) ]
    ii[to_elim] <- i
    i <- i + 1
  }
  # convert iterations into ranks
  i - ii
  
}