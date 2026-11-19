EPS <- 1e-9

# sum of x by integer index into a length-n vector
.rsum <- function(x, idx, n) {
  out <- numeric(n)
  if (!length(x)) return(out)
  s <- rowsum(x, idx)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

.say <- function(verbose, t0, ...) {
  if (isTRUE(verbose)) message(sprintf("[edgeicg %5.0f s] ", proc.time()[[3]] - t0), sprintf(...))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
