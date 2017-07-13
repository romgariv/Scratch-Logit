NumericHessian <- function(func, x0, ...) {
  epsilon <- 1e-8
  x0.rows <- nrow(x0)
  Hess <- matrix(NA, nrow = x0.rows, ncol = x0.rows )
  for (i in 1:x0.rows) {
    x1 <- x0
    x2 <- x0
    x1[i,1] <- x0[i,1] - epsilon
    x2[i,1] <- x0[i,1] + epsilon
    
    NJ1 <- NumericJacobian(func, x1, ...)
    NJ2 <- NumericJacobian(func, x2, ...)
    
    Hess[i, ] <- t((NJ2 - NJ1) / (2 * epsilon))
  }
  return(Hess)
}