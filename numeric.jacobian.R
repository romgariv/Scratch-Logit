NumericJacobian <- function(func, x0,...) {
  epsilon <- 1e-8
  x0.rows <- nrow(x0)
  Jac <- matrix(NA, nrow = x0.rows, ncol = 1)
  
   func0 <- func(x0, ...)

 for (i in 1:x0.rows) {
   dx <- c(rep(0,i-1), epsilon, rep(0,x0.rows - i))
   Jac[i,1] <- (func((x0 + dx), ...) - func0)/epsilon
 }
 return (Jac)
}
