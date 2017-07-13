sigmoid <- function(x,theta) {
  # if infinity check it
  return (1 / (1 + exp(-1 * x %*% theta)))
}

LLReg <- function(x,y,theta,l, reg.intercept) {
  return (
    1/nrow(y) * ( (-1 * t(y) %*% log(sigmoid(x,theta))) - (t((1 - y)) %*% log(1 - sigmoid(x,theta)))
                  + l/2 * (t(theta[(2 - reg.intercept):nrow(theta), 1]) %*% theta[(2 - reg.intercept):nrow(theta), 1]) )
  )
}

GradientReg <- function(x,y,theta,l, reg.intercept) {
  return (
    1/nrow(y) * ((t(x) %*% (sigmoid(x,theta) - y)) 
                 + l * rbind((1 - reg.intercept), matrix(1, nrow = (ncol(x) - 1), ncol = 1)))
  )
}

HessianReg <- function(x,y,theta,l, reg.intercept) {
  return (
    1/nrow(y) * ( (t(x) %*% diag(as.vector(sigmoid(X,theta) * (1 - sigmoid(X, theta)))) %*% x) 
                  + l * diag(as.vector(rbind((1 - reg.intercept), matrix(1, nrow = (ncol(x) - 1), ncol = 1)))) )
  )
}

NMReg <-function(x,y,l,tolerance, reg.intercept) {
  theta <- matrix(0, ncol = 1, nrow = ncol(x))
  LL.old <- 1000000000
  LL.new <- LLReg(x,y,theta,l, reg.intercept)
  while(LL.old - LL.new > tolerance) {
    LL.old <- LL.new
    theta <- theta - solve(HessianReg(x,y,theta,l, reg.intercept), GradientReg(x,y,theta,l,reg.intercept))
    LL.new <- LLReg(x,y,theta,l, reg.intercept)
    print(LL.new)
  }
  return (theta)
}
