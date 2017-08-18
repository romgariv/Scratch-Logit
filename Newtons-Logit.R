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
    1/nrow(y) * ( (t(x) %*% diag(as.vector(sigmoid(x,theta) * (1 - sigmoid(x, theta)))) %*% x) 
                  + l * diag(as.vector(rbind((1 - reg.intercept), matrix(1, nrow = (ncol(x) - 1), ncol = 1)))) )
  )
}

NMReg <-function(x,y,l,tolerance, reg.intercept, print.LL = FALSE) {
  theta <- matrix(0, ncol = 1, nrow = ncol(x))
  LL.old <- 1000000000
  LL.new <- LLReg(x,y,theta,l, reg.intercept)
  while(LL.old - LL.new > tolerance) {
    LL.old <- LL.new
    theta <- theta - solve(HessianReg(x,y,theta,l, reg.intercept), GradientReg(x,y,theta,l,reg.intercept))
    LL.new <- LLReg(x,y,theta,l, reg.intercept)
    if (print.LL) print(LL.new)
  }
  return (theta)
}

NMLogit <- function(x,y) {
  return(NMReg(x, y, 0, 1e-10, FALSE))
}

CompareLogits <- function(x,y) {
  R.Done <- as.matrix(glm.fit(x, y, family = binomial(link="logit"))$coefficients)
  NM <- NMLogit(x, y)
  comp <- as.data.frame(cbind(R.Done, NM))
  names(comp) <- c("GLM Logit", "Newton's Method")
  comp$Error <- (comp$`Newton's Method` - comp$`GLM Logit`) / comp$`GLM Logit`
  return(comp)
  
}

Logit.LL <- function(theta, x,y, as.scalar = TRUE) {
  return(LLReg(x,y,theta,0,0, as.scalar))
}

Standard.Errors <- function(theta, x) {
  p <- sigmoid(x, theta)
  v <- diag(as.vector(p * (1-p)))
  return(sqrt(diag(solve(t(x) %*% v %*% x))))
}

Numeric.Standard.Errors <- function(theta, x, y, epp = 1e-6) {
  d <- abs(theta * epp)
  h <- matrix(NA, nrow = nrow(theta), ncol = 1)
  for (i in 1:nrow(theta)) {
    h[i] <- max(epp, d[i])
  }
  
  Gi <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  
  for (i in 1:nrow(theta)) {
    t1 <- theta
    t2 <- theta
    t1[i,] <- theta[i,] + h[i]
    t2[i,] <- theta[i,] - h[i]
    Gi[,i] <- (LLReg(x,y,t1,0,0,FALSE) - LLReg(x,y,t2,0,0,FALSE))/(2*h[i])
  }
  return( sqrt(diag(solve(t(Gi) %*% Gi))) )
}

Z.Score <- function(theta, SE) {
  return(theta/SE)
}

P.Value <- function(z) {
  return (1.96 * pnorm(-abs(z)))
}

Logit.Summary <- function(x,y) {
  theta <- NMLogit(x, y)
  se <- Standard.Errors(theta, x)
  z <- Z.Score(theta, se)
  p <- P.Value(z)
  df <- data.frame(Estimates = round(theta,8),
                   `Std. Error` = as.vector(round(se,7)),
                   `Z-Score` = as.vector(round(z,3)),
                   `P-Value` = as.vector(round(p,4)))
  colnames(df) <- c("Estimates", "Std. Error", "Z-Score", "P-Value")
  return(df)
}