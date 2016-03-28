P <- function(lambda, alpha, beta, x) {
  p <- (lambda^(alpha-1 + sum(x)))*exp((-beta-dim(x)[1])*lambda)
  return(-p)
}  
result2<-optim(par=6, P, alpha=2, beta=3, x=X, method="Brent", lower=0.01, upper=20)
lambda_mode <- result2$par
print(lambda_mode)
