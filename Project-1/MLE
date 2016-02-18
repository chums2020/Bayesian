X=data.frame(rpois(40,6))
L <- function(lambda,x) {
  l <- 0
  for(i in 1:dim(x)[1]) {
    l <-l +x[i,1]*log(lambda)-lambda-log(factorial(x[i,1]))
  }
  return(-l) 
  #default setting of optim is minimization, 
  #thus we negate the function to make it a maximization problem
}
result <- optim(par = 6, L, x = X,method = "Brent", lower = 0.01, upper = 100)
lambda_hat <- result$par
print(lambda_hat)


%%
X=data.frame(rexp(50,3))

% method 1
L <- function(lambda,x) {
  l <- 0
  for(i in 1:dim(x)[1]) {
    l <-l +log(lambda)-lambda*x[i,1]
  }
  return(-l) 
  #default setting of optim is minimization, 
  #thus we negate the function to make it a maximization problem
}

%method 2
L <- function(lambda,x) {
    l <-dim(x)[1]*log(lambda)-lambda*sum(x)
    return(-l) 
  #default setting of optim is minimization, 
  #thus we negate the function to make it a maximization problem
}  
    
result <- optim(par = 3, L, x = X,method = "Brent", lower = 0.01, upper = 100)
lambda_hat <- result$par
print(lambda_hat)
