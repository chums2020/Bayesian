##########################Random Walk Metropolis-Hastings#########################
  #generate some "observed" data
  #choose the starting values alpha & beta
  #generate a new alpha & beta (assume independence? and sample from univariate Gaussian)
  #compute acceptance probability
  #accept or not
#proposal_density, working_cond_den,k, n, alpha, beta
GammaMHExample <- function(n.sim, n.burnin){
  alpha <- 3; beta <-3 #pick the parameter values to generate a random "observed" data
  N <- 30 #set observed sample size
  x <- rgamma(N, shape = alpha, scale = 1/beta)#generate observed data of size N
  
  #defines a function that returns unnormalized posterior density
  post_den <-function(alpha, beta, x, N){
    prior <- sqrt(alpha*trigamma(alpha)-1)/beta#Jeffery's prior
    likelihood <- prod(x^(alpha-1))*prod(exp(-beta*x))*((beta^alpha)/gamma(alpha))^N   #gamma distr
    working_cond_den <- prior*likelihood # unnormalized target distr
  }
  
  theta.mh <- matrix(NA, nrow = n.sim, ncol = 2) #empty matrix to store sampled parameters
  theta.current <- rnorm(n = 2, mean = 3, sd = 0.5) #generate inital parameters
  theta.update <- function(index, theta.current,x,N) {
    #sample a new parameter value from Gaussian
    theta.star <- rnorm(n = 1, mean = theta.current[index], sd = 0.5)
    if (index == 1) #if it's sample for alpha, update the first element of the parameter vector
      theta.temp <- c(theta.star, theta.current[2])
    else theta.temp <- c(theta.current[1], theta.star)#for beta, update the second element
    #compute MH ratio
    r <- post_den(theta.temp[1],theta.temp[2],x,N)/post_den(theta.current[1],theta.current[2],x,N)
    r <- min(r, 1, na.rm = TRUE) # r is the accpetance probability
    if (runif(1) < r)
      #runif(1) generates a uniform random number bewteen 0 and 1
      #if runif(1) < r, accept the new value; else reject the new value and keep the old value
      theta.star
    else theta.current[index]
  }
    for (i in 1:n.sim) {
    theta.current[1] <- theta.mh[i, 1] <- theta.update(1, theta.current,x, N) #iteration for alpha
                                                       
    theta.current[2] <- theta.mh[i, 2] <- theta.update(2, theta.current,x, N)#iteration for beta
                                                       
  }
  theta.mh <- theta.mh[(n.burnin + 1):n.sim, ] #discard burn-in
} 

mh.draws <- GammaMHExample(n.sim = 10000, n.burnin = 1000)


##########################Diagnostics #########################
library(coda)

mh.draws <- mcmc(mh.draws)#turn the chain into an mcmc object
summary(mh.draws)

#Traceplot
plot(mh.draws)

#Autocorrelation function plots
autocorr.plot(mh.draws)

#Gelman and Rubin

#Geweke 
geweke.diag(mh.draws)


#Raftery and Lewis
raftery.diag(mh.draws,q=0.025,r=0.005,s=0.95)

#Heidelberg and Welch
heidel.diag(mh.draws)


#Plot the MCMC estimate of the marginal posterior density for each parameter with an 95% HPD interval for each parameter
library(TeachingDemos)
d1 <- density(mh.draws[,1])
plot(d1, main="Marginal Posterior Density of alpha")
abline(v=emp.hpd(mh.draws[,1], conf=0.95), col = "red")
emp.hpd(mh.draws[,1], conf=0.95) #95% HPD interval for alpha

d2 <- density(mh.draws[,2])
plot(d2, main="Marginal Posterior Density of beta")
abline(v=emp.hpd(mh.draws[,2], conf=0.95), col = "red")
emp.hpd(mh.draws[,2], conf=0.95) #95% HPD interval for beta
