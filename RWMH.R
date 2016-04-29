#Generate 25 observations from a univariate normal distribution with mean 2 and std 2
#Given the simulated data, 
#perform a random walk Metropolis Hastings algorithms for the paramters of a normal distribution, mu & sigma2
#assume joint prior p(mu, sigma2) is proportional to sigma2^(-1)
#tune the algorithm to have an acceptance rate between 0.25 and 0.7

library(coda)
NormalMH <- function(n.sim, n.burnin){
  N <- 25
  #generate 25 obs from a univariate normal (2, 4)
  x <- rnorm(N, mean = 2, sd = 2)
  #defines a function that returns unnormalized posterior density
  post_den <-function(mu, sigma2, x, N){
    #prior
    prior <- sigma2^(-1)
    #normal distr (covariance matrix is sigma2*I)
    likelihood <- ((2*pi)^(-N/2))*sigma2^{-N/2}*exp(-sum((x-mu)^2)/(sigma2*2))
    # unnormalized target distr
    working_cond_den <- prior*likelihood
  }
  #empty matrix to store sampled parameters
  theta.mh <- matrix(NA, nrow = n.sim, ncol = 2) 
  
  rho <-0.2
  #generate inital parameters
  theta.current <- mvrnorm(n=1, mu = c(2,4), Sigma = matrix(c(1,rho,rho,1), nrow=2, ncol=2))
  theta.update <- function(index, theta.current,x,N) {
    
    if (index == 1){ #if it's sample for mu, update the first element of the parameter vector
      #sample a new parameter value from Gaussian
      theta.star <- rnorm(n = 1, mean = theta.current[index], sd =2)
      theta.temp <- c(theta.star, theta.current[2])
    }
    else {
      theta.star <- abs(rnorm(n = 1, mean = theta.current[index], sd =2))
      theta.temp <- c(theta.current[1], theta.star)#for beta, update the second element
      }
      #compute MH ratio
      r <- post_den(theta.temp[1],theta.temp[2],x,N)/
                                      post_den(theta.current[1],theta.current[2],x,N)
      r <- min(r, 1, na.rm = TRUE) # r is the accpetance probability, NA values are ignored
      if (runif(1) < r)
      #runif(1) generates a uniform random number bewteen 0 and 1
      #if runif(1) < r, accept the new value; else reject the new value and keep the old value
        theta.star
      else theta.current[index]
    
  }
    for (i in 1:n.sim) {
    #iteration for mu
    theta.current[1] <- theta.mh[i, 1] <- theta.update(1, theta.current,x, N) 
    #iteration for sigma                                                   
    theta.current[2] <- theta.mh[i, 2] <- theta.update(2, theta.current,x, N)
  }
  theta.mh <- theta.mh[(n.burnin + 1):n.sim, ] #discard burn-in
} 

mh.draws <- NormalMH (n.sim = 10000, n.burnin = 1000)

#turn the result into  an MCMC object
mh.draws <- mcmc(mh.draws)

#output the acceptance rates for both parameters
#fine tune by adjusting the variance of proposal densities
cat("Acceptance rate:", 1- rejectionRate(mh.draws), "\n")

#traceplot
plot(mh.draws)

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


#Gelman and Rubin diagnostics
mh.draws1 <- mcmc(NormalMH(n.sim = 10000, n.burnin = 1000))
mh.draws2 <- mcmc(NormalMH(n.sim = 10000, n.burnin = 1000))
mh.draws3 <- mcmc(NormalMH(n.sim = 10000, n.burnin = 1000))
mh.draws4 <- mcmc(NormalMH(n.sim = 10000, n.burnin = 1000))
mh.draws5 <- mcmc(NormalMH(n.sim = 10000, n.burnin = 1000))
mh.list <- mcmc.list(list(mh.draws1, mh.draws2, mh.draws3,mh.draws4, mh.draws5))
gelman.diag(mh.list)  
gelman.plot(mh.list)
