#Perform Bayesian linear regression analysis for the dataset stackloss
#Estimate a 95% HPD interval for each parameter
#Parameters: beta's for intercept and 4 variables, and sigma2 (error variance)
#Prior assumptions: beta has distribution N(b, inv(B))
#sigma-2 has distribution Gamma(c/2, d/2), and beta, sigma2 a priori independent
#Default: b=0, B=0, i.e., the default prior of beta is the uniform distribution on an infinite interval
#Default: c= 0.001, d= 0.001

data <-stackloss
library(MCMCpack)
Bfit <- MCMCregress(stack.loss ~ Air.Flow + Water.Temp + Acid.Conc. ,
                    data = data, burning = 1000, mcmc= 25000, thin =  25)
summary(Bfit)
plot(Bfit)

#95% HPD interval for each paramter
library(TeachingDemos)
par(mfrow = c(3,2))
d1 <- density(Bfit[,1])
plot(d1, main="Marginal Posterior Density of Intercept")
emp.hpd(Bfit[,1], conf=0.95) 
abline(v=emp.hpd(Bfit[,1], conf=0.95), col = "red")

d2 <- density(Bfit[,2])
plot(d2, main="Marginal Posterior Density of Air Flow Parameter")
emp.hpd(Bfit[,2], conf=0.95) 
abline(v=emp.hpd(Bfit[,2], conf=0.95), col = "red")

d3 <- density(Bfit[,3])
plot(d3, main="Marginal Posterior Density of Water Temp Parameter")
emp.hpd(Bfit[,3], conf=0.95) 
abline(v=emp.hpd(Bfit[,3], conf=0.95), col = "red")

d4 <- density(Bfit[,4])
plot(d4, main="Marginal Posterior Density of Acid Conc Parameter")
emp.hpd(Bfit[,4], conf=0.95) 
abline(v=emp.hpd(Bfit[,4], conf=0.95), col = "red")

d5 <- density(Bfit[,5])
plot(d5, main="Marginal Posterior Density of sigma2 ")
emp.hpd(Bfit[,5], conf=0.95) 
abline(v=emp.hpd(Bfit[,5], conf=0.95), col = "red")
