#Assume the data have pareto distribution with parameter eta, and eta has a prior distribution of Gamma(α,β)
#Use the function rpareto in the VGAM package to generate a r.s. of size 100
#Construct a 95% Bayesian equal-tailed credible interval for eta, assume α = 2 and β = 0.5.

library(VGAM)
set.seed(1)
n <- 100
x <- rpareto(n, scale=1, shape=1)
alpha <- 2
beta <- 0.5

#compute the quantile function  
vals1 <-qgamma(p=seq(from =0.025, to = 0.975, length.out= 2), shape=n+alpha,scale =1/(sum(log(x))+beta))

#credible intervals
print(vals1)

#verify
pgamma(q=vals1[1], shape=n+alpha,scale =1/(sum(log(x))+beta))
pgamma(q=vals1[2], shape=n+alpha,scale =1/(sum(log(x))+beta))



#Find a 95% credible interval with length approximately 1, using any hyperparameter values. 
alpha<-400
beta <-0.01
  
quantile <-qgamma(p=seq(from =0.01, to = 0.99, length.out= 99), shape=n+alpha,scale =1/(sum(log(x))+beta))
post_density<-dgamma(x=quantile, shape=n+alpha,scale =1/(sum(log(x))+beta))
plot(quantile,post_density)

ends <-qgamma(p=seq(from =0.025, to = 0.975, length.out= 2), shape=n+alpha,scale =1/(sum(log(x))+beta))
print(ends)

#length of the 95% CI
print(ends[2]-ends[1])
