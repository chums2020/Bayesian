#Compare Bayesian Binomial logit and probit regression models
#Using MCMCpack, we can compute the marginal likelihoods and Bayes factor directly by Laplace approximation.


library(MCMCpack)
library(ncvreg)
library(rstanarm)
data(heart)

posterior_logit <- MCMClogit(chd ~
                  sbp + tobacco +ldl + adiposity + famhist + typea + obesity + alcohol + age, 
                  data = heart, burnin = 1000, mcmc = 10000, b0 = 0, B0 =  0.1, 
                  user.prior.density=NULL, marginal.likelihood="Laplace")

posterior_probit <- MCMCprobit(chd ~
                   sbp + tobacco +ldl + adiposity + famhist + typea + obesity + alcohol + age, 
                   data = heart, burnin = 1000, mcmc = 10000, b0 = 0, B0 = 0.1, 
                   user.prior.density=NULL, marginal.likelihood="Laplace")


#Using MCMCpack, we produce a Bayes factor matrix. 
#The (1,2) element is the Bayes factor for logit model relative to probit model, which equals to 45.9. 
#Hence, we have a very strong evidence in favor of the logit model.(Jefferys scale, Kass \& Raftery)

BayesFactor(posterior_logit, posterior_probit )
