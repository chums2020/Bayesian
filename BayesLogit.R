#Bayesian Binomial (Logistic) Regression using heart data from package ncvreg
#Estimate marginal likelihood by Laplace approximation/MCMC and estimate the Bayes factor
#Reference: http://www.magesblog.com/2015/09/bayesian-regression-models-using-stan.html

library(ncvreg)
library(brms)
data(heart)

#need to install Rtools from http://cran.r-project.org/bin/windows/Rtools/
#assume the number of trials is constant across all observation 
bin.mod <- brm(chd| trials(10) ~ 
              sbp + tobacco +ldl + adiposity + famhist + typea + obesity + alcohol + age, 
              data = heart , family = "binomial")
