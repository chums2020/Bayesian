#Fit both Bayesian Binomial probit and logistic regression models
#Use the data 'heart' in package ncvreg. The response variable is chd. chd = 1: with coronary heart disease; chd = 0: not

library(ncvreg)
library(rstanarm)

data(heart)

#logit
fit1 <- stan_glm(chd ~sbp + tobacco +ldl + adiposity + famhist + typea + obesity + alcohol + age, data = heart , family = binomial(link = "logit"), prior = normal(0, 1), prior_intercept = normal(0, 1))

coef(fit1)
posterior_interval(fit1, prob = 0.95)

#probit
fit2<- stan_glm(chd ~sbp + tobacco +ldl + adiposity + famhist + typea + obesity + alcohol + age, data = heart , family = binomial(link = "probit"), prior = normal(0, 1), prior_intercept = normal(0, 1))

coef(fit2)
posterior_interval(fit2, prob = 0.95)

##
#Use MCMC to estimate marginal likelihood for each model and hence estimate Bayes factor


##
#Remove one data point from the original dataset. Fit the probit and logistic models again.
#Provide a HPD interval for the predicted value of the response using the predictor values
#for the one observation left out from the original dataset. This means you need to simulate
#from the posterior predictive distribution. Give a 95% HPD interval, as well as an estimate
#of the posterior predictive mean. Compare this estimate to the actual observed value of
#that response variable.

library(TeachingDemos)

heart2 <- heart[2:nrow(heart),]

#logit
fit1_omit <- stan_glm(chd ~sbp + tobacco +ldl + adiposity + famhist + typea + obesity + alcohol + age, 
                      data = heart2 , 
                      family = binomial(link = "logit"), 
                      prior = normal(0, 1), prior_intercept = normal(0, 1))

#simulate from posterior predictive distribution, using the predictor values for the omitted observation
y_rep_logit <- posterior_predict(fit1_omit, newdata = heart[1,])
table(y_rep_logit)
posterior_predictive_mean_logit <- sum(y_rep_logit)/length(y_rep_logit)
draw1 <- mcmc(y_rep_logit)

#probit
fit2_omit<- stan_glm(chd ~sbp + tobacco +ldl + adiposity + famhist + typea + obesity + alcohol + age,
                      data = heart2 , 
                      family = binomial(link = "probit"), 
                      prior = normal(0, 1), prior_intercept = normal(0, 1))

y_rep_probit <- posterior_predict(fit2_omit, newdata = heart[1,])
table(y_rep_probit)
posterior_predictive_mean_probit <- sum(y_rep_probit)/length(y_rep_probit)
draw2 <- mcmc(y_rep_probit)

#95% HPD interval (logit)
emp.hpd(draw1, conf = 0.95)

#posterior predictive mean (logit)
posterior_predictive_mean_logit

#95% HPD interval (probit)
emp.hpd(draw2, conf = 0.95)

#posterior predictive mean (probit)
posterior_predictive_mean_probit

#actual observed value
heart[1,10]
