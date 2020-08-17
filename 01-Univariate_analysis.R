rm(list = ls())

library(dplyr)
library(caTools)
library(rstanarm)
library(bayesplot)
library(janitor)
library(aod)
library(loo)
library(ggplot2)
library(skimr)

reboa_database <- read.csv("data/database_reboa.csv")
###### MODELS
names(reboa_database)
reboa_database <- reboa_database %>%
  mutate(sbp_i = case_when(
    sbp == 0 ~ 1,
    sbp > 0 & sbp <= 60 ~ 2,
    sbp > 60 & sbp <= 90 ~ 3,
    sbp > 90  ~ 4))
reboa_database$sbp_i <- as.factor(reboa_database$sbp_i)
levels(reboa_database$sbp_i) <- c("0","1-60","61-90",">90")
reboa_database %>%
  tabyl(sbp_i)

### ISS
summary(reboa_database$iss)
reboa_database <- reboa_database %>%
  mutate(iss_i = case_when(
    iss >= 1 & iss <= 25 ~ 1,
    iss >= 26 & iss <= 34 ~ 2,
    iss >= 35 ~ 3
  ))

reboa_database$iss_i[is.na(reboa_database$iss)] <- 2
reboa_database %>%
  tabyl(iss_i)
reboa_database$iss_i <- as.factor(reboa_database$iss_i)




#########################################################
##########   MORTALIDAD 24 H ############################


##### SBP continous
# Uniform
reboa_database <- reboa_database %>%
  mutate(ntca_binomial = ifelse(tca_binomial == 1,0,1))

fit1 <- stan_glm(tca_binomial ~ sbp, 
                 data = reboa_database, 
                 family = binomial(link = "logit"),
                 prior_intercept = normal(0,5),
                 prior = normal(0,5),
                 chains = 4, iter = 10000,
                 seed = 12345)

loo(fit1)
summary(fit1, digits = 3,probs = c(0.025, 0.975))
wald.test(b = coef(fit1), Sigma = vcov(fit1), Terms = 2)

exp(fit1$coefficients)



### Mortalidad predictivea

# Normal (0,5)
fit1 <- stan_glm(tca_binomial ~ sbp, 
                 data = reboa_database, 
                 family = binomial(link = "logit"),
                 prior_intercept = normal(0,5),
                 prior = normal(0,5),
                 chains = 4, iter = 10000)

wald.test(b = coef(fit1), Sigma = vcov(fit1), Terms = 2)
loo(fit1)
base <- data.frame(sbp = reboa_database$sbp_i,lin = exp(fit1$linear.predictors))
base %>%
  group_by(sbp) %>%
  skim_without_charts(lin)

# Cauchy (0,3)
fit1 <- stan_glm(tca_binomial ~ sbp, 
                 data = reboa_database, 
                 family = binomial(link = "logit"),
                 prior_intercept = cauchy(0,3),
                 prior = cauchy(0,3),
                 chains = 4, iter = 10000)

loo(fit1)


##### REBOA
names(reboa_database)
# Uniform
fit1 <- stan_glm(tca_binomial ~ reboa_binomial, 
                 data = reboa_database, 
                 family = binomial(link = "logit"),
                 prior_intercept = NULL,
                 prior = NULL,
                 chains = 4, iter = 10000,
                 seed = 12345)


loo(fit1)

summary(fit1, probs = c(0.025, 0.975), digits = 3)

# Normal (0,5)
fit1 <- stan_glm(tca_binomial ~ reboa_binomial, 
                 data = reboa_database, 
                 family = binomial(link = "logit"),
                 prior_intercept = normal(0,5),
                 prior = normal(0,5),
                 chains = 4, iter = 10000,
                 seed = 12345)

loo(fit1)
### Selected
summary(fit1, probs = c(0.025, 0.975), digits = 3)

ci95 <- posterior_interval(fit1, prob = 0.95)
data.frame(RR = coefficients(fit1), ci95)
data.frame(RR = exp(coefficients(fit1)), exp(ci95))

prior_summary(fit1)

# Cauchy (0,3)
fit1 <- stan_glm(mortalidad_24h ~ reboa_binomial, 
                 data = reboa_database, 
                 family = binomial(link = "logit"),
                 prior_intercept = cauchy(0,3),
                 prior = cauchy(0,3),
                 chains = 4, iter = 10000,
                 seed = 12345)


##### Mechanims of Trauma
# Uniform
fit1 <- stan_glm(tca_binomial ~ penetrating, 
                 data = reboa_database, 
                 family = binomial(link = "logit"),
                 prior_intercept = NULL,
                 prior = NULL,
                 chains = 4, iter = 10000,
                 seed = 12345)

#### Evaluation of each P
loo(fit1)

summary(fit1, probs = c(0.025, 0.975), digits = 3)

# Normal (0,5)


fit1 <- stan_glm(tca_binomial ~ penetrating, 
                 data = reboa_database, 
                 family = binomial(link = "logit"),
                 prior_intercept = normal(0,5),
                 prior = normal(0,5),
                 chains = 4, iter = 10000,
                 seed = 12345)

loo(fit1)
wald.test(b = coef(fit1), Sigma = vcov(fit1), Terms = 2)


### Selected
summary(fit1, probs = c(0.025, 0.975), digits = 3)

ci95 <- posterior_interval(fit1, prob = 0.95)
data.frame(RR = coefficients(fit1), ci95)
data.frame(RR = exp(coefficients(fit1)), exp(ci95))

prior_summary(fit1)

# Cauchy (0,3)
fit1 <- stan_glm(mortalidad_24h ~ penetrating, 
                 data = reboa_database, 
                 family = binomial(link = "log"),
                 prior_intercept = cauchy(0,3),
                 prior = cauchy(0,3),
                 chains = 4, iter = 10000,
                 seed = 12345)

wald.test(b = coef(fit1), Sigma = vcov(fit1), Terms = 2)
loo(fit1)
mcmc_areas(fit1)

