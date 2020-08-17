rm(list = ls())

library(dplyr)
library(caTools)
library(rstanarm)
library(bayesplot)
library(janitor)
library(aod)
library(loo)
library(ggplot2)

reboa_database <- read.csv("data/database_reboa.csv")
###### MODELS
names(reboa_database)
reboa_database <- reboa_database %>%
  mutate(sbp_i = case_when(
    sbp == 0 ~ 1,
    sbp > 0 & sbp <= 60 ~ 2,
    sbp > 60 & sbp <= 90 ~ 3,
    sbp > 90 & sbp <= 120 ~ 4,
    sbp > 120 ~ 5
  ))
reboa_database$sbp_i <- as.factor(reboa_database$sbp_i)
levels(reboa_database$sbp_i) <- c("0","1-60","61-90",">90")
reboa_database %>%
  tabyl(sbp_i)

#########################################################
##########   MORTALIDAD 24 H ############################


##### SBP continous
# Uniform
reboa_database2 <- reboa_database %>%
  filter(penetrating == 0)


fit1 <- stan_glm(mortalidad_24h ~ sbp + tca_binomial + reboa_binomial + penetrating, 
                 data = reboa_database, 
                 family = binomial(link = "log"),
                 prior_intercept = normal(0,5),
                 prior = normal(0,5),
                 chains = 4, iter = 10000,
                 seed = 12345)

summary(fit1, digits = 3)
wald.test(b = coef(fit1), Sigma = vcov(fit1), Terms = 2)
loo(fit1)

ci95 <- posterior_interval(fit1, prob = 0.90)
data.frame(RR = coefficients(fit1), ci95)
data.frame(RR = exp(coefficients(fit1)), exp(ci95))


fit1 <- stan_glm(mortalidad_24h ~ sbp_i + tca_binomial + reboa_binomial + penetrating, 
                 data = reboa_database, 
                 family = binomial(link = "log"),
                 prior_intercept = normal(0,5),
                 prior = normal(0,5),
                 chains = 4, iter = 10000,
                 seed = 12345)

wald.test(b = coef(fit1), Sigma = vcov(fit1), Terms = 5)
loo(fit1)

ci95 <- posterior_interval(fit1, prob = 0.95)
data.frame(RR = coefficients(fit1), ci95)
data.frame(RR = exp(coefficients(fit1)), exp(ci95))




mcmc_areas(fit1)



fit1 <- stan_glm(mortalidad_24h ~ sbp_i + tca_binomial, 
                 data = reboa_database, 
                 family = binomial(link = "log"),
                 prior_intercept = normal(0,5),
                 prior = normal(0,5),
                 chains = 4, iter = 10000,
                 seed = 12345)


summary(fit1, probs = c(0.025, 0.975), digits = 3)


fit1 <- stan_glm(mortalidad_24h ~ sbp + reboa_binomial, 
                 data = reboa_database, 
                 family = binomial(link = "log"),
                 prior_intercept = normal(0,5),
                 prior = normal(0,5),
                 chains = 4, iter = 10000,
                 seed = 12345)


summary(fit1, probs = c(0.025, 0.975), digits = 3)
wald.test(b = coef(fit1), Sigma = vcov(fit1), Terms = 2)

fit1 <- stan_glm(mortalidad_24h ~ sbp + penetrating, 
                 data = reboa_database, 
                 family = binomial(link = "log"),
                 prior_intercept = normal(0,5),
                 prior = normal(0,5),
                 chains = 4, iter = 10000,
                 seed = 12345)


summary(fit1, probs = c(0.025, 0.975), digits = 3)





















# Normal (0,5)
fit1 <- stan_glmer(mortalidad_24h ~  sbp_i + (1  | tca_binomial), 
                 data = reboa_database, 
                 family = binomial(link = "log"),
                 chains = 4, iter = 10000,
                 seed = 12345)

wald.test(b = coef(fit1), Sigma = vcov(fit1), Terms = 2)
loo(fit1)
### Selected
summary(fit1, probs = c(0.025, 0.975), digits = 3)


ci95 <- posterior_interval(fit1, prob = 0.95)
data.frame(RR = fit1$coefficients, ci95[1:6,])
data.frame(RR = exp(fit1$coefficients), exp(ci95[1:6,]))






# Cauchy (0,3)
fit1 <- stan_glm(mortalidad_24h ~ sbp, 
                 data = reboa_database, 
                 family = binomial(link = "log"),
                 prior_intercept = cauchy(0,3),
                 prior = cauchy(0,3),
                 chains = 4, iter = 10000)

wald.test(b = coef(fit1), Sigma = vcov(fit1), Terms = 2)
loo(fit1)

summary(fit1, probs = c(0.025, 0.975), digits = 3)
prior_summary(fit1)








wald.test(b = coef(fit1), Sigma = vcov(fit1), Terms = 2)



color_scheme_set("viridis")
mcmc_trace(fit1, n_warmup = 1000)
mcmc_hist(fit1)


loo1 <- loo(fit1)


#### Normal (0,3)

fit_normal <- stan_glm(mortalidad_24h ~ sbp_i + tca_binomial, 
                       data = reboa_database, 
                       family = binomial(link = "log"),
                       prior_intercept = normal(0,3),
                       prior = normal(0,3),
                       chains = 4, iter = 10000)

summary(fit_normal, probs = c(0.025, 0.975), digits = 3)
pplot <- plot(fit_normal, "areas", prob_outer = 1, prob = 0.95)
pplot + geom_vline(xintercept = 0)

exp(fit_normal$coefficients)

ci95 <- exp(posterior_interval(fit_normal, prob = 0.95))
c(exp(coef(fit_normal)[2]),exp(ci95))
wald.test(b = coef(fit_normal), Sigma = vcov(fit_normal), Terms = 5)

color_scheme_set("viridis")
mcmc_trace(fit1, n_warmup = 1000)
mcmc_hist(fit1)


plot(reboa_database$sbp, fit_normal$residuals)
lines(reboa_database$sbp, fit_normal$fitted.values)

loo2 <- loo(fit_normal)
loo_compare(loo1,loo2)

summary(fit1$residuals)
summary(fit_normal$residuals)


posterior_vs_prior(fit_normal)

#### auc
predpr <- predict(fit_normal, type = c("response"))
library(pROC)

roccurve <- roc(reboa_database$tca_binomial ~ reboa_database$sbp)
plot(roccurve)
auc(roccurve)


# Predicted probabilities
linpred <- posterior_linpred(fit_normal)
preds <- posterior_linpred(fit_normal, transform=TRUE)
pred <- colMeans(preds)
pr <- as.integer(pred >= 0.5)

# posterior classification accuracy
round(mean(xor(pr,as.integer(sbp==0))),2)


cv_vars
print(fit1, digits = 3)


coef(fit1)
vcov(fit1)

plot(reboa_database$sbp, fit1$linear.predictors)
loo(fit1)

model1$
  fit1$
  fit1$coefficients
summary(fit1$residuals)
summary(model1)
mcmc_areas(fit1,
           pars = c("sbp"))

plot(fit1)
BIC(fit1)



prior_summary(fit1)
fit1$coefficients
hist(fit1$fitted.values)
fit1$model
fit1$stan_function

pred




##### intervalos de confianza



plot(fit1)
launch_shinystan(fit1, ppd = FALSE)
y_rep <- posterior_predict(fit1) #### matrix de los datos
dim(y_rep)

############# LOO
loo1 <- loo(fit1)
plot(loo1)

pairs(fit1)
