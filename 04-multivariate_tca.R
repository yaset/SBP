rm(list = ls())

library(dplyr)
library(caTools)
library(rstanarm)
library(bayesplot)
library(janitor)
library(aod)
library(loo)
library(ggplot2)
library(cowplot)
library(MASS)

reboa_database <- read.csv("data/database_reboa.csv")
###### MODELS
names(reboa_database)
reboa_database <- reboa_database %>%
  mutate(sbp_i = case_when(
    sbp == 0 ~ 1,
    sbp > 0 & sbp <= 40 ~ 2,
    sbp > 40 & sbp <= 60 ~ 3,
    sbp > 60 & sbp <= 80 ~ 4,
    sbp > 80 & sbp <= 100 ~ 5,
    sbp > 100 & sbp <= 120 ~ 6,
    sbp > 120 ~ 7
  ))
reboa_database$sbp_i <- as.factor(reboa_database$sbp_i)
levels(reboa_database$sbp_i) <- c("0","1-40","41-60","61-80","81-100",
                                  "101-120",">120")
reboa_database %>%
  tabyl(sbp_i)

#########################################################
##########   MORTALIDAD 24 H ############################


##### SBP continous
# Uniform
model1 <- glm(mortalidad_24h ~ tca_binomial, 
            data = reboa_database, 
            family = poisson(link = "log"))

vcov(model1)
model1
exp(coef(model1))
confint(model1)


foca1 <- glm(mortalidad_24h ~ tca_binomial, 
             data = reboa_database, 
             family = poisson(link = "log"))

summary(foca1)

ci95 <- confint(foca1)
coef(foca1)

data.frame(RR = coef(foca1), ci95)
data.frame(RR = exp(coef(foca1)), exp(ci95))

summary(foca1)
summ()






fit1 <- stan_glm(mortalidad_24h ~ sbp_i, 
                 data = reboa_database, 
                 family = binomial(link = "log"),
                 prior_intercept = normal(0,5),
                 prior = normal(0,5),
                 chains = 4, iter = 10000,
                 seed = 12345)

summary(fit1)
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

wald.test(b = coef(fit1), Sigma = vcov(fit1), Terms = 2)
loo(fit1)





fit1$model
exp(fit1$linear.predictors)

analysis <- data.frame(sbp_i = fit1$model$sbp, morta = exp(fit1$linear.predictors))

plot(analysis$sbp_i,analysis$morta)

analysis <- analysis %>%
  mutate(sbp_i2 = case_when(
    sbp_i <= 40 ~ 1,
    sbp_i > 40 & sbp_i <= 60 ~ 2,
    sbp_i > 60 & sbp_i <= 80 ~ 3,
    sbp_i > 80 & sbp_i <= 100 ~ 4,
    sbp_i > 100 & sbp_i <= 120 ~ 5,
    sbp_i > 120 ~ 6,
  ))

##### RR

# Uniform

fit2 <- stan_glm(mortalidad_24h ~ sbp_i + tca_binomial + reboa_binomial + penetrating, 
                 data = reboa_database, 
                 family = binomial(link = "log"),
                 prior_intercept = normal(0,5),
                 prior = normal(0,5),
                 chains = 4, iter = 10000,
                 seed = 12345)

wald.test(b = coef(fit1), Sigma = vcov(fit1), Terms = 2)
loo(fit1)

ci95 <- posterior_interval(fit1, prob = 0.90)

data2 <- data.frame(RR = exp(coefficients(fit1)), exp(ci95))
data2 <- cbind(names = row.names(data2), data2)


data2 <- data2[2:7,]
colnames(data2)[3] <- "low"
colnames(data2)[4] <- "high"
data2$names <- as.character(data2$names)


ci95 <- posterior_interval(fit2, prob = 0.90)

data3 <- data.frame(RR = exp(coefficients(fit2)), exp(ci95))
data3 <- cbind(names = row.names(data3), data3)


data3 <- data3[2:7,]
colnames(data3)[3] <- "low"
colnames(data3)[4] <- "high"
data3$names <- as.character(data3$names)

data2
data3


ggplot(data = data2, aes(x = names, y = RR))+
  geom_errorbar(aes(x = names, ymin = low , ymax = high),width = 0.5, color = "Black", size = 1, inherit.aes = FALSE)+
  geom_point( aes(x = names, y = RR), size = 2, color = "Black", stroke = 1.5, inherit.aes = FALSE)+
  
  geom_errorbar(data = data3, aes(x = names, ymin = low , ymax = high),width = 0.5, color = "Blue", size = 1, inherit.aes = FALSE)+
  geom_point(data = data3, aes(x = names, y = RR),size = 2, color = "Blue", stroke = 1.5, inherit.aes = FALSE)+
  
  
  scale_x_discrete("Systolic Blood Pressure mmHg",labels = c("1-40","41-60","61-80","81-100",
                                                             "101-120",">120"))+
  scale_y_continuous("Risk Ratio of Traumatic Cardiac Arrest", breaks = seq(0,1.4,by = 0.2),
                     limits = c(0,1.6))+
  geom_hline(yintercept = 1, linetype = 2)+
  theme_cowplot()
  
  






fit1$model
exp(fit1$linear.predictors)

analysis <- data.frame(sbp_i = fit1$model$sbp, morta = exp(fit1$linear.predictors))

plot(analysis$sbp_i,analysis$morta)

analysis <- analysis %>%
  mutate(sbp_i2 = case_when(
    sbp_i <= 40 ~ 1,
    sbp_i > 40 & sbp_i <= 60 ~ 2,
    sbp_i > 60 & sbp_i <= 80 ~ 3,
    sbp_i > 80 & sbp_i <= 100 ~ 4,
    sbp_i > 100 & sbp_i <= 120 ~ 5,
    sbp_i > 120 ~ 6,
  ))



analysis %>%
  select(sbp_i2, morta) %>%
  group_by(sbp_i2) %>%
  skimr::skim(morta)


plot(analysis$sbp_i2, analysis$morta)

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
