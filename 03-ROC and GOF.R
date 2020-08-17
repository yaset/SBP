rm(list = ls())

library(dplyr)
library(janitor)
library(ggplot2)
library(skimr)
library(pROC)
library(cutpointr)


reboa_database <- read.csv("data/database_reboa.csv")

names(reboa_database)

a1 <- roc(mortalidad_24h ~ sbp, data = reboa_database)

a2 <- cutpointr(data = reboa_database,
                class = mortalidad_24h,
                method = oc_youden_normal,
                x = sbp,
                direction = "<=")

a2 <- cutpointr(data = reboa_database,
                class = mortalidad_24h,
                method = maximize_metric,
                x = sbp,
                direction = "<=")


summary(a2)
plot(a2)


a3 <- cutpointr(data = reboa_database,
                class = tca_binomial,
                method = maximize_metric,
                x = sbp,
                direction = "<=")


summary(a3)
plot(a3)

reboa_database %>%
  skim(age)

reboa_database$sbp <- as.numeric(reboa_database$sbp)
table(reboa_database$reboa_binomial)

model1 <- glm(mortalidad_24h ~ sbp + penetrating + reboa_binomial,
              data = reboa_database,
              family = binomial(link = "logit"))

summary(model1)



model1 <- glm(tca_binomial ~ sbp + penetrating + reboa_binomial,
              data = reboa_database,
              family = binomial(link = "logit"))


lrtest(model1)
#### Goodness of fit
ResourceSelection::hoslem.test(reboa_database$tca_binomial, fitted(model1))

##### Test of individual predictors
summary(model1) ### P value is the same of Wald Test





library(pROC)
# Compute AUC for predicting Class with the variable CreditHistory.Critical
a1 <- roc(mortalidad_24h ~ sbp, data = reboa_database)
plot(a1, col="red")

library(ROCR)
library(caTools)
# Compute AUC for predicting Class with the model

model <- glm(mortalidad_24h ~ sbp + penetrating + reboa_binomial, family = binomial,
             data = reboa_database)

summary(model)
predict <- predict(model, type = "response")
table(reboa_database$mortalidad_24h, predict > 0.5)


ROCpred <- prediction(predict, reboa_database$mortalidad_24h)
ROCperf <- performance(ROCpred, 'tpr', 'fpr')
plot(ROCperf, colorize = TRUE, text.adj = c(-0.2,1.7))

auc <- performance(ROCpred, measure = "auc")
auc@y.values[[1]]



##### Mortalidad 24 hours

model1 <- glm(tca_binomial ~ sbp + penetrating + reboa_binomial,
              data = reboa_database,
              family = binomial(link = "logit"))


lrtest(model1)
#### Goodness of fit
ResourceSelection::hoslem.test(reboa_database$tca_binomial, fitted(model1))

##### Test of individual predictors
summary(model1) ### P value is the same of Wald Test





library(pROC)
# Compute AUC for predicting Class with the variable CreditHistory.Critical
a1 <- roc(mortalidad_24h ~ sbp, data = reboa_database)
plot(a1, col="red")

library(ROCR)
library(caTools)
# Compute AUC for predicting Class with the model

model <- glm(mortalidad_24h ~ sbp + penetrating + reboa_binomial + tca_binomial, family = binomial,
             data = reboa_database)

summary(model)

#### Goodness of fit
ResourceSelection::hoslem.test(reboa_database$mortalidad_24h, fitted(model))



predict <- predict(model, type = "response")
table(reboa_database$mortalidad_24h, predict > 0.5)


ROCpred <- prediction(predict, reboa_database$mortalidad_24h)
ROCperf <- performance(ROCpred, 'tpr', 'fpr')
plot(ROCperf, colorize = TRUE, text.adj = c(-0.2,1.7))

auc <- performance(ROCpred, measure = "auc")
auc@y.values[[1]]
