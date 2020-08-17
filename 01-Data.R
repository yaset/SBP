rm(list = ls())
library(readxl)
library(dplyr)

## Load data

reboa_data <- read_excel("data/BD_REBOA_Toracotomia 780y1131_Revisores.xlsx")
names(reboa_data)
reboa_data <- reboa_data %>%
  rename(age = c1,
         penetrating_mechanism = c4,
         sbp = c6, ## systolic blood pressure
         hr = c7, ## heart rate
         rr = c8, ## respiratory rate
         glasgow_score = c9,
         iss = c19, ## ISS score
         tca = c42, ## traumatic cardiac arrest
         reboa = c43,
         survival_time = c196,
         sex = sexo2,
         type_trauma = mecanismo,
         mortality_28 = mortalidad_28d,
         mortality_7 = mortalidad_7d,
         mortalidad_24h = mortalidad_24
  )
str(reboa_data)

### Variables adjusted
reboa_data$penetrating_mechanism <- as.factor(reboa_data$penetrating_mechanism)
reboa_data$tca <- as.factor(reboa_data$tca)
reboa_data <- reboa_data %>%
  mutate(tca_binomial = ifelse(tca == "Paro",1,0))
reboa_data <- reboa_data %>%
  mutate(reboa_binomial = ifelse(reboa == "REBOA",1,0))
reboa_data <- reboa_data %>%
  mutate(penetrating = ifelse(type_trauma == "Penetrante",1,0))

write.csv(reboa_data, "data/database_reboa.csv")
