library(dplyr)
library(tidyr)
library(DescTools)
library(tidyr)
library(broom)

# Read the data
data <- read.delim("hillburns_metadata_Clean.tsv", stringsAsFactors = FALSE)

#AGE
data$age_imputed <- ifelse(is.na(data$age), median(data$age, na.rm = TRUE), data$age)
age_results <- data %>%
  group_by(alcohol_frequent) %>%
  summarise(mean_age = mean(age_imputed, na.rm = TRUE))
t_test <- t.test(age_imputed ~ alcohol_frequent, data = data)
age_results <- age_results %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean_age) %>%
  mutate(
    p_value = t_test$p.value
  )
age_results <- bind_rows(age_results, c("No" = NA, "Yes" = NA, p_value = t_test$p.value))
print("age")
print(age_results)


# SeX
sex_results <- data %>%
  group_by(sex, alcohol_frequent) %>%
  summarise(mean = mean(age, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean)
p_value <- chisq.test(table(data$sex, data$alcohol_frequent))$p.value
sex_results <- bind_rows(sex_results, c("p_value" = p_value))
print(sex_results)

# Race
race_results <- data %>%
  group_by(race, alcohol_frequent) %>%
  summarise(mean = mean(age, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean)
p_value <- chisq.test(table(data$race, data$alcohol_frequent))$p.value
race_results <- bind_rows(race_results, c("p_value" = p_value))
print(race_results)

#HEIGHT_IN
data$height_in_imputed <- ifelse(is.na(data$height_in), median(data$height_in, na.rm = TRUE), data$height_in)
height_in_results <- data %>%
  group_by(alcohol_frequent) %>%
  summarise(mean_height_in = mean(height_in_imputed, na.rm = TRUE))
t_test <- t.test(height_in_imputed ~ alcohol_frequent, data = data)
height_in_results <- height_in_results %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean_height_in) %>%
  mutate(
    p_value = t_test$p.value
  )
height_in_results <- bind_rows(height_in_results, c("No" = NA, "Yes" = NA, p_value = t_test$p.value))
print ("height_in")
print(height_in_results)

#WEIGHT_LBS
data$weight_lbs_imputed <- ifelse(is.na(data$weight_lbs), median(data$weight_lbs, na.rm = TRUE), data$weight_lbs)
weight_lbs_results <- data %>%
  group_by(alcohol_frequent) %>%
  summarise(mean_weight_lbs = mean(weight_lbs_imputed, na.rm = TRUE))
t_test <- t.test(weight_lbs_imputed ~ alcohol_frequent, data = data)
weight_lbs_results <- weight_lbs_results %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean_weight_lbs) %>%
  mutate(
    p_value = t_test$p.value
  )
weight_lbs_results <- bind_rows(weight_lbs_results, c("No" = NA, "Yes" = NA, p_value = t_test$p.value))
print("weight_lbs")
print(weight_lbs_results)


#BMI
data$bmi_imputed <- ifelse(is.na(data$bmi), median(data$bmi, na.rm = TRUE), data$bmi)
bmi_results <- data %>%
  group_by(alcohol_frequent) %>%
  summarise(mean_bmi = mean(bmi_imputed, na.rm = TRUE))
t_test <- t.test(bmi_imputed ~ alcohol_frequent, data = data)
bmi_results <- bmi_results %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean_bmi) %>%
  mutate(
    p_value = t_test$p.value
  )
bmi_results <- bind_rows(bmi_results, c("No" = NA, "Yes" = NA, p_value = t_test$p.value))
print("BMI")
print(bmi_results)

#LATITUDE
data$latitude_imputed <- ifelse(is.na(data$latitude), median(data$latitude, na.rm = TRUE), data$latitude)
latitude_results <- data %>%
  group_by(alcohol_frequent) %>%
  summarise(mean_latitude = mean(latitude_imputed, na.rm = TRUE))
t_test <- t.test(latitude_imputed ~ alcohol_frequent, data = data)
latitude_results <- latitude_results %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean_latitude) %>%
  mutate(
    p_value = t_test$p.value
  )
latitude_results <- bind_rows(latitude_results, c("No" = NA, "Yes" = NA, p_value = t_test$p.value))
print("latitude")
print(latitude_results)

# Location
location_results <- data %>%
  group_by(location, alcohol_frequent) %>%
  summarise(mean = mean(age, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean)
p_value <- chisq.test(table(data$location, data$alcohol_frequent))$p.value
location_results <- bind_rows(location_results, c("p_value" = p_value))
print(location_results)

#LONGITUDE
data$longitude_imputed <- ifelse(is.na(data$longitude), median(data$longitude, na.rm = TRUE), data$longitude)
longitude_results <- data %>%
  group_by(alcohol_frequent) %>%
  summarise(mean_longitude = mean(longitude_imputed, na.rm = TRUE))
t_test <- t.test(longitude_imputed ~ alcohol_frequent, data = data)
longitude_results <- longitude_results %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean_longitude) %>%
  mutate(
    p_value = t_test$p.value
  )
longitude_results <- bind_rows(longitude_results, c("No" = NA, "Yes" = NA, p_value = t_test$p.value))
print("longitude")
print(longitude_results)

# Smoke100
smoke100_results <- data %>%
  group_by(smoke100, alcohol_frequent) %>%
  summarise(mean = mean(age, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean)
p_value <- chisq.test(table(data$smoke100, data$alcohol_frequent))$p.value
smoke100_results <- bind_rows(smoke100_results, c("p_value" = p_value))
print(smoke100_results)

#COFFEE_AMOUNT
data$coffee_amount_imputed <- ifelse(is.na(data$coffee_amount), median(data$coffee_amount, na.rm = TRUE), data$coffee_amount)
coffee_amount_results <- data %>%
  group_by(alcohol_frequent) %>%
  summarise(mean_coffee_amount = mean(coffee_amount_imputed, na.rm = TRUE))
t_test <- t.test(coffee_amount_imputed ~ alcohol_frequent, data = data)
coffee_amount_results <- coffee_amount_results %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean_coffee_amount) %>%
  mutate(
    p_value = t_test$p.value
  )
coffee_amount_results <- bind_rows(coffee_amount_results, c("No" = NA, "Yes" = NA, p_value = t_test$p.value))
print("coffee_amount")
print(coffee_amount_results)

# Fruits_or_vegetables
fruits_or_vegetables_results <- data %>%
  group_by(fruits_or_vegetables, alcohol_frequent) %>%
  summarise(mean = mean(age, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean)
p_value <- chisq.test(table(data$fruits_or_vegetables, data$alcohol_frequent))$p.value
fruits_or_vegetables_results <- bind_rows(fruits_or_vegetables_results, c("p_value" = p_value))
print(fruits_or_vegetables_results)

# Grains
grains_results <- data %>%
  group_by(grains, alcohol_frequent) %>%
  summarise(mean = mean(age, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean)
p_value <- chisq.test(table(data$grains, data$alcohol_frequent))$p.value
grains_results <- bind_rows(grains_results, c("p_value" = p_value))
print(grains_results)

# Meats
meats_results <- data %>%
  group_by(meats, alcohol_frequent) %>%
  summarise(mean = mean(age, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean)
p_value <- chisq.test(table(data$meats, data$alcohol_frequent))$p.value
meats_results <- bind_rows(meats_results, c("p_value" = p_value))
print(meats_results)

# Nuts
nuts_results <- data %>%
  group_by(nuts, alcohol_frequent) %>%
  summarise(mean = mean(age, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean)
p_value <- chisq.test(table(data$nuts, data$alcohol_frequent))$p.value
nuts_results <- bind_rows(nuts_results, c("p_value" = p_value))
print(nuts_results)

# Yogurt
yogurt_results <- data %>%
  group_by(yogurt, alcohol_frequent) %>%
  summarise(mean = mean(age, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean)
p_value <- chisq.test(table(data$yogurt, data$alcohol_frequent))$p.value
yogurt_results <- bind_rows(yogurt_results, c("p_value" = p_value))
print(yogurt_results)

# Other_neuro_problems
other_neuro_problems_results <- data %>%
  group_by(other_neuro_problems, alcohol_frequent) %>%
  summarise(mean = mean(age, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = alcohol_frequent, values_from = mean)
p_value <- chisq.test(table(data$other_neuro_problems, data$alcohol_frequent))$p.value
other_neuro_problems_results <- bind_rows(other_neuro_problems_results, c("p_value" = p_value))
print(other_neuro_problems_results)

