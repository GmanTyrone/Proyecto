library(dplyr)
library(tidyverse)

metadata <- read.csv("hillburns_metadata_O.csv")
filtered_files <- read.table("filtered_files.txt", header = FALSE, col.names = "ID", stringsAsFactors = FALSE)
metadata <- metadata %>%
  filter(ID %in% filtered_files$ID)
# 
# # Remove rows that have NA in 'alcohol_how_much'
# metadata <- metadata %>%
#   filter(!is.na(alcohol_how_much))

# Remove rows that have something in the exclude column
metadata <- metadata %>%
  filter(is.na(excluded))

# Remove rows that have antiobiotics information
metadata <- metadata %>%
  filter(antibiotics_bool == "No")

# metadata <- metadata %>%
#   filter(!(parkinson_disease == "Yes" & parkinson_disease_duration > 5))


metadata <- metadata %>%
  select(ID, age, alcohol_amount, alcohol_how_much, anxiety, bmi, cesarean, coffee_amount, fruits_or_vegetables,
         grains, height_in, impulse, insomnia, latitude, location, longitude, loss_10lbs, meats, medication_digestive_problems_bool, nuts,
         obsessive, other_neuro_problems, parkinson_disease, parkinson_disease_duration, race, sex, sleep_aid, smoke100,
         ss_bloating, ss_constipation, ss_excess_gas, ssdigest_prob, weight_lbs, yogurt)

# Update the values in selected columns
metadata <- metadata %>%
  mutate(
    fruits_or_vegetables = case_when(
      fruits_or_vegetables %in% c("Less than once a month or never", "Few times a month", "Few times a week") ~ "No",
      fruits_or_vegetables == "At least once a day" ~ "Yes",
      TRUE ~ fruits_or_vegetables
    ),
    grains = case_when(
      grains %in% c("Less than once a month or never", "Few times a month", "Few times a week") ~ "No",
      grains == "At least once a day" ~ "Yes",
      TRUE ~ grains
    ),
    meats = case_when(
      meats %in% c("Less than once a month or never", "Few times a month", "Few times a week") ~ "No",
      meats == "At least once a day" ~ "Yes",
      TRUE ~ meats
    ),
    medication_digestive_problems_bool = case_when(
      medication_digestive_problems_bool == "Don't Know" ~ NA_character_,
      TRUE ~ medication_digestive_problems_bool
      ),
    nuts = case_when(
      nuts %in% c("Less than once a month or never", "Few times a month", "Few times a week") ~ "No",
      nuts == "At least once a day" ~ "Yes",
      TRUE ~ nuts
    ),
    yogurt = case_when(
      yogurt %in% c("Less than once a month or never", "Few times a month", "Few times a week") ~ "No",
      yogurt == "At least once a day" ~ "Yes",
      TRUE ~ yogurt
    ),
    sleep_aid = case_when(
      sleep_aid %in% c("Sometimes", "Regularly") ~ "Yes",
      sleep_aid == "Missing" ~ NA_character_,
      TRUE ~ sleep_aid
    ),
    ss_bloating = case_when(
      ss_bloating == "Missing" ~ NA_character_,
      TRUE ~ ss_bloating
    ),
    ssdigest_prob = case_when(
      ssdigest_prob == "Missing" ~ NA_character_,
      TRUE ~ ssdigest_prob
    ),
    alcohol_frequent = ifelse(alcohol_how_much %in% c("None"), "No", "Yes")
  )

# Write the new data frame to a CSV file
write_csv(metadata, "hillburns_metadata_Clean.csv")

#Unique values in alcohol_how_much
unique_values <- unique(metadata$alcohol_how_much)

# Print the unique values
print(unique_values)

# Count of unique values
count_unique_values <- table(metadata$parkinson_disease)

# Print the count of unique values
print(count_unique_values)