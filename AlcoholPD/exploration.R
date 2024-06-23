library(dplyr)
library(tidyverse)

metadata <- read.csv("hillburns_metadata_O.csv")
filtered_files <- read.table("filtered_files.txt", header = FALSE, col.names = "ID", stringsAsFactors = FALSE)
metadata <- metadata %>%
  filter(ID %in% filtered_files$ID)

# Remove rows that have NA in 'alcohol_how_much'
metadata <- metadata %>%
  filter(!is.na(alcohol_how_much))

# Remove rows that have something in the exclude column
metadata <- metadata %>%
  filter(is.na(excluded))

# Remove rows that have antiobiotics information
metadata <- metadata %>%
  filter(antibiotics_bool == "No")

# Remove rows that have pd information
metadata <- metadata %>%
  filter(parkinson_disease == "Yes")

metadata <- metadata %>%
  select(ID, age, alcohol_amount, alcohol_how_much, bmi, coffee_amount, fruits_or_vegetables,
         grains, height_in, latitude, location, longitude, meats, nuts,
         other_neuro_problems, parkinson_disease, race, sex, smoke100,
         ss_constipation, weight_lbs, yogurt)

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
    alcohol_frequent = ifelse(alcohol_how_much %in% c("None"), "No", "Yes")
  )

# Write the new data frame to a CSV file
write_csv(metadata, "hillburns_metadata_Clean.csv")

# #Unique values in alcohol_how_much
# unique_values <- unique(metadata$alcohol_how_much)
# 
# # Print the unique values
# print(unique_values)

# Count of unique values
count_unique_values <- table(metadata$alcohol_frequent)

# Print the count of unique values
print(count_unique_values)