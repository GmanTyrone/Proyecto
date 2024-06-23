library(tidyverse)

# Read the CSV file
metadata <- read_csv("hillburns_metadata_O.csv")

# Convert all columns to character
metadata <- metadata %>%
  mutate(across(everything(), as.character))

# Select the desired columns
selected_columns <- c("ID", "age", "alcohol_amount", "bmi", "coffee_amount", 
                      "fruits_or_vegetables", "grains", "height_in", 
                      "latitude", "location", "longitude", "meats", 
                      "nuts", "other_neuro_problems", "parkinson_disease", 
                      "race", "sex", "smoke100", "ss_constipation", 
                      "weight_lbs", "yogurt")

# Filter the metadata to keep only the selected columns
clean_metadata <- metadata %>%
  select(all_of(selected_columns)) %>%
  # Filter rows where 'ss_constipation' is not NA
  filter(!is.na(ss_constipation))

# Replace "NA" values with empty fields
clean_metadata <- clean_metadata %>%
  mutate(across(everything(), na_if, "NA")) %>%
  mutate(across(everything(), replace_na, ""))

# Write the new data frame to a CSV file
write_csv(clean_metadata, "hillburns_metadata_Clean.csv")
