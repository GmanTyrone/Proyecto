# Load required libraries
library(tidyverse)
library(kableExtra)

# Read the TSV file
data <- read.delim("phyloseq/hillburns_metadata_Clean.tsv", stringsAsFactors = FALSE)

# Create a summary table
summary_table <- data %>%
  group_by(alcohol_frequent) %>%
  summarise(
    N = n(),
    mean_age = mean(age),
    sd_age = sd(age),
    female = sum(sex == "female"),
    male = sum(sex == "male"),
    white = sum(race == "White"),
    more_than_one_race = sum(race == "More than one race"),
    black_or_african_american = sum(race == "Black or African American"),
    mean_height = mean(height_in),
    sd_height = sd(height_in),
    mean_weight = mean(weight_lbs),
    sd_weight = sd(weight_lbs)
  )

# Generate the table
kable(summary_table, "html") %>%
  kable_styling(full_width = FALSE) %>%
  add_header_above(c(" ", "No" = 8, "Yes" = 8), bold = TRUE) %>%
  add_header_above(c(" ", "N", "Age", "", "Sex", "", "Race", "", "Height", "", "Weight", ""), bold = TRUE) %>%
  add_header_above(c("Alcohol frequency", "", "Mean (SD)", "Mean (SD)", "Female", "Male", "White", "More than one race", "Black or African American", "Mean (SD)", "Mean (SD)"), bold = TRUE)

# Export the table as an image
ggsave("summary_table.png", device = "png", width = 12, height = 6)
