library(dplyr)
library(tidyr)
library(DescTools)
library(tidyr)
library(broom)

# Read the data
data <- read.delim("phyloseq/hillburns_metadata_Clean.tsv", stringsAsFactors = FALSE)

# Function to compute means and p-values for numerical variables
compute_numerical_results <- function(data, variable_name, variable) {
  # Impute missing values with the median
  data[[paste0(variable, "_imputed")]] <- ifelse(is.na(data[[variable]]), median(data[[variable]], na.rm = TRUE), data[[variable]])
  
  # Compute the means for each combination of the variable and 'alcohol_frequent'
  results <- data %>%
    group_by(alcohol_frequent) %>%
    summarise(mean = mean({{variable}}_imputed, na.rm = TRUE))
  
  # Perform a t-test for the variable between 'Yes' and 'No' alcohol_frequent
  t_test <- t.test({{variable}}_imputed ~ alcohol_frequent, data = data)
  
  # Add the t-test results to the results
  results <- results %>%
    pivot_wider(names_from = alcohol_frequent, values_from = mean) %>%
    mutate(
      p_value = t_test$p.value
    )
  
  # Add missing 'p_value' row
  results <- bind_rows(results, c("No" = NA, "Yes" = NA, p_value = t_test$p.value))
  
  # Print the variable name and results
  cat("\n", variable_name, "\n")
  print(results)
}

# Compute results for numerical variables
numerical_variables <- c("age", "height_in", "weight_lbs", "bmi", "latitude", "longitude", "coffee_amount")
numerical_variable_names <- c("Age", "Height (inches)", "Weight (lbs)", "BMI", "Latitude", "Longitude", "Coffee Amount")
for (i in seq_along(numerical_variables)) {
  compute_numerical_results(data, numerical_variable_names[i], numerical_variables[i])
}
