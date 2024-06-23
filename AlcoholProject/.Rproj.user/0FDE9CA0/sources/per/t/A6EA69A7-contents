# Load necessary libraries
library(dplyr)
library(readr)
library(mediation)

# Read in the diversity metrics
observed_features <- read_tsv("exported_alpha_diversity/observed_features.tsv")
chao1 <- read_tsv("exported_alpha_diversity/chao1.tsv")
shannon <- read_tsv("exported_alpha_diversity/shannon.tsv")
simpson <- read_tsv("exported_alpha_diversity/simpson.tsv")

# Read in metadata
metadata <- read_tsv("phyloseq/hillburns_metadata_Clean.tsv")

# Rename the 'ID' column to 'sample-id'
metadata <- metadata %>%
  rename("sample-id" = ID)

metadata <- metadata %>%
  filter(!(parkinson_disease == "Yes" & parkinson_disease_duration > 5))

# Combine alpha diversity metrics with metadata
combined_data <- metadata %>%
  left_join(observed_features, by = "sample-id") %>%
  left_join(chao1, by = "sample-id") %>%
  left_join(shannon, by = "sample-id") %>%
  left_join(simpson, by = "sample-id")

# Prepare data for mediation analysis
combined_data <- combined_data %>%
  mutate(parkinson_disease = ifelse(parkinson_disease == "Yes", 1, 0), 
         alcohol_frequent = ifelse(alcohol_frequent == "Yes", 1, 0))

# Function to check mediation significance
check_mediation <- function(mediation_result) {
  summary_result <- summary(mediation_result)
  acme_p_value <- summary_result$d0.p
  ade_p_value <- summary_result$z0.p
  
  if (acme_p_value < 0.05) {
    return("Mediation is significant")
  } else {
    return("No significant mediation")
  }
}

# Mediation model for observed
mediator_model <- lm(observed_features ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ observed_features + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_observed <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "observed_features", boot = TRUE, sims = 5000)
observed_mediation_result <- check_mediation(mediation_observed)

# Mediation model for shannon
mediator_model <- lm(shannon_entropy ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ shannon_entropy + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_shannon <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "shannon_entropy", boot = TRUE, sims = 5000)
shannon_mediation_result <- check_mediation(mediation_shannon)

# Mediation model for chao1
mediator_model <- lm(chao1 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ chao1 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_chao1 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "chao1", boot = TRUE, sims = 5000)
chao1_mediation_result <- check_mediation(mediation_chao1)

# Mediation model for simpson
mediator_model <- lm(simpson ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ simpson + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_simpson <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "simpson", boot = TRUE, sims = 5000)
simpson_mediation_result <- check_mediation(mediation_simpson)

# Print mediation results
print(paste("Observed features mediation result:", observed_mediation_result))
print(paste("Shannon mediation result:", shannon_mediation_result))
print(paste("Chao1 mediation result:", chao1_mediation_result))
print(paste("Simpson mediation result:", simpson_mediation_result))
print(summary(mediation_chao1))
