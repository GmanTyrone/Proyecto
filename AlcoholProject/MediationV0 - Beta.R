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
  rename(`sample-id` = ID)

# Combine alpha diversity metrics with metadata
combined_data <- metadata %>%
  left_join(observed_features, by = "sample-id") %>%
  left_join(chao1, by = "sample-id") %>%
  left_join(shannon, by = "sample-id") %>%
  left_join(simpson, by = "sample-id")

# Read in PCoA coordinates (assuming files are in a specific format)
unweighted_unifrac_pcoa <- read_tsv("exported_beta_diversity/unweighted_unifrac_pcoa.tsv")
weighted_unifrac_pcoa <- read_tsv("exported_beta_diversity/weighted_unifrac_pcoa.tsv")
canberra_pcoa <- read_tsv("exported_beta_diversity/canberra_pcoa.tsv")

# Extract the first five principal coordinates (PC1 to PC5) from each PCoA table

# For unweighted Unifrac
unweighted_unifrac_PC1 <- unweighted_unifrac_pcoa %>%
  dplyr::select(`sample-id`, PC1) %>%
  rename(unweighted_unifrac_PC1 = PC1)

unweighted_unifrac_PC2 <- unweighted_unifrac_pcoa %>%
  dplyr::select(`sample-id`, PC2) %>%
  rename(unweighted_unifrac_PC2 = PC2)

unweighted_unifrac_PC3 <- unweighted_unifrac_pcoa %>%
  dplyr::select(`sample-id`, PC3) %>%
  rename(unweighted_unifrac_PC3 = PC3)

unweighted_unifrac_PC4 <- unweighted_unifrac_pcoa %>%
  dplyr::select(`sample-id`, PC4) %>%
  rename(unweighted_unifrac_PC4 = PC4)

unweighted_unifrac_PC5 <- unweighted_unifrac_pcoa %>%
  dplyr::select(`sample-id`, PC5) %>%
  rename(unweighted_unifrac_PC5 = PC5)

# For weighted Unifrac
weighted_unifrac_PC1 <- weighted_unifrac_pcoa %>%
  dplyr::select(`sample-id`, PC1) %>%
  rename(weighted_unifrac_PC1 = PC1)

weighted_unifrac_PC2 <- weighted_unifrac_pcoa %>%
  dplyr::select(`sample-id`, PC2) %>%
  rename(weighted_unifrac_PC2 = PC2)

weighted_unifrac_PC3 <- weighted_unifrac_pcoa %>%
  dplyr::select(`sample-id`, PC3) %>%
  rename(weighted_unifrac_PC3 = PC3)

weighted_unifrac_PC4 <- weighted_unifrac_pcoa %>%
  dplyr::select(`sample-id`, PC4) %>%
  rename(weighted_unifrac_PC4 = PC4)

weighted_unifrac_PC5 <- weighted_unifrac_pcoa %>%
  dplyr::select(`sample-id`, PC5) %>%
  rename(weighted_unifrac_PC5 = PC5)

# For Canberra
canberra_PC1 <- canberra_pcoa %>%
  dplyr::select(`sample-id`, PC1) %>%
  rename(canberra_PC1 = PC1)

canberra_PC2 <- canberra_pcoa %>%
  dplyr::select(`sample-id`, PC2) %>%
  rename(canberra_PC2 = PC2)

canberra_PC3 <- canberra_pcoa %>%
  dplyr::select(`sample-id`, PC3) %>%
  rename(canberra_PC3 = PC3)

canberra_PC4 <- canberra_pcoa %>%
  dplyr::select(`sample-id`, PC4) %>%
  rename(canberra_PC4 = PC4)

canberra_PC5 <- canberra_pcoa %>%
  dplyr::select(`sample-id`, PC5) %>%
  rename(canberra_PC5 = PC5)


# Combine beta diversity PC1s with metadata
# For unweighted_unifrac
combined_data <- combined_data %>%
  left_join(unweighted_unifrac_PC1, by = "sample-id") %>%
  left_join(unweighted_unifrac_PC2, by = "sample-id") %>%
  left_join(unweighted_unifrac_PC3, by = "sample-id") %>%
  left_join(unweighted_unifrac_PC4, by = "sample-id") %>%
  left_join(unweighted_unifrac_PC5, by = "sample-id")

# For weighted_unifrac
combined_data <- combined_data %>%
  left_join(weighted_unifrac_PC1, by = "sample-id") %>%
  left_join(weighted_unifrac_PC2, by = "sample-id") %>%
  left_join(weighted_unifrac_PC3, by = "sample-id") %>%
  left_join(weighted_unifrac_PC4, by = "sample-id") %>%
  left_join(weighted_unifrac_PC5, by = "sample-id")

# For canberra
combined_data <- combined_data %>%
  left_join(canberra_PC1, by = "sample-id") %>%
  left_join(canberra_PC2, by = "sample-id") %>%
  left_join(canberra_PC3, by = "sample-id") %>%
  left_join(canberra_PC4, by = "sample-id") %>%
  left_join(canberra_PC5, by = "sample-id")

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

# Mediation model for unweighted_unifrac_PC1
mediator_model <- lm(unweighted_unifrac_PC1 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ unweighted_unifrac_PC1 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_unweighted_unifrac_PC1 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "unweighted_unifrac_PC1", boot = TRUE, sims = 1000)
unweighted_unifrac_mediation_result_PC1 <- check_mediation(mediation_unweighted_unifrac_PC1)

# Mediation model for unweighted_unifrac_PC2
mediator_model <- lm(unweighted_unifrac_PC2 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ unweighted_unifrac_PC2 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_unweighted_unifrac_PC2 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "unweighted_unifrac_PC2", boot = TRUE, sims = 1000)
unweighted_unifrac_mediation_result_PC2 <- check_mediation(mediation_unweighted_unifrac_PC2)

# Mediation model for unweighted_unifrac_PC3
mediator_model <- lm(unweighted_unifrac_PC3 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ unweighted_unifrac_PC3 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_unweighted_unifrac_PC3 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "unweighted_unifrac_PC3", boot = TRUE, sims = 1000)
unweighted_unifrac_mediation_result_PC3 <- check_mediation(mediation_unweighted_unifrac_PC3)

# Mediation model for unweighted_unifrac_PC4
mediator_model <- lm(unweighted_unifrac_PC4 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ unweighted_unifrac_PC4 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_unweighted_unifrac_PC4 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "unweighted_unifrac_PC4", boot = TRUE, sims = 1000)
unweighted_unifrac_mediation_result_PC4 <- check_mediation(mediation_unweighted_unifrac_PC4)

# Mediation model for unweighted_unifrac_PC5
mediator_model <- lm(unweighted_unifrac_PC5 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ unweighted_unifrac_PC5 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_unweighted_unifrac_PC5 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "unweighted_unifrac_PC5", boot = TRUE, sims = 1000)
unweighted_unifrac_mediation_result_PC5 <- check_mediation(mediation_unweighted_unifrac_PC5)


# Print mediation results for beta diversity
print(paste("Unweighted UniFrac PC1 mediation result:", unweighted_unifrac_mediation_result_PC1))
print(paste("Unweighted UniFrac PC2 mediation result:", unweighted_unifrac_mediation_result_PC2))
print(paste("Unweighted UniFrac PC3 mediation result:", unweighted_unifrac_mediation_result_PC3))
print(paste("Unweighted UniFrac PC4 mediation result:", unweighted_unifrac_mediation_result_PC4))
print(paste("Unweighted UniFrac PC5 mediation result:", unweighted_unifrac_mediation_result_PC5))
print(" ")


# Mediation model for weighted_unifrac_PC1
mediator_model <- lm(weighted_unifrac_PC1 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ weighted_unifrac_PC1 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_weighted_unifrac_PC1 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "weighted_unifrac_PC1", boot = TRUE, sims = 1000)
weighted_unifrac_mediation_result_PC1 <- check_mediation(mediation_weighted_unifrac_PC1)

# Mediation model for weighted_unifrac_PC2
mediator_model <- lm(weighted_unifrac_PC2 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ weighted_unifrac_PC2 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_weighted_unifrac_PC2 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "weighted_unifrac_PC2", boot = TRUE, sims = 1000)
weighted_unifrac_mediation_result_PC2 <- check_mediation(mediation_weighted_unifrac_PC2)

# Mediation model for weighted_unifrac_PC3
mediator_model <- lm(weighted_unifrac_PC3 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ weighted_unifrac_PC3 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_weighted_unifrac_PC3 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "weighted_unifrac_PC3", boot = TRUE, sims = 1000)
weighted_unifrac_mediation_result_PC3 <- check_mediation(mediation_weighted_unifrac_PC3)

# Mediation model for weighted_unifrac_PC4
mediator_model <- lm(weighted_unifrac_PC4 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ weighted_unifrac_PC4 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_weighted_unifrac_PC4 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "weighted_unifrac_PC4", boot = TRUE, sims = 1000)
weighted_unifrac_mediation_result_PC4 <- check_mediation(mediation_weighted_unifrac_PC4)

# Mediation model for weighted_unifrac_PC5
mediator_model <- lm(weighted_unifrac_PC5 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ weighted_unifrac_PC5 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_weighted_unifrac_PC5 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "weighted_unifrac_PC5", boot = TRUE, sims = 1000)
weighted_unifrac_mediation_result_PC5 <- check_mediation(mediation_weighted_unifrac_PC5)

# Mediation model for canberra_PC1
mediator_model <- lm(canberra_PC1 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ canberra_PC1 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_canberra_PC1 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "canberra_PC1", boot = TRUE, sims = 1000)
canberra_mediation_result_PC1 <- check_mediation(mediation_canberra_PC1)

# Mediation model for canberra_PC2
mediator_model <- lm(canberra_PC2 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ canberra_PC2 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_canberra_PC2 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "canberra_PC2", boot = TRUE, sims = 1000)
canberra_mediation_result_PC2 <- check_mediation(mediation_canberra_PC2)

# Mediation model for canberra_PC3
mediator_model <- lm(canberra_PC3 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ canberra_PC3 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_canberra_PC3 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "canberra_PC3", boot = TRUE, sims = 1000)
canberra_mediation_result_PC3 <- check_mediation(mediation_canberra_PC3)

# Mediation model for canberra_PC4
mediator_model <- lm(canberra_PC4 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ canberra_PC4 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_canberra_PC4 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "canberra_PC4", boot = TRUE, sims = 1000)
canberra_mediation_result_PC4 <- check_mediation(mediation_canberra_PC4)

# Mediation model for canberra_PC5
mediator_model <- lm(canberra_PC5 ~ alcohol_frequent, data = combined_data)
outcome_model <- glm(parkinson_disease ~ canberra_PC5 + alcohol_frequent, data = combined_data, family = binomial(link = "logit"))
mediation_canberra_PC5 <- mediate(mediator_model, outcome_model, treat = "alcohol_frequent", mediator = "canberra_PC5", boot = TRUE, sims = 1000)
canberra_mediation_result_PC5 <- check_mediation(mediation_canberra_PC5)

# Print mediation results for beta diversity PCs 1 to 5
print(paste("Weighted UniFrac PC1 mediation result:", weighted_unifrac_mediation_result_PC1))
print(paste("Weighted UniFrac PC2 mediation result:", weighted_unifrac_mediation_result_PC2))
print(paste("Weighted UniFrac PC3 mediation result:", weighted_unifrac_mediation_result_PC3))
print(paste("Weighted UniFrac PC4 mediation result:", weighted_unifrac_mediation_result_PC4))
print(paste("Weighted UniFrac PC5 mediation result:", weighted_unifrac_mediation_result_PC5))
print(" ")
print(paste("Canberra PC1 mediation result:", canberra_mediation_result_PC1))
print(paste("Canberra PC2 mediation result:", canberra_mediation_result_PC2))
print(paste("Canberra PC3 mediation result:", canberra_mediation_result_PC3))
print(paste("Canberra PC4 mediation result:", canberra_mediation_result_PC4))
print(paste("Canberra PC5 mediation result:", canberra_mediation_result_PC5))

      