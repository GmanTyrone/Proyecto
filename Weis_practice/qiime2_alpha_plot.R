# if(!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R")
library(tidyverse)
library(magrittr)
library(qiime2R)
library(ggh4x)
metadata <- read.table("phyloseq/weis_metadata_1210.tsv", sep = "\t", header = TRUE)
observed_features <- "metrics-results/observed_features_vector.qza"
observed_features <- qiime2R::read_qza(observed_features) %>% .$data
chao1 <- "metrics-results/chao1_vector.qza"
chao1 <- qiime2R::read_qza(chao1) %>% .$data
shannon <- "metrics-results/shannon_vector.qza"
shannon <- qiime2R::read_qza(shannon) %>% .$data
simpson <- "metrics-results/simpson_vector.qza"
simpson <- qiime2R::read_qza(simpson) %>% .$data
alphaobj <- cbind(observed_features, chao1, shannon, simpson)
colnames(alphaobj) <- c("Observed", "Chao1", "Shannon", "Simpson")
alphaobj <- cbind(metadata, alphaobj)

print(
  alphaobj %>%
    pivot_longer(c("Observed", "Chao1", "Shannon", "Simpson"), names_to = "metric", values_to = "Alpha Diversity Measure") %>%
    filter(!is.na(Entacapone)) %>%
    ggplot() +
    geom_boxplot(aes(y = `Alpha Diversity Measure`, fill = Entacapone), outlier.shape = 21, width = 0.5) +
    facet_grid(. ~ metric, scales = "free_y") +
    theme_minimal() +
    theme(axis.text=element_text(size=10), 
          axis.title=element_text(size=10)) +
    scale_x_discrete(labels = element_blank(), name = "") +
    scale_fill_manual("Entacapone", values=c("#E69F00", "#0072B2"))
)
