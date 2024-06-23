# if(!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R")
library(ggplot2)
library(magrittr)
library(qiime2R)
metadata <- read.table("phyloseq/weis_metadata_1210.tsv", sep = "\t", header = TRUE)

# Canberra
canberra_poca <- "metrics-results/canberra_pcoa_matrix.qza"
canberra_poca <- qiime2R::read_qza(canberra_poca) %>% .$data
p_canberra <- ggplot(cbind(canberra_poca$Vectors, metadata), aes(x = PC1, y = PC2, color = Entacapone)) +
  geom_point() +
  stat_ellipse(type = "t", geom="point", size=0.5, aes(group = Entacapone), level=0.95, linetype=2, colour="black") +
  scale_colour_manual("Entacapone", values=c("#E69F00", "#0072B2")) +
  labs(title = "Canberra PCoA",
       x = paste0("Axis.1 [", (canberra_poca$ProportionExplained[[1]] * 100) %>% round(1), "%]"),
       y = paste0("Axis.2 [", (canberra_poca$ProportionExplained[[2]] * 100) %>% round(1), "%]")) +
  theme_classic()

print(p_canberra)


# Unweighted Unifrac
unweighted_unifrac_poca <- "metrics-results/unweighted_unifrac_pcoa_matrix.qza"
unweighted_unifrac_poca <- qiime2R::read_qza(unweighted_unifrac_poca) %>% .$data
p_unweighted_unifrac <- ggplot(cbind(unweighted_unifrac_poca$Vectors, metadata), aes(x = PC1, y = PC2, color = Entacapone)) +
  geom_point() +
  stat_ellipse(type = "t") +
  scale_colour_manual("Entacapone", values=c("#E69F00", "#0072B2")) +
  labs(title = "Unweighted Unifrac PCoA",
       x = paste0("Axis.1 [", (unweighted_unifrac_poca$ProportionExplained[[1]] * 100) %>% round(1), "%]"),
       y = paste0("Axis.2 [", (unweighted_unifrac_poca$ProportionExplained[[2]] * 100) %>% round(1), "%]")) +
  theme_classic()

print(p_unweighted_unifrac)

# Weighted Unifrac
weighted_unifrac_poca <- "metrics-results/weighted_unifrac_pcoa_matrix.qza"
weighted_unifrac_poca <- qiime2R::read_qza(weighted_unifrac_poca) %>% .$data
p_weighted_unifrac <- ggplot(cbind(weighted_unifrac_poca$Vectors, metadata), aes(x = PC1, y = PC2, color = Entacapone)) +
  geom_point() +
  stat_ellipse(type = "t") +
  scale_colour_manual("Entacapone", values=c("#E69F00", "#0072B2")) +
  labs(title = "Weighted Unifrac PCoA",
       x = paste0("Axis.1 [", (weighted_unifrac_poca$ProportionExplained[[1]] * 100) %>% round(1), "%]"),
       y = paste0("Axis.2 [", (weighted_unifrac_poca$ProportionExplained[[2]] * 100) %>% round(1), "%]")) +
  theme_classic()

print(p_weighted_unifrac)
