# library(phyloseq)
library(tidyverse)
library(magrittr)
library(openxlsx)
library(MicrobiotaProcess)
library(treeio)
library(ape)
library(stringr)
library(coin)

### Wilcoxon in asv
## prepare data
otu <- "phyloseq/rel-feature-table.qza"
tax <- "phyloseq/taxonomy.qza"
metadata <- read.table("phyloseq/weis_metadata_1210.tsv", sep = "\t", header = TRUE)

asv <- otu %>% MicrobiotaProcess::read_qza(.) %>% .$otutab %>% data.frame
tax %<>% MicrobiotaProcess::read_qza(.) %>% .$taxatab  %>% data.frame
Genus <-  str_split(tax$Genus, "g__", simplify = TRUE)[,2]
Genus[Genus  == ""] <- "Unidentified"
Phylum <- str_split(tax$Phylum, "p__", simplify = TRUE)[,2]
Phylum[Phylum  == ""] <- "Unidentified"

Genus %<>% str_replace_all(c("-" = "_", "\\[" = "", "\\]" = ""))
Phylum %<>% str_replace_all(c("-" = "_", "\\[" = "", "\\]" = ""))

rownames(asv) <- paste0("asv", 1:nrow(asv))
Group <- ifelse(metadata$Entacapone == "yes", 0L, 1L)
few <- apply(asv, 1, function(x) mean(x > 0) < 0.1)
Genus <- Genus[!few]
Phylum <- Phylum[!few]
asv_new <- asv[!few,]
asv_new <- asv_new[,colnames(asv_new) %in% metadata$ID]


new_data <- cbind(Group, t(asv_new)) %>% as.data.frame
new_data <- new_data[!is.na(new_data$Group),]
new_data$Group <- as.factor(new_data$Group)
k <- aggregate(new_data[2:ncol(new_data)], list(new_data$Group), mean)
tmp <- apply(k, 2, function(x) all(x==0))
k <- k[, !tmp]
new_data <- new_data[, !tmp]
Genus <- Genus[!tmp[-1]]
Phylum <- Phylum[!tmp[-1]]
k <- k[,-1]

result <- data.frame()
for (i in 1:ncol(k)){
  FoldChange <- (k[2,i] + .Machine$double.eps)/(k[1,i] + .Machine$double.eps)
  result[i,1] <- FoldChange
}

## wilcoxon test
for (i in 1:ncol(k)){
  formula <- reformulate("Group", response =  colnames(new_data)[i+1])
  p <- coin::wilcox_test(formula,data=new_data) %>% pvalue
  # equal_var <- bartlett.test(formula, data=new_data)$p.value
  # if(equal_var > 0.05){
  #   p <- t.test(formula, data=new_data, var.equal=TRUE)$p.value
  # }else{
  #   p <- t.test(formula, data=new_data, var.equal=FALSE)$p.value
  # }
  result[i,2] <- p
}

for (i in 1:ncol(k)){
  criteria <- ifelse(k[1,i] == 0 | k[2,i] == 0, "1", "0")
  result[i,3] <- criteria
}

result$Phylum <- Phylum
result$Genus <- Genus
colnames(result) <- c("FoldChange","pvalue","Criteria","Phylum","Genus")
result <- result[order(result$pvalue),]


## volcano plot
# colored by phylum
result$sig_p <- result$Phylum %>%
  ifelse(is.na(.), "Unidentified", .) %>%
  ifelse(result$pvalue < 0.05, ., NA)
result$sig_p %<>% factor(levels = {
  result$sig_p %>% .[. != "Unidentified"] %>% factor %>% levels %>% c(., "Unidentified")
})

print(
result %>%
  ggplot() +
  geom_point(aes(x = log2(FoldChange), y = -log10(pvalue)), color = "gray", size = 1) +
  geom_point(data = result[result$pvalue < 0.05,], aes(x = log2(FoldChange), y = -log10(pvalue), color = sig_p), size = 1) +
  labs(x = "log2 fold change", y = expression(-log10*paste(" ","p-value")), color="Phylum") +
  scale_y_continuous(breaks = seq(from = 0, to = 3, by = 0.5),limits = c(0,3)) +
  scale_x_continuous(breaks = seq(from = -45, to = 45, by = 15),limits = c(-50,50)) +
  geom_hline(yintercept = -log10(0.05),color='red',linetype = "dashed")+
  geom_vline(xintercept = max(log2(result[result$Criteria==0,1])),color='red')+
  geom_vline(xintercept = min(log2(result[result$Criteria==0,1])),color='red')+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "right",
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5)),
        legend.text = element_text(size=12),
        legend.title  = element_text(size=13),
        legend.margin = margin(17, 2, -25, 2),
        plot.margin = margin(0.6, 0.6 , 0.6, 0.3, "cm"),
        axis.title.x = element_text(vjust = -0.8),
        axis.title.y = element_text(vjust = 1.3))
)

# colored by genus
result$sig_g <- result$Genus %>%
  ifelse(is.na(.), "Unidentified", .) %>%
  ifelse(result$pvalue < 0.05, ., NA)
result$sig_g %<>% factor(levels = {
  result$sig_g %>% .[. != "Unidentified"] %>% factor %>% levels %>% c(., "Unidentified")
})
print(
  result %>%
    ggplot() +
    geom_point(aes(x = log2(FoldChange), y = -log10(pvalue)), color = "gray", size = 1) +
    geom_point(data = result[result$pvalue < 0.05,], aes(x = log2(FoldChange), y = -log10(pvalue), color = sig_g), size = 1) +
    labs(x = "log2 fold change", y = expression(-log10*paste(" ","p-value")), color="Genus") +
    scale_y_continuous(breaks = seq(from = 0, to = 3, by = 0.5),limits = c(0,3)) +
    scale_x_continuous(breaks = seq(from = -45, to = 45, by = 15),limits = c(-50,50)) +
    geom_hline(yintercept = -log10(0.05),color='red',linetype = "dashed")+
    geom_vline(xintercept = max(log2(result[result$Criteria==0,1])),color='red')+
    geom_vline(xintercept = min(log2(result[result$Criteria==0,1])),color='red')+
    theme(panel.background = element_blank(),
          panel.border = element_rect(colour = "black",fill=NA),
          legend.position = "bottom",  # Place legend at the bottom
          axis.title = element_text(size = rel(2)),
          axis.text = element_text(size = rel(1.5)),
          legend.text = element_text(size=12),
          legend.title  = element_text(size=13),
          legend.margin = margin(17, 2, -25, 2),
          plot.margin = margin(0.6, 0.6 , 1.5, 0.3, "cm"),  # Adjust bottom margin here
          axis.title.x = element_text(vjust = -0.8),
          axis.title.y = element_text(vjust = 1.3))
  
)

### Wilcoxon in Genus level
otu <- "phyloseq/rel-genus-table.qza"
metadata <- read.table("phyloseq/weis_metadata_1210.tsv", sep = "\t", header = TRUE)

asv <- otu %>% MicrobiotaProcess::read_qza(.) %>% .$otutab %>% data.frame
Genus <- str_split(rownames(asv), ";g__", simplify = TRUE)[,2] %>% paste0("g__", .)
Genus[Genus  == "g__"] <- paste0("Unidentified_", 1:sum(Genus == "g__"))
for(name in Genus[duplicated(Genus)] %>% unique) Genus[Genus == name] <- paste0(name, 1:sum(Genus == name))
Genus %<>% str_replace_all(c("-" = "_", "\\[" = "", "\\]" = ""))
rownames(asv) <- Genus
Group <- ifelse(metadata$Entacapone == "yes", 0L, 1L)
few <- apply(asv, 1, function(x) mean(x > 0) < 0.1)
Genus <- Genus[!few]
asv_new <- asv[!few,]
asv_new <- asv_new[,colnames(asv_new) %in% metadata$ID]


new_data <- cbind(Group, t(asv_new)) %>% as.data.frame
new_data <- new_data[!is.na(new_data$Group),]
new_data$Group <- as.factor(new_data$Group)
k <- aggregate(new_data[2:ncol(new_data)], list(new_data$Group), mean)
k <- k[, !(apply(k, 2, function(x) all(x==0)))]
k <- k[,-1]


result <- data.frame()
for (i in 1:ncol(k)){
  FoldChange <- (k[2,i] + .Machine$double.eps)/(k[1,i] + .Machine$double.eps)
  result[i,1] <- FoldChange
}

# wilcoxon test
for (i in 1:ncol(k)){
  formula <- as.formula(paste(colnames(new_data)[i+1],"~","Group"))
  p <- coin::wilcox_test(formula,data=new_data) %>% pvalue
  result[i,2]<-p
}

for (i in 1:ncol(k)){
  criteria <- ifelse(k[1,i] == 0 | k[2,i] == 0, "1", "0")
  result[i,3] <- criteria
}

result$Genus <- Genus
colnames(result) <- c("FoldChange","pvalue","Criteria","Genus")
result <- result[order(result$pvalue),]

# save
# wb <- createWorkbook()
# addWorksheet(wb, "All")
# writeData(wb, "All", result %>% select(FoldChange, pvalue, Genus))
# saveWorkbook(wb, file = "Wilcoxon_Genus.xlsx", overwrite = TRUE)