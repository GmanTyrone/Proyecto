library(dplyr)
library(table1)
library(Hmisc)
library(Exact)
library(DescTools)
library(rlang)

###整理資料變數
metadata <- read.table("phyloseq/weis_metadata_1210.tsv", sep = "\t", header = TRUE)

#----------------------------------------------------------------------------
##### Specific Variables
# 論文裡面放的table1 output

# colnames(metadata_new)

data_table <- metadata %>%
  select(-ID) %>%
  upData(
    labels = c(
      Calprotectin_greater_than_50 = "Calprotectin greater than 50 mcg/g",
      Gastrointestinal_symptoms = "Gastrointestinal symptoms",
      Disease_duration = "Disease duration (months)",
      Hoehn_Yahr_stage = "Hoehn-Yahr stage",
      Average_L_dopa_dose = "Average L-dopa dose (last 2 years)",
      Family_history_for_neurodegenerative_disorders = "Family history for neurodegenerative disorders"
    )
  )

####
compared_group <- "Entacapone"
compared_levels <- c("yes", "no")
compared_labels <- c("Yes", "No")
####

data_table1 <- data_table
data_table1[[gsub("`", "", compared_group)]] <- factor(data_table1 %>% pull(!!sym(compared_group)), levels = c(compared_levels, Inf), labels = c(compared_labels, "P-value"))

data_table1 <- data_table1[!is.na(data_table1 %>% pull(!!sym(compared_group))), ]

pvalue <- function(x, ...) {
  y <- unlist(x)
  g <- factor(rep(1:length(x), times = sapply(x, length)))
  if (is.numeric(y)) {
    p <- wilcox.test(y ~ g, exact = FALSE)$p.value
  } else if (table(y, g) %>% dim %>% identical(c(2, 2))) {
    p <- exact.test(table(y, g), method = "CSM")$p.value
  } else {
    p <- GTest(table(y, g), correct = "williams")$p.value
  }
  c("", sub("<", "&lt;", format.pval(p, digits = 3, eps = 0.001)))
}

my.render.cont <- function(x) {
  with(
    stats.apply.rounding(stats.default(x), digits = 3), 
    c("", "Mean (SD)" = sprintf("%s (%s)", MEAN, SD))
  )
}
##  table 1

(tb1 <- colnames(data_table1) %>%
    .[. != gsub("`", "", compared_group)] %>%
    paste0("`", ., "`", collapse = " + ") %>%
    paste0("~ ", ., " | ", compared_group) %>%
    as.formula %>%
    table1(data = data_table1, 
           overall = F, 
           render.continuous = my.render.cont,
           extra.col = list(`P-value` = pvalue)))
