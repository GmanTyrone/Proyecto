library(xml2)
library(dplyr)
library(XML)
library(stringr)
library(rvest)
library(magrittr)


report <- read.table("filereport_read_run_PRJEB30615_tsv.txt", sep = "\t", header = TRUE)

metadata0 <- list()
for (i in 1:nrow(report)){
  url <- paste("https://www.ebi.ac.uk/ena/browser/api/xml/", report$sample_accession[[i]], sep = "")
  aqi <- read_xml(url)
  primary_id <- aqi %>%
    xml_find_all(xpath = "//SAMPLE_SET/SAMPLE/IDENTIFIERS/PRIMARY_ID") %>%
    xml_text()
  external_id <- aqi %>%
    xml_find_all(xpath = "//SAMPLE_SET/SAMPLE/IDENTIFIERS/EXTERNAL_ID") %>%
    xml_text()
  submitter_id <- aqi %>%
    xml_find_all(xpath = "//SAMPLE_SET/SAMPLE/IDENTIFIERS/SUBMITTER_ID") %>%
    xml_text()
  title <- aqi %>%
    xml_find_all(xpath = "//SAMPLE_SET/SAMPLE/TITLE") %>%
    xml_text()
  sample <- aqi %>%
    xml_find_all(xpath = "//SAMPLE_SET/SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE") %>%
    xml_text() %>%
    str_replace_all(c("\\(" = "_", "\\)" = "_"))
  names <- aqi %>%
    xml_find_all(xpath = "//SAMPLE_SET/SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE/TAG") %>%
    xml_text() %>%
    str_replace_all(c("\\(" = "_", "\\)" = "_"))
  tmp <- matrix(c(report$run_accession[[i]],
                  report$sample_accession[[i]],
                  primary_id,
                  external_id,
                  submitter_id,
                  title,
                  sample), nrow = 1)
  tmp %<>% as.data.frame()
  colnames(tmp) <- c("ID", 
                     "sample_accession",
                     "primary_id",
                     "external_id",
                     "submitter_id",
                     "title",
                     names)
  metadata0[[i]] <- tmp
}
metadata <- bind_rows(metadata0)

sup_data <- readxl::read_xlsx("supplementary data.xlsx")
colnames(sup_data) %<>% str_replace_all(c("\\(" = "_", "\\)" = "_"))
sup_data %<>% apply(2, \(x) ifelse(x == "-", NA, x))
sup_data %<>% as.data.frame()

metadata$submitter_id %<>% str_split("_Illumina", simplify = TRUE) %>% .[,1]
for(i in 7:ncol(metadata)) {
  metadata[,i] <- metadata[,i] %>% str_split(colnames(metadata)[i], simplify = TRUE) %>% .[[2]] %>% {ifelse(`==`(.,""), NA_character_, .)}
}

metadata %<>% right_join(sup_data, by = "submitter_id")
metadata %<>% apply(2, \(x) ifelse(x == "Missing: Not provided", NA_character_, x)) %>%
  apply(2, \(x) ifelse(x == "not applicable", NA_character_, x)) %>%
  apply(2, \(x) ifelse(x == "NA", NA_character_, x)) %>%
  apply(2, \(x) ifelse(x == "-", NA_character_, x)) %>%
  apply(2, \(x) ifelse(x == "", NA_character_, x)) %>%
  as.data.frame()

metadata2 <- metadata %>% 
  mutate(ID = submitter_id,
         Sex = ifelse(Sex == "f", "Female",
               ifelse(Sex == "m", "Male", NA)),
         Phenotype = ifelse(Phenotype == "HR", "hypokinetic rigid", 
                     ifelse(Phenotype == "T", "tremor_dominant",
                     ifelse(Phenotype == "E", "equivalent", NA)))) %>% 
  rename("Gastrointestinal_symptoms" = "Gastrointestinal symptoms",
         "Disease_duration" = "Disease duration _months_",
         "Hoehn_Yahr_stage" = "Hoehn-Yahr stage",
         "Average_L_dopa_dose" = "Average L-dopa dose _last 2 years_",
         "Family_history_for_neurodegenerative_disorders" = "Family history for neurodegenerative disorders") %>%
  filter(Group == "PD",
         as.numeric(Average_L_dopa_dose) > 0) %>%
  select(ID,
         Entacapone,
         Sex,
         Age,
         Smoker,
         Calprotectin_greater_than_50,
         Constipation,
         Appendectomy,
         Gastrointestinal_symptoms,
         Disease_duration,
         Hoehn_Yahr_stage,
         Phenotype,
         Average_L_dopa_dose,
         Family_history_for_neurodegenerative_disorders) %>%
  mutate(across(c(Sex,
                  Smoker,
                  Calprotectin_greater_than_50,
                  Constipation,
                  Appendectomy,
                  Gastrointestinal_symptoms,
                  Hoehn_Yahr_stage,
                  Phenotype,
                  Average_L_dopa_dose,
                  Family_history_for_neurodegenerative_disorders), as.factor),
         across(c(Age, Disease_duration), as.numeric)) %>%
  distinct()

  
write.table(metadata %>% select(ID, submitter_id) %>% .[.$submitter_id %in% metadata2$ID,], "id_list.csv", sep = ",", col.names = FALSE, row.names = FALSE)
write.csv(metadata2, "weis_metadata_1210.csv", row.names = FALSE)
