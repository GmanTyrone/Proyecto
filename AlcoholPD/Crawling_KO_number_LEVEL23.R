library(rvest)
library(magrittr)
library(stringr)
library(purrr)
library(readr)
library(dplyr)


url <- "https://www.genome.jp/kegg/pathway.html"

sections <- url %>%
  read_html() %>%
  html_nodes(".list") %>%
  .[1:13]

table <- lapply(1:13,
                \(i) data.frame(
                  number = paste0("ko",
                                  html_nodes(sections[i], "dt") %>%
                                    html_text %>%
                                    str_sub(1,5)),
                  name = paste0(LETTERS[i], 
                                "_",
                                html_nodes(sections[i], "dd") %>% html_text) %>%
                    str_split("Including:", simplify = TRUE) %>% .[,1] %>%
                    str_replace_all(c("/" = "_"))
                )) %>% bind_rows()
write_tsv(table, file = "KO_LEVEL23_2022_1214.tsv", col_names = FALSE)
