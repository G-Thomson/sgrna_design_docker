suppressPackageStartupMessages(library(tidyverse))

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  action <- args[1]
  filename <- args[-1]
  
  if(action != "--in"){
    stop("Format_sgRNA_table.R needs an input", call. = FALSE)
  }
  if(length(filename) != 1){
    stop("Format_sgRNA_table.R can only take 1 input",call. = FALSE)
  }
  fexists = file.exists(filename)
  if(fexists != TRUE){
    stop("Format_sgRNA_table.R can't find that file")
  }
  
  process(filename)
}

process <- function(filename){
  
  output <- readr::read_csv(filename,
    col_names = c("sgRNA", "sgRNA_score", "result"),
    skip = 1,
    col_types = cols()) 
  
  output %>% 
    dplyr::mutate(sgRNA_ind = 1:n(), .before = sgRNA) %>% 
    splitstackshape::cSplit('result', ')') %>% 
    tidyr::pivot_longer(starts_with("result"), names_to = "gene_ind", values_to = "result") %>% 
    tidyr::separate(result, c("gene", "result"), " ", extra = "merge") %>% 
    dplyr::mutate(result = str_remove(result, '\\(')) %>% 
    splitstackshape::cSplit('result', ';') %>% 
    dplyr::rename(gene_score = result_1, pos_match1 = result_3, pos_match2 = result_5) %>%
    dplyr::mutate(result_2 = as.character(result_2),
                  result_4 = as.character(result_4)) %>% 
    dplyr::mutate(mm_match1 = case_when(nchar(result_2) == 0 ~ 0,
                                        TRUE ~ str_count(result_2, " ") + 1),
                  mm_match2 = case_when(nchar(result_4) == 0 ~ 0,
                                        TRUE ~ str_count(result_4, " ") + 1),
                  pos_match1 = str_replace(pos_match1, "pos: ", ""),
                  pos_match2 = str_replace(pos_match2, "pos: ", "")) %>% 
    dplyr::select(-result_2, -result_4) %>% 
    dplyr::mutate_at(vars(matches("match")), as.character) %>% 
    tidyr::pivot_longer(cols = contains("match"),
                        names_to = "variable",
                        values_to = "dat") %>% 
    tidyr::separate(variable, c("variable", "match_ind"), -1) %>% 
    tidyr::pivot_wider(names_from = variable, values_from = dat) %>% 
    dplyr::filter(!is.na(pos_match),
                  mm_match == 0) %>% 
    dplyr::mutate(strand = if_else(str_detect(pos_match, "R"), "reverse", "forward"), 
                  pos_match = str_replace(pos_match, "R", "")) %>% 
    dplyr::group_by(sgRNA) %>% 
    dplyr::mutate(min_pos = min(pos_match)) %>% 
    readr::write_csv("Formated_sgRNA_table.csv")
}

main()