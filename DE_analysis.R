rm(list = ls())
library(tidyverse)
library(sleuth)
(WD <- getwd())
if (!is.null(WD)) setwd(WD)

myfun <- function(i,j,list,df_info){
  
  df_transcript_info <- df_info %>%
    pull(path) %>%
    head(1) %>%
    file.path(., "abundance.tsv") %>%
    read_tsv() %>%
    separate(
      target_id, 
      c("ens_transcript", "ens_gene", NA, NA, NA, "ext_gene", NA, "role", 
        NA), 
      sep = "\\|", remove = FALSE
    ) %>%
    select(target_id, ens_transcript, ens_gene, ext_gene, role) 
  
  so <- sleuth_prep(
    df_info,
    target_mapping = df_transcript_info, 
    extra_bootstrap_summary = TRUE)
  
  so <- sleuth_fit(so, ~condition, "full")
  so <- sleuth_fit(so, ~1, "reduced")
  so <- sleuth_lrt(so, "reduced", "full")
  models(so)
  tests(so)
  
  sleuth_table <- sleuth_results(
    so,
    "reduced:full",
    "lrt", 
    # "reduced",
    # "lrt", 
    show_all = FALSE
  ) %>%
    filter(qval <= 0.05)
  
  if (nrow(sleuth_table) == 0){
    # do nothing
  }
  else{
    sleuth_results(
      so,
      "reduced:full",
      "lrt", 
      # "reduced",
      # "lrt", 
      show_all = FALSE
    )
    
    jpeg(sprintf("%s_%s.jpg",list[i],list[j]), width = 1000, height = 500)
    plt <- cowplot::plot_grid(
      plot_bootstrap(so, sleuth_table %>% pull(target_id) %>% head(1)) + ggtitle("top DE"),
      plot_bootstrap(so, sleuth_table %>% pull(target_id) %>% tail(1)) + ggtitle("bottom DE")
    )
    print(plt)
    dev.off()
    
    
  }
  # browser()
  save(list = ls(all.names = TRUE), file = sprintf("B:/OneDrive - ETH Zurich/ETHZ/2021 spring semester/Lab/Mingzhe/results/DE/%s_%s.RData",list[i],list[j]))
  return(nrow(sleuth_table))
}


  # kallisto_quant_root = "C:/Users/Mingzhe/Desktop/Rworkspace/quant"
  # kallisto_quant_root = "B:/OneDrive - ETH Zurich/ETHZ/2021 spring semester/Lab/NGS/backup"
  kallisto_quant_root = "B:/OneDrive - ETH Zurich/ETHZ/2021 spring semester/Lab/NGS/Kallisto"
  list.files(kallisto_quant_root, full.names = TRUE) %>%
    map_dfr(function(path) {
      sample <- basename(path)
      # Here change the index according to the way you slicing the name of your folder
      condition <- strsplit(sample, "_")[[1]][[2]]
      # condition <- strsplit(condition_temp, "\\.")[[1]][[1]]
      data.frame(sample = sample, condition = condition, path = path)
    }) -> df_info_all
 
  list = c("AUntreated","LPS","C5a","LPS+C5a")
  list[1]
  list[2]
  
  table <- c()
  for (i in 1:3){
          for (j in (i+1):4){
  #for (i in 1){
    #for (j in 2){
      df_info_all %>% filter(condition %in% c(list[i],list[j])) -> df_info
      number <- myfun(i,j,list,df_info)
      table <- c(table,number)
    }
  }
  
  save.image("Total_local_k_.RData")