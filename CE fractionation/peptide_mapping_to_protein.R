library(seqinr)
library(stringr)
library(readr)
library(dplyr)
library(ggplot2)

peptide_df_gen <- function(protein_fasta, peptide_file, group = 'default'){
  
  protein_seq_list <- read.fasta(protein_fasta, seqtype = 'AA', seqonly = TRUE)

  protein_seq <- protein_seq_list[[1]]
  
  peptides <- read_csv(peptide_file, 
                       col_names = FALSE)
  
  return_df <- data.frame(peptide_seq=character(),
                          position=integer())
  
  for (i in 1:nrow(peptides)){
    temp_peptide_df <- single_peptide_df(peptides$X1[i], protein_seq)
    return_df <- bind_rows(return_df, temp_peptide_df)
  }
  
  return_df$group <- group
  
  return(return_df)
}


#generates all the indices of a single peptide
single_peptide_df <- function(peptide, protein_seq){
  peptide_location <- str_locate_all(protein_seq, peptide)
  matches <- str_count(protein_seq, peptide)

  peptide_pos_df <- data.frame(peptide_seq=character(),
                             position=integer())
  for (j in 1:matches){
  
    position <- rep(NA, nchar(peptide))
  
    pos <- peptide_location[[1]][j]
  
    for (i in 1:nchar(peptide)){
      position[i] <- pos
      pos <- pos+1
    }
    peptide_seq <- rep(peptide, nchar(peptide))
    temp_pep_df <- data.frame(peptide_seq,position)
    #print(temp_pep_df)
    peptide_pos_df <- bind_rows(peptide_pos_df, temp_pep_df)
  }  
  return(peptide_pos_df)
}

remove_mods_str <- function(string){
  return_string <- str_replace_all(string, "[a-z1-9()+-:.]", "")
  return(return_string)
}

remove_mods_unique_file <- function(text_file){
  strings <- read_csv(text_file, 
                       col_names = FALSE)
  
  basefile <- tools::file_path_sans_ext(text_file)
  
  strings$X1 <- remove_mods_str(strings$X1)
  
  return_strings <- unique(strings)
  
  write.table(return_strings, paste0(basefile, "_unique.txt"),col.names = FALSE, row.names = FALSE)
}


