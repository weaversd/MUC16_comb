#import packages
suppressPackageStartupMessages({
  library(stringr)
  library(readxl)
  library(seqinr)
  library(dplyr)
  library(ggplot2)
})

setwd("~/03_Notre_Dame/01_Dissertation/Chapter-B Ovarian Cancer Method Development and Data analyisis/Code/nanopore")
####import sequences####
#repeats 8-12 for ovcar3 and 6sample consensus
novel_repeat_seq_ovcar3 <- "TAGPLLVPFTLNFTITNLQYEEDMHRPGSRKFNATERVLQGLLSPIFKNSSVGPLYSGCRLTSLRPEKDGAATGMDAVCLYHPNPKRPGLDREQLYWELSQLTHNITELGPYSLDRDSLYVNGFTHQNSVPTTSTPGTSTVYWATTGTPSSFPGHTEPGPLLIPFTFNFTITNLHYEENMQHPGSRKFNTTERVLQGLLTPLFKNTSVGPLYSGCRLTLLRPEKQEAATGVDTICTHRVDPIGPGLDRERLYWELSQLTNSITELGPYTLDRDSLYVNGFNPWSSVPTTSTPGTSTVHLATSGTPSSLPGHTAPVPLLIPFTLNFTITNLHYEENMQHPGSRKFNTTERVLQGLLKPLFKSTSVGPLYSGCRLTLLRPEKHGAATGVDAICTLRLDPTGPGLDRERLYWELSQLTNSVTELGPYTLDRDSLYVNGFTHRSSVPTTSIPGTSAVHLETSGTPASLPGHTAPGPLLVPFTLNFTITNLQYEEDMRHPGSRKFNTTERVLQGLLKPLFKSTSVGPLYSGCRLTLLRPEKRGAATGVDTICTHRLDPLNPGLDREQLYWELSKLTRGIIELGPYLLDRGSLYVNGFTHRNFVPITSTPGTSTVHLGTSETPSSLPRPIVPGPLLVPFTLNFTITNLQYEEAMRHPGSRKFNTTERVLQGLLRPLFKNTSIGPLYSSCRLTLLRPEKDKAATRVDAICTHHPDPQSPGLNREQLYWELSQLTHGITELGPYTLDRDSLYVDGFTHWSPIPTTSTPGTSIVNLGTSGIPPSLPETT"
novel_repeat_seq_6samp_concensus <- "TAGPLLVPFTLNFTITNLQYEEDMHRPGSRKFNTTERVLQGLLSPIFKNSSVGPLYSGCRLTSLRPEKDGAATGMDAVCLYHPNPKRPGLDREQLYWELSQLTHNITELGPYSLDRDSLYVNGFTHQNSVPTTSTPGTSTVYWATTGTPSSFPGHTEPGPLLIPFTFNFTITNLHYEENMQHPGSRKFNTTERVLQGLLTPLFKNTSVGPLYSGCRLTLLRPEKHEAATGVDTICTHRVDPIGPGLDRERLYWELSQLTNSITELGPYTLDRDSLYVNGFNPWSSVPTTSTPGTSTVHLATSGTPSSLPGHTAPVPLLIPFTLNFTITNLHYEENMQHPGSRKFNTTERVLQGLLKPLFKSTSVGPLYSGCRLTLLRPEKHGAATGVDAICTLRLDPTGPGLDRERLYWELSQLTNSVTELGPYTLDRDSLYVNGFTHRSSVPTTSIPGTSAVHLETSGTPASLPGHTAPGPLLVPFTLNFTITNLQYEEDMRHPGSRKFNTTERVLQGLLKPLFKSTSVGPLYSGCRLTLLRPEKRGAATGVDTICTHRLDPLNPGLDREQLYWELSKLTRGIIELGPYLLDRGSLYVNGFTHRNFVPITSTPGTSTVHLGTSETPSSLPRPIVPGPLLVPFTLNFTITNLQYEEAMRHPGSRKFNTTERVLQGLLRPLFKNTSIGPLYSSCRLTLLRPEKDKAATRVDAICTHHPDPQSPGLNREQLYWELSQLTHGITELGPYTLDRDSLYVDGFTHWSPIPTTSTPGTSIVNLGTSGIPPSLPETT"
#import entire MUC16 sequence, first as dataframe and then convert to string
muc16_file_ovcar3 <- read.csv("MUC16_OVCAR3_v1.txt", header = F)
muc16_seq_ovcar3 <- muc16_file_ovcar3[1,1]
muc16_file_6samp_concensus <- read.csv("MUC16_consensus_from 6_sequences.txt", header = F)
muc16_seq_6samp_concensus <- muc16_file_ovcar3[1,1]
#import proteome. This is the human proteome w/o MUC16 entry (every other protein)
proteome <- read.fasta(file = "Homo Sapiens UP000005640 no MUC16.fasta", seqtype = "AA", as.string = T)

####import peptide files####
####For each sample, create new column that is peptide sequence w/o mods.
####For each sample, create a dataframe with only MUC16 peptides.
Ovcar_ov3 <- read.csv("IP_nanopore_SDW.OVCAR3_db_search_ovcar3_and_fitz_IPs/db.OVCAR3_IP.peptides.csv")
Ovcar_ov3$sequence <- str_remove_all(Ovcar_ov3$Peptide, "[()\\.0-9\\-+]")
Ovcar_ov3_MUC16 <- Ovcar_ov3[str_detect(Ovcar_ov3$Accession, "MUC16") ,]
Ovcar_ov3_MUC16_peps <- Ovcar_ov3_MUC16$sequence
#Fitzgerald sample searched against 6sample consensus seq cont. database
Fitz_6s <- read.csv("IP_nanopore_SDW.6sampConcensus_db_Search_ovcar3_and_Fitz_IPs/db.Fitzgerald_IP.peptides.csv")
Fitz_6s$sequence <- str_remove_all(Fitz_6s$Peptide, "[()\\.0-9\\-+]")
Fitz_6s_MUC16 <- Fitz_6s[str_detect(Fitz_6s$Accession, "MUC16") ,]
Fitz_6s_MUC16_peps <- Fitz_6s_MUC16$sequence


####convert proteome to a string####
#create empty list, fill it with each entry in the proteome object
proteome_2 <- list()
for (i in 1:length(proteome)) {
  proteome_2[[i]] <- proteome[[i]][1]
}

#turn the list into a vector, and collapse it into one long string
proteome_vector <- unlist(proteome_2)
proteome_string <- paste(proteome_vector, sep = '', collapse = '')
proteome_string2 <- paste(proteome, sep = '', collapse = '')

####ANALYSIS FOR FITZGERALD####
Fitz_6s_MUC16$in_repeat <- F
Fitz_6s_MUC16$count <- NA
Fitz_6s_MUC16$count_proteome <- NA
Fitz_6s_MUC16$count_in_repeat <- NA
for (i in 1:nrow(Fitz_6s_MUC16)) {
  Fitz_6s_MUC16$in_repeat[i] <- grepl(Fitz_6s_MUC16$sequence[i], novel_repeat_seq_6samp_concensus)
  Fitz_6s_MUC16$count[i] <- str_count(muc16_seq_6samp_concensus, Fitz_6s_MUC16$sequence[i])
  Fitz_6s_MUC16$count_proteome[i] <- str_count(proteome_string2, Fitz_6s_MUC16$sequence[i])
  Fitz_6s_MUC16$count_in_repeat[i] <- str_count(novel_repeat_seq_6samp_concensus, Fitz_6s_MUC16$sequence[i])
}

Fitz_6s_MUC16$unique_to_repeat <- F
for (i in 1:nrow(Fitz_6s_MUC16)) {
    if (Fitz_6s_MUC16$in_repeat[i] == T && Fitz_6s_MUC16$count[i] == Fitz_6s_MUC16$count_in_repeat[i] && Fitz_6s_MUC16$count_proteome[i] == 0) {
      Fitz_6s_MUC16$unique_to_repeat[i] <- T
    }
}
peptides_of_interest_Fitz_6s <- Fitz_6s_MUC16[Fitz_6s_MUC16$unique_to_repeat == T,]
peptide_list_Fitz_6s <- peptides_of_interest_Fitz_6s$sequence
peptide_list_final_Fitz_6s <- unique(peptide_list_Fitz_6s)

#calculate percent coverage of repeats 8-12
#pull out peptides that mapped to R8-12
R8_12 <- Fitz_6s_MUC16[Fitz_6s_MUC16$in_repeat == T,]
AA_vec <- str_split(novel_repeat_seq_6samp_concensus, pattern = "")[[1]]
AA_df <- data.frame(AA = AA_vec,
                    covered = F)
for (i in 1:nrow(R8_12)) {
  matches <- data.frame(str_locate_all(novel_repeat_seq_6samp_concensus, R8_12$sequence[i]))
  for (j in 1:nrow(matches)) {
    start <- matches$start[j]
    end <- matches$end[j]
    sequence <- start:end
    for (k in sequence) {
      AA_df$covered[k] <- T
    }
  }
}

percent_coverage_8_12 <- (sum(AA_df$covered) / nrow(AA_df)) * 100
#print results
cat(paste0("Results for Fitzgerald Sample mapped to 6sample consensus:", "\n", 
           "MUC16 peptides: ", nrow(Fitz_6s_MUC16), "\n",
           "MUC16 sequences: ", length(unique(Fitz_6s_MUC16$sequence)), "\n",
           "R8 - R12 peptides: ", nrow(R8_12), "\n",
           "R8 - R12 sequences: ", length(unique(R8_12$sequence)), "\n",
           "R8 - R12 percent coverage: ", round(percent_coverage_8_12, 2),"%", "\n",
           "sequences unique to R8-12: ", length(peptide_list_final_Fitz_6s), "\n"))
peptide_list_final_Fitz_6s


#save peptides to .txt file
write.table(peptide_list_final_Fitz_6s, "unique_repeat_peptides_Fitz_IP_6s_db.txt",
            row.names = F, col.names = F, quote = F)

####ANALYSIS FOR OVCAR3####
#same analysis as above, but for the OVCAR3 IP sample
#ovcar3 db
Ovcar_ov3_MUC16$in_repeat <- F
Ovcar_ov3_MUC16$count <- NA
Ovcar_ov3_MUC16$count_proteome <- NA
Ovcar_ov3_MUC16$count_in_repeat <- NA
for (i in 1:nrow(Ovcar_ov3_MUC16)) {
  Ovcar_ov3_MUC16$in_repeat[i] <- grepl(Ovcar_ov3_MUC16$sequence[i], novel_repeat_seq_ovcar3)
  Ovcar_ov3_MUC16$count[i] <- str_count(muc16_seq_ovcar3, Ovcar_ov3_MUC16$sequence[i])
  Ovcar_ov3_MUC16$count_proteome[i] <- str_count(proteome_string2, Ovcar_ov3_MUC16$sequence[i])
  Ovcar_ov3_MUC16$count_in_repeat[i] <- str_count(novel_repeat_seq_ovcar3, Ovcar_ov3_MUC16$sequence[i])
}
Ovcar_ov3_MUC16$unique_to_repeat <- F
for (i in 1:nrow(Ovcar_ov3_MUC16)) {
  if (Ovcar_ov3_MUC16$in_repeat[i] == T && Ovcar_ov3_MUC16$count[i] == Ovcar_ov3_MUC16$count_in_repeat[i] && Ovcar_ov3_MUC16$count_proteome[i] == 0) {
    Ovcar_ov3_MUC16$unique_to_repeat[i] <- T
  } 
}
peptides_of_interest_Ovcar3_ov3 <- Ovcar_ov3_MUC16[Ovcar_ov3_MUC16$unique_to_repeat == T,]
peptide_list_Ovcar3_ov3 <- peptides_of_interest_Ovcar3_ov3$sequence
peptide_list_final_Ovcar3_ov3 <- unique(peptide_list_Ovcar3_ov3)

#calculate percent coverage of repeats 8-12
#pull out peptides that mapped to R8-12
R8_12 <- Ovcar_ov3_MUC16[Ovcar_ov3_MUC16$in_repeat == T,]
AA_vec <- str_split(novel_repeat_seq_ovcar3, pattern = "")[[1]]
AA_df <- data.frame(AA = AA_vec,
                    covered = F)
for (i in 1:nrow(R8_12)) {
  matches <- data.frame(str_locate_all(novel_repeat_seq_ovcar3, R8_12$sequence[i]))
  for (j in 1:nrow(matches)) {
    start <- matches$start[j]
    end <- matches$end[j]
    sequence <- start:end
    for (k in sequence) {
      AA_df$covered[k] <- T
    }
  }
}

percent_coverage_8_12 <- (sum(AA_df$covered) / nrow(AA_df)) * 100
#print results
cat(paste0("Results for OVCAR3 Sample mapped to OVCAR3 consensus:", "\n", 
           "MUC16 peptides: ", nrow(Ovcar_ov3_MUC16), "\n",
           "MUC16 sequences: ", length(unique(Ovcar_ov3_MUC16$sequence)), "\n",
           "R8 - R12 peptides: ", nrow(R8_12), "\n",
           "R8 - R12 sequences: ", length(unique(R8_12$sequence)), "\n",
           "R8 - R12 percent coverage: ", round(percent_coverage_8_12, 2),"%", "\n",
           "sequences unique to R8-12: ", length(peptide_list_final_Ovcar3_ov3), "\n"))
peptide_list_final_Ovcar3_ov3

write.table(peptide_list_final_Ovcar3_ov3, "unique_repeat_peptides_OVCAR3_IP_ov3_db.txt",
            row.names = F, col.names = F, quote = F)

####plot creation####
####Creates the plots to map peptides onto MUC16 repeats. Color R8-12 a different
####color, and color peptides unique to R8-12 a different color
#fitz mapped to 6s concensus
repeats_6s <- read.csv("Consensus sequences by line v2.txt", header = F)
repeats_6s <- repeats_6s$V1
Fitz_6s_MUC16 <- Fitz_6s_MUC16[Fitz_6s_MUC16$count_proteome == 0,]

sequence <- list()
repeat_n <- list()
start <- list()
end <- list()
unique_8to12 <- list()

for (i in 1:nrow(Fitz_6s_MUC16)) {
  for (j in 1:length(repeats_6s)) {
    matches <- data.frame(str_locate_all(repeats_6s[j], Fitz_6s_MUC16$sequence[i]))
    if (nrow(matches > 0)) {
      sequence[[length(sequence) + 1]] <- Fitz_6s_MUC16$sequence[i]
      repeat_n[[length(repeat_n) + 1]] <- j
      start[[length(start) + 1]] <- matches$start[1]
      end[[length(end) + 1]] <- matches$end[1]
      unique_8to12[[length(unique_8to12) + 1]] <- Fitz_6s_MUC16$unique_to_repeat[i]
    }
  }
}

matches_df <- data.frame(sequence = unlist(sequence),
                         repeat_n = unlist(repeat_n),
                         start = unlist(start),
                         end = unlist(end),
                         unique = unlist(unique_8to12))
matches_df$inv_repeat <- (19 - matches_df$repeat_n) + 1

aa_count <- rep(NA, 19)
for (i in 1:length(repeats_6s)) {
  aa_count[i] <- nchar(repeats_6s[i])
}

repeat_category <- c(rep("old", 7), rep("new", 5), rep("old", 7))

repeat_df <- data.frame(repeat_number = 1:19, 
                       aa_count = aa_count,
                       repeat_category = repeat_category)
repeat_df_0 <- repeat_df
repeat_df_0$aa_count <- 0

repeat_df <- bind_rows(repeat_df, repeat_df_0)
repeat_df$repeat_number <- factor(repeat_df$repeat_number, levels = 19:1)

fitz_6s_plot <- ggplot() + 
  geom_line(data = repeat_df, 
            aes(x = aa_count, y = repeat_number, color = repeat_category, group = repeat_number),
            size = 2) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(y = "Repeat Number", x = "Amino Acid Position in Repeat",
       title = element_blank()) +
  scale_color_manual(values = c("#04d635", "grey30")) +
  geom_rect(data = matches_df, 
            aes(xmin = start, xmax = end, ymin = inv_repeat - 0.3, ymax = inv_repeat + 0.2,
                fill = unique)) +
  scale_fill_manual(values = c("#6abfeb", "#fc03d3"))
fitz_6s_plot 
ggsave("FitzCoveragePlot.png", width = 8, height = 6)


#ovcar3 mapped to Ovcar3 concensus
repeats_ov3 <- read.csv("Ovcar3 seq fasta by line v2.txt", header = F)
repeats_ov3 <- repeats_ov3$V1
Ovcar_ov3_MUC16 <- Ovcar_ov3_MUC16[Ovcar_ov3_MUC16$count_proteome == 0,]

sequence <- list()
repeat_n <- list()
start <- list()
end <- list()
unique_8to12 <- list()

for (i in 1:nrow(Ovcar_ov3_MUC16)) {
  for (j in 1:length(repeats_ov3)) {
    matches <- data.frame(str_locate_all(repeats_ov3[j], Ovcar_ov3_MUC16$sequence[i]))
    if (nrow(matches > 0)) {
      sequence[[length(sequence) + 1]] <- Ovcar_ov3_MUC16$sequence[i]
      repeat_n[[length(repeat_n) + 1]] <- j
      start[[length(start) + 1]] <- matches$start[1]
      end[[length(end) + 1]] <- matches$end[1]
      unique_8to12[[length(unique_8to12) + 1]] <- Ovcar_ov3_MUC16$unique_to_repeat[i]
    }
  }
}

matches_df <- data.frame(sequence = unlist(sequence),
                         repeat_n = unlist(repeat_n),
                         start = unlist(start),
                         end = unlist(end),
                         unique = unlist(unique_8to12))
matches_df$inv_repeat <- (19 - matches_df$repeat_n) + 1

aa_count <- rep(NA, 19)
for (i in 1:length(repeats_ov3)) {
  aa_count[i] <- nchar(repeats_ov3[i])
}

repeat_category <- c(rep("old", 7), rep("new", 5), rep("old", 7))

repeat_df <- data.frame(repeat_number = 1:19, 
                        aa_count = aa_count,
                        repeat_category = repeat_category)
repeat_df_0 <- repeat_df
repeat_df_0$aa_count <- 0

repeat_df <- bind_rows(repeat_df, repeat_df_0)
repeat_df$repeat_number <- factor(repeat_df$repeat_number, levels = 19:1)

Ovcar_ov3_plot <- ggplot() + 
  geom_line(data = repeat_df, 
            aes(x = aa_count, y = repeat_number, color = repeat_category, group = repeat_number),
            size = 2) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(y = "Repeat Number", x = "Amino Acid Position in Repeat",
       title = element_blank()) +
  scale_color_manual(values = c("#04d635", "grey30")) +
  geom_rect(data = matches_df, 
            aes(xmin = start, xmax = end, ymin = inv_repeat - 0.3, ymax = inv_repeat + 0.2,
                fill = unique)) +
  scale_fill_manual(values = c("#6abfeb", "#fc03d3"))
Ovcar_ov3_plot  
ggsave("ovcar3CoveragePlot.png", width = 8, height = 6)
