setwd("~/03_Notre_Dame/01_Dissertation/Chapter-B Ovarian Cancer Method Development and Data analyisis/Code/CE fractionation")

MUC16_seq_list <- read.fasta("MUC16_22K.fasta", seqtype = 'AA', seqonly = TRUE)

MUC16_seq <- MUC16_seq_list[[1]]

MUC16_seq <- read.csv("MUC16_consensus_from 6_sequences.txt", header = F)$V1
MUC16_AA_count <- nchar(MUC16_seq)

#combined Fractions
CCM_C_peptides <- read.csv("SDW_MUC16_Deglyco.CE_fractionation_CCM_C/db.peptides.csv")
CCM_C_MUC16_peps <- CCM_C_peptides[grep("Q8WXI7", CCM_C_peptides$Accession),]
CCM_C_MUC16_peps_list <- CCM_C_MUC16_peps$Peptide
for (i in 1:length(CCM_C_MUC16_peps_list)) {
  CCM_C_MUC16_peps_list[i] <- remove_mods_str(CCM_C_MUC16_peps_list[i])
}
CCM_C_MUC16_peps_list_unique <- unique(CCM_C_MUC16_peps_list)




CCM_C_MUC16_peps$Fraction_sum_PSMs <- rowSums(CCM_C_MUC16_peps[,grepl("X.Spec.*F[1-5]", names(CCM_C_MUC16_peps))])
CCM_C_MUC16_peps$Raw_sum_PSMs <- rowSums(CCM_C_MUC16_peps[,grepl("X.Spec.*Raw", names(CCM_C_MUC16_peps))])
CCM_C_MUC16_peps$Plug_sum_PSMs <- rowSums(CCM_C_MUC16_peps[,grepl("X.Spec.*Plug", names(CCM_C_MUC16_peps))])



FractionPeptides <- CCM_C_MUC16_peps$Peptide[CCM_C_MUC16_peps$Fraction_sum_PSMs > 0]
for (i in 1:length(FractionPeptides)) {
  FractionPeptides[i] <- remove_mods_str(FractionPeptides[i])
}
FractionPeptides <- unique(FractionPeptides)
write.table(FractionPeptides, "Fractionpeptides.txt", col.names = FALSE, row.names = FALSE)


CCM_C_fracs <- peptide_df_gen('MUC16_consensus_from 6_sequences.fasta', 'Fractionpeptides.txt', 'Fractions')
CCM_C_fracs_positions <- CCM_C_fracs$position
count_CCM_C_fracs_positions <- unique(CCM_C_fracs_positions)
covered_AA_count <- length(count_CCM_C_fracs_positions)


CCM_C_percent_cov <- round((covered_AA_count / MUC16_AA_count) * 100, 2)



RawPeptides <- CCM_C_MUC16_peps$Peptide[CCM_C_MUC16_peps$Raw_sum_PSMs > 0]
for (i in 1:length(RawPeptides)) {
  RawPeptides[i] <- remove_mods_str(RawPeptides[i])
}
RawPeptides <- unique(RawPeptides)
write.table(RawPeptides, "Rawpeptides.txt", col.names = FALSE, row.names = FALSE)


CCM_C_peptides_raw <- peptide_df_gen('MUC16_consensus_from 6_sequences.fasta', 'Rawpeptides.txt', 'Raw')
CCM_C_raw_positions <- CCM_C_peptides_raw$position
count_CCM_C_raw_positions <- unique(CCM_C_raw_positions)
covered_AA_count <- length(CCM_C_raw_positions)


CCM_C_percent_cov_raw <- round((covered_AA_count / MUC16_AA_count) * 100, 2)

PlugPeptides <- CCM_C_MUC16_peps$Peptide[CCM_C_MUC16_peps$Plug_sum_PSMs > 0]
for (i in 1:length(PlugPeptides)) {
  PlugPeptides[i] <- remove_mods_str(PlugPeptides[i])
}
PlugPeptides <- unique(PlugPeptides)
write.table(PlugPeptides, "Plugpeptides.txt", col.names = FALSE, row.names = FALSE)


CCM_C_peptides_plug <- peptide_df_gen('MUC16_consensus_from 6_sequences.fasta', 'Plugpeptides.txt', 'Plug')
CCM_C_plug_positions <- CCM_C_peptides_plug$position
count_CCM_C_plug_positions <- unique(CCM_C_plug_positions)
covered_AA_count <- length(CCM_C_plug_positions)


CCM_C_percent_cov_plug <- round((covered_AA_count / MUC16_AA_count) * 100, 2)



# 
# 
# #Raw
# CCM_C_peptides_raw <- read.csv("CCM_C/db.CCM_C_RAW.peptides.csv")
# CCM_C_MUC16_peps_raw <- CCM_C_peptides_raw[grep("Q8WXI7", CCM_C_peptides_raw$Accession),]
# CCM_C_MUC16_peps_list_raw <- CCM_C_MUC16_peps_raw$Peptide
# for (i in 1:length(CCM_C_MUC16_peps_list_raw)) {
#   CCM_C_MUC16_peps_list_raw[i] <- remove_mods_str(CCM_C_MUC16_peps_list_raw[i])
# }
# CCM_C_MUC16_peps_list_unique_raw <- unique(CCM_C_MUC16_peps_list_raw)
# write.table(CCM_C_MUC16_peps_list_unique_raw, "CCM_C/CCM_C_MUC16_unique_peps_raw.txt",col.names = FALSE, row.names = FALSE)
# 
# 
# CCM_C_fracs_raw <- peptide_df_gen('MUC16_22K.fasta', 'CCM_C/CCM_C_MUC16_unique_peps_raw.txt', 'CCM_C_combined_RAW')
# CCM_C_fracs_positions_raw <- CCM_C_fracs_raw$position
# count_CCM_C_fracs_positions_raw <- unique(CCM_C_fracs_positions_raw)
# covered_AA_count_raw <- length(count_CCM_C_fracs_positions_raw)
# 
# 
# CCM_C_percent_cov_raw <- round((covered_AA_count_raw / MUC16_AA_count) * 100, 2)
# print(paste0("CCM C MUC16 percent coverage (raw): ", CCM_C_percent_cov_raw), quote = F)
# 
# 
# 
# #Plug
# CCM_C_peptides_plug <- read.csv("CCM_C/db.CCM_C_PLUG.peptides.csv")
# CCM_C_MUC16_peps_plug <- CCM_C_peptides_plug[grep("Q8WXI7", CCM_C_peptides_plug$Accession),]
# CCM_C_MUC16_peps_list_plug <- CCM_C_MUC16_peps_plug$Peptide
# for (i in 1:length(CCM_C_MUC16_peps_list_plug)) {
#   CCM_C_MUC16_peps_list_plug[i] <- remove_mods_str(CCM_C_MUC16_peps_list_plug[i])
# }
# CCM_C_MUC16_peps_list_unique_plug <- unique(CCM_C_MUC16_peps_list_plug)
# write.table(CCM_C_MUC16_peps_list_unique_plug, "CCM_C/CCM_C_MUC16_unique_peps_plug.txt",col.names = FALSE, row.names = FALSE)
# 
# 
# CCM_C_fracs_plug <- peptide_df_gen('MUC16_22K.fasta', 'CCM_C/CCM_C_MUC16_unique_peps_plug.txt', 'CCM_C_combined_PLUG')
# CCM_C_fracs_positions_plug <- CCM_C_fracs_plug$position
# count_CCM_C_fracs_positions_plug <- unique(CCM_C_fracs_positions_plug)
# covered_AA_count_plug <- length(count_CCM_C_fracs_positions_plug)
# 
# 
# CCM_C_percent_cov_plug <- round((covered_AA_count_plug / MUC16_AA_count) * 100, 2)
# print(paste0("CCM C MUC16 percent coverage (plug): ", CCM_C_percent_cov_plug), quote = F)



CCM_C_plot_df <- bind_rows(CCM_C_fracs, CCM_C_peptides_raw, CCM_C_peptides_plug)

labels_df <- data.frame(group = c("Fractions",
                                  "Plug",
                                  "Raw"),
                        coverage = c(CCM_C_percent_cov,
                                     CCM_C_percent_cov_plug, 
                                     CCM_C_percent_cov_raw))
labels_df$coverage_str <- as.character(labels_df$coverage)
for (i in 1:nrow(labels_df)) {
  if (nchar(labels_df$coverage_str[i] < 4)) {
    number_char <- nchar(labels_df$coverage_str[i])
    zeros_needed <- 4 - number_char               
    labels_df$coverage_str[i] <- paste0(labels_df$coverage_str[i], rep("0", zeros_needed))
  }
  labels_df$coverage_str[i] <- paste0(labels_df$coverage_str[i], "%")
}


plot <- ggplot(data = CCM_C_plot_df) +
  geom_point(aes(x = position, y = group, color = group), shape = 124,
             size = 5) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none") +
  labs(y = element_blank(), x = "MUC 16 Residue Number") + 
  scale_y_discrete(labels = c("Fractions", "Plug", "Raw")) +
  scale_x_continuous(limits = c(0,17000), 
                     breaks = c(0, 5000, 10000, 15000)) +
  geom_vline(xintercept = 0, color = "red", linetype = 2) +
  geom_vline(xintercept = nchar(MUC16_seq), color = "red", linetype = 2) +
  geom_label(data = labels_df, aes(x = 16500, y = group, label = coverage_str),
             size = 6)
plot

ggsave("CoveragePlot.png", width = 10, height = 3)




repeatFile <- read.csv("Consensus sequences by line v2.txt", header = F)

repeatList <- list()
for (i in 1:nrow(repeatFile)) {
  repeatString <- repeatFile$V1[i]
  for (j in 1:length(FractionPeptides)) {
    peptide <- FractionPeptides[j]
    match <- gregexpr(peptide, repeatString)
    print(peptide)
    startIndex <- match[[1]]
    endIndex <- startIndex - nchar(peptide) - 1
    print(c(startIndex, endIndex))
    if (startIndex >= 0) {
      repeatList[[length(repeatList) + 1]] <- c(i, startIndex, endIndex)
    }
  }
}

coverageDF_fractions <- as.data.frame(do.call(rbind, repeatList))
names(coverageDF_fractions) <- c("repeat", "start", "end")
coverageDF_fractions$sample <- "Fractions"



repeatList <- list()
for (i in 1:nrow(repeatFile)) {
  repeatString <- repeatFile$V1[i]
  for (j in 1:length(PlugPeptides)) {
    peptide <- PlugPeptides[j]
    match <- gregexpr(peptide, repeatString)
    print(peptide)
    startIndex <- match[[1]]
    endIndex <- startIndex - nchar(peptide) - 1
    print(c(startIndex, endIndex))
    if (startIndex >= 0) {
      repeatList[[length(repeatList) + 1]] <- c(i, startIndex, endIndex)
    }
  }
}

coverageDF_plug <- as.data.frame(do.call(rbind, repeatList))
names(coverageDF_plug) <- c("repeat", "start", "end")
coverageDF_plug$sample <- "Plug"



repeatList <- list()
for (i in 1:nrow(repeatFile)) {
  repeatString <- repeatFile$V1[i]
  for (j in 1:length(RawPeptides)) {
    peptide <- RawPeptides[j]
    match <- gregexpr(peptide, repeatString)
    print(peptide)
    startIndex <- match[[1]]
    endIndex <- startIndex - nchar(peptide) - 1
    print(c(startIndex, endIndex))
    if (startIndex >= 0) {
      repeatList[[length(repeatList) + 1]] <- c(i, startIndex, endIndex)
    }
  }
}

coverageDF_raw <- as.data.frame(do.call(rbind, repeatList))
names(coverageDF_raw) <- c("repeat", "start", "end")
coverageDF_raw$sample <- "Raw"


RepeatCoverage <- bind_rows(coverageDF_raw, coverageDF_fractions, coverageDF_plug)
RepeatCoverage$inv_repeat <- (19 - RepeatCoverage$`repeat`) + 1

RepeatCoverage$yoffset <- 0
RepeatCoverage$yoffset[RepeatCoverage$sample == "Fractions"] <- -0.3
RepeatCoverage$yoffset[RepeatCoverage$sample == "Raw"] <- 0.3

repeatlabeldf <- data.frame(y = 1:19, value = 19:1)

repeatLengthDF <- data.frame(`repeat` = 1:19,
                             inv_repeat = 19:1,
                             length = nchar(repeatFile$V1))

ggplot() +
  geom_rect(data = repeatLengthDF, aes(ymin = inv_repeat - 0.4, xmin = 1, xmax = length, ymax = inv_repeat + 0.4),
               fill = "grey", alpha = 0.6, size = 7) +
  geom_rect(data = RepeatCoverage,
            aes(xmin = start, xmax = end, ymin = inv_repeat + yoffset - 0.1, ymax = inv_repeat + 0.1 + yoffset,
                fill = sample)) +
  theme_bw(base_size = 15) +
  scale_y_continuous(breaks = 1:19) +
  coord_cartesian(ylim = c(1,19), xlim = c(.2,155)) +
  geom_hline(yintercept = seq(0.5, 19.5, 1)) +
  geom_label(data = repeatlabeldf, x = -4, aes(y = y, label = value)) +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(fill = element_blank(), y = "Repeat Number",
       x = "Amino Acid in Repeat")

ggsave("CEfractionationCoverageMapofRepeats.png", width = 7, height = 9)

