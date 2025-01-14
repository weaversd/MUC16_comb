library(seqinr)
library(ggplot2)
library(stringr)
library(dplyr)


setwd("~/03_Notre_Dame/01_Dissertation/Chapter-B Ovarian Cancer Method Development and Data analyisis/Code/Deglyco")

repeatFasta <- read.fasta("Consensus seq fasta by repeat v2 for deglyco.fasta", seqtype = "AA")

file <- "Area_R04_intensity_by_aminoacid_two_sample.csv"




dfRepeat <- read.csv(file)
repeatTitle <- str_extract(file, "R[0-9]+")
repeatSequence <- c2s(repeatFasta[[repeatTitle]])

dfRepeat$RepeatN <- as.numeric(str_remove(str_extract(file, "R[0-9]+"), "R"))
dfRepeat$Deamidation <- "none"

for (sampleType in unique(dfRepeat$sample)) {
  print(sampleType)
  uniquePepsInSample <- unique(dfRepeat$origin_pep[dfRepeat$sample == sampleType &
                                                     dfRepeat$Area > 0])
  uniquePepsInSample <- uniquePepsInSample[uniquePepsInSample != ""]
  uniquePepsInSample <- unique(unlist(strsplit(uniquePepsInSample, ";")))
  for (peptide in uniquePepsInSample) {
    deamidationPeptide <- str_replace_all(peptide, "\\(\\+0\\.98\\)", "{deam}")
    deamidationPeptide <- str_remove_all(deamidationPeptide, "\\([+-][0-9]+\\.[0-9]+\\)")
    deamidationCount <- str_count(deamidationPeptide, "\\{deam\\}")
    if (deamidationCount > 0) {
      if (deamidationCount > 1) {
        print("multiple deamidations")
        print(deamidationPeptide)
        print(repeatTitle)
      }
      deamidationIndexInPeptide <- gregexpr("\\{deam\\}", deamidationPeptide)[[1]][1] -1
      peptideIndexInRepeat <- gregexpr(str_remove_all(peptide, "\\([+-][0-9]+\\.[0-9]+\\)"), repeatSequence)[[1]][1]
      dfRepeat$Deamidation[dfRepeat$sample == sampleType &
                             dfRepeat$AA_index == peptideIndexInRepeat + deamidationIndexInPeptide - 1] <- "deamidation"
    }
  }
}

glycosites <- data.frame(index = gregexpr("[N].[ST]", repeatSequence)[[1]],
                         repeatN = as.numeric(str_remove(str_extract(file, "R[0-9]+"), "R")))





ggplot(dfRepeat) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_line(aes(x = AA_index, y = Area, color = sample)) +
  geom_point(data = dfRepeat[dfRepeat$Deamidation == "deamidation",],
             aes(x = AA_index, y = Area), size = 3, color = "hotpink") +
  geom_vline(data = glycosites, aes(xintercept = index), color = "green4") +
  facet_grid(RepeatN~.)

























fileList <- list.files(pattern = "Area_R.*.csv")
repeatDFList <- list()
glycositeDFList <- list()

for (file in fileList) {
  dfRepeat <- read.csv(file)
  repeatTitle <- str_extract(file, "R[0-9]+")
  repeatSequence <- c2s(repeatFasta[[repeatTitle]])
  
  dfRepeat$RepeatN <- as.numeric(str_remove(str_extract(file, "R[0-9]+"), "R"))
  dfRepeat$Deamidation <- "none"
  
  for (sampleType in unique(dfRepeat$sample)) {
    uniquePepsInSample <- unique(dfRepeat$origin_pep[dfRepeat$sample == sampleType &
                                                       dfRepeat$Area > 0])
    uniquePepsInSample <- uniquePepsInSample[uniquePepsInSample != ""]
    uniquePepsInSample <- uniquePepsInSample[uniquePepsInSample != ""]
    uniquePepsInSample <- unique(unlist(strsplit(uniquePepsInSample, ";")))
    for (peptide in uniquePepsInSample) {
      deamidationPeptide <- str_replace_all(peptide, "\\(\\+0\\.98\\)", "{deam}")
      deamidationPeptide <- str_remove_all(deamidationPeptide, "\\([+-][0-9]+\\.[0-9]+\\)")
      deamidationCount <- str_count(deamidationPeptide, "\\{deam\\}")
      if (deamidationCount > 0) {
        if (deamidationCount > 1) {
          print("multiple deamidations")
          print(deamidationPeptide)
          print(repeatTitle)
        }
        deamidationIndexInPeptide <- gregexpr("\\{deam\\}", deamidationPeptide)[[1]][1] -1
        peptideIndexInRepeat <- gregexpr(str_remove_all(peptide, "\\([+-][0-9]+\\.[0-9]+\\)"), repeatSequence)[[1]][1]
        dfRepeat$Deamidation[dfRepeat$sample == sampleType &
                               dfRepeat$AA_index == peptideIndexInRepeat + deamidationIndexInPeptide - 1] <- "deamidation"
      }
    }
  }
  
  glycosites <- data.frame(index = gregexpr("[N].[ST]", repeatSequence)[[1]],
                           RepeatN = as.numeric(str_remove(str_extract(file, "R[0-9]+"), "R")))
  
  repeatDFList[[file]] <- dfRepeat
  glycositeDFList[[file]] <- glycosites
  
}


dfRepeat <- bind_rows(repeatDFList)
glycosites <- bind_rows(glycositeDFList)

dfRepeat$AreaLog <- log(dfRepeat$Area, base = 10)
dfRepeat$AreaLog2 <- log(dfRepeat$Area, base = 2)
dfRepeat$AreaLog[!is.finite(dfRepeat$AreaLog)] <- 0
dfRepeat$AreaLog2[!is.finite(dfRepeat$AreaLog2)] <- 0

ymax <- max(dfRepeat$Area,na.rm = T)

ggplot() +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  geom_vline(data = glycosites, aes(xintercept = index), color = "purple", size = 1) +
  geom_line(data = dfRepeat, aes(x = AA_index, y = AreaLog, color = sample), linewidth = 1) +
  geom_point(data = dfRepeat[dfRepeat$Deamidation == "deamidation",],
             aes(x = AA_index, y = AreaLog, fill = sample),
             size = 3, shape = 21, stroke = 2, alpha = 0.8)+
  facet_grid(RepeatN~.) +
  labs(x = "AA position in Repeat", y = "Peptide Area (Log Transformed)",
       fill = element_blank(), color = element_blank())

ggsave("RepeatDeglycoPlotOverlay.png", width = 8, height = 14)



ggplot() +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none") +
  geom_vline(data = glycosites, aes(xintercept = index), color = "purple", size = 1) +
  geom_line(data = dfRepeat, aes(x = AA_index, y = AreaLog, color = sample), linewidth = 1) +
  geom_point(data = dfRepeat[dfRepeat$Deamidation == "deamidation",],
             aes(x = AA_index, y = AreaLog, fill = sample),
             size = 3, shape = 21, stroke = 2, alpha = 0.8)+
  facet_grid(RepeatN~sample) +
  labs(x = "AA position in Repeat", y = "Peptide Area (Log Transformed)",
       fill = element_blank(), color = element_blank())


ggsave("RepeatDeglycoPlotOverlaySplit.png", width = 8, height = 14)




dfRepeatGrouped <- dfRepeat %>% group_by(RepeatN, AA, AA_index) %>%
  summarise(zeroCount = sum(Area == 0),
            LFC = AreaLog2[sample == "PNGaseF"] - AreaLog2[sample == "CTRL"],
            PNGaseDEAM = Deamidation[sample == "PNGaseF"],
            CTRLDEAM = Deamidation[sample == "CTRL"])

dfRepeatGrouped$infiniteYpos <- NA
dfRepeatGrouped$infiniteYpos[dfRepeatGrouped$zeroCount == 1 &
                               dfRepeatGrouped$LFC > 0] <- dfRepeatGrouped$LFC[dfRepeatGrouped$zeroCount == 1 &
                                                                                 dfRepeatGrouped$LFC > 0]
dfRepeatGrouped$infiniteYneg <- NA
dfRepeatGrouped$infiniteYneg[dfRepeatGrouped$zeroCount == 1 &
                               dfRepeatGrouped$LFC < 0] <- dfRepeatGrouped$LFC[dfRepeatGrouped$zeroCount == 1 &
                                                                                 dfRepeatGrouped$LFC < 0]
                                                                                 

ggplot() +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  geom_vline(data = glycosites, aes(xintercept = index), color = "purple", size = 1) +
  geom_line(data = dfRepeatGrouped, aes(x = AA_index, y = LFC), linewidth = 1) +
  facet_grid(RepeatN~.) +
  labs(x = "AA position in Repeat", y = "Peptide Area Log 2 Fold Change (PNGaseF / CTRL)",
       fill = element_blank(), color = element_blank()) +
  geom_hline(yintercept = 0, color = "pink", linetype = 3, linewidth = 0.5) +
  geom_line(data = dfRepeatGrouped, aes(x = AA_index, y = infiniteYpos),
            linewidth = 1, color = "green") +
  geom_line(data = dfRepeatGrouped, aes(x = AA_index, y = infiniteYneg),
            linewidth = 1, color = "red") +
  geom_point(data = dfRepeatGrouped[dfRepeatGrouped$PNGaseDEAM == "deamidation",],
             aes(x = AA_index, y = LFC),
             size = 3, shape = 21, stroke = 2, fill = "#00BFC4") +
  geom_point(data = dfRepeatGrouped[dfRepeatGrouped$CTRLDEAM == "deamidation",],
             aes(x = AA_index, y = LFC),
             size = 3, shape = 21, stroke = 2, fill = "#F8766D")


ggsave("RepeatDeglycoPlotLFC.png", width = 8, height = 14)
