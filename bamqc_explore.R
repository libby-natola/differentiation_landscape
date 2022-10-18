### explore the bam file quality, coverage 

setwd('Documents/UBC/Bioinformatics/wgs/')
data <- read.csv("bamqc_fixed.csv")

library(ggplot2)
library(tidyr)
library(dplyr)

data$contig_3prop <- data$contig_3/data$MeanCoverage
data$contig_30prop <- data$contig_30/data$MeanCoverage
data$contig_93prop <- data$contig_93/data$MeanCoverage
data$contig_138prop <- data$contig_138/data$MeanCoverage
data$contig_135prop <- data$contig_135/data$MeanCoverage
data$contig_190prop <- data$contig_190/data$MeanCoverage
data$contig_79prop <- data$contig_79/data$MeanCoverage
data$contig_76prop <- data$contig_76/data$MeanCoverage


ggplot(data, aes(x=contig_3prop)) + geom_histogram(aes(color=Sex))
ggplot(data, aes(x=contig_30prop)) + geom_histogram(aes(color=Sex))
ggplot(data, aes(x=contig_93prop)) + geom_histogram(aes(color=Sex))
ggplot(data, aes(x=contig_138prop)) + geom_histogram(aes(color=Sex))
ggplot(data, aes(x=contig_135prop)) + geom_histogram(aes(color=Sex))
ggplot(data, aes(x=contig_190prop)) + geom_histogram(aes(color=Sex))



genomic_F <- data %>% filter(contig_3prop < 0.6)
genomic_M <- data %>% filter(contig_3prop > 0.81)
genomic_huh <- data %>% filter(contig_3prop > 0.6 & contig_3prop < 0.8)

genomic_F_93 <- data %>% filter(contig_93prop > 0.5)
genomic_M_93 <- data %>% filter(contig_93prop < 0.2)
genomic_huh_93 <- data %>% filter(contig_93prop < 0.5 & contig_93prop > 0.2)

genomic_F_30 <- data %>% filter(contig_30prop < 0.6)
genomic_M_30 <- data %>% filter(contig_30prop > 0.8)
genomic_huh_30 <- data %>% filter(contig_30prop < 0.8 & contig_30prop > 0.6)

genomic_F_135 <- data %>% filter(contig_135prop > 1)
genomic_M_135 <- data %>% filter(contig_135prop < 1)
genomic_huh_135 <- data %>% filter(contig_135prop > 10)

genomic_F_79 <- data %>% filter(contig_79prop > 0.6)
genomic_M_79 <- data %>% filter(contig_79prop < 0.6)

genomic_F_76 <- data %>% filter(contig_76prop < 0.6)
genomic_M_76 <- data %>% filter(contig_76prop > 0.6)


### write an if statement to make a new column where if the proportion on 135 is more than 1 write F and if it is less write M 


sex_palette <- c("black", "sky blue")
species_palette <- c("darkred", "red", "purple", "blue", "green", "brown", "white", "yellow", "orange")
ggplot(data2, aes(x=tig00025419)) + geom_histogram(aes(color=genomic_sex))

ggplot(data2, aes(x=tig00007935)) + geom_histogram(aes(fill=Species)) + scale_fill_manual(values=species_palette)

#list of Z chromosomes: 3, 30, 44(weird), 48, 59, 73, 76, 80, 82, 95, 97, 99, 100, 144, 227(?), tig00029968
#list of W chromosomes: 50, 64(?), 71(?), 79, 81, 84, 93, 102, 104, 108(?), 117, 119, 135, 138, 154, 190, 193, 246, 247, 265, 288(WEIRD TRIMODAL), 315, 346, 363, tig00034263, 00007238, 00006789, tig00007920, 00034315, tig00007744, tig00007083, tig00007976(WEIRD TRIMODAL), tig00007935 (WEIRD TRIMODAL) , tig00034407, tig00034174, tig 00161802,
# 196, 257 really weird, check that one out (suspect WISA)
### want to write a script that gets the standardized coverage for all contigs
# take column i/column 5 

data2 <- data
for(i in 8:276){
  data2[, i] <- data2[, i] / data2[, 5]
}

### write a script that makes a new column (genomic sex) that gives genomic sex of
data2$genomic_sex <- if_else(data2$contig_135 < 1, "M", "F", missing = NULL) 


### change the colors of the wisa male and femal
data2[64,278] = "wM"
data2[68,278] = "wF"
ggplot(data2, aes(x=contig_1)) + geom_histogram(aes(color=genomic_sex))

### now I want to loop all the contigs through ggplot to make histograms colored by sex to ID the sex chromosomes

### loop the ggplot for columnns. So I want ggplot to save a histogram of coverage color coded by genomic sex for each column 8-276 in data2, make a loop that cycles through each numbered column
contigs <- print(colnames(data2[,8:276]))
contigs <- noquote(contigs)

for (i in 8:276){
 hist <- ggplot(data2, aes(x=i)) + geom_histogram(aes(color=genomic_sex))
 ggsave(hist, file=paste0("hist_", i,".png"), width = 14, height = 10, units = "cm")
}

# eh I prefer one to have the contig names than the column names
for (i in contigs) {
  x <- print(i, quote = FALSE)
    hist <- ggplot(data2, aes_string(x=x)) + geom_histogram(aes(color=genomic_sex))
      ggsave(hist, file=paste0("hist_", x,".png"), width = 14, height = 10, units = "cm")
}

# remove the wisa to get better bins on the really wild ones
data3 <- data2[-c(64, 68),]

for (i in contigs) {
  x <- print(i, quote = FALSE)
  hist <- ggplot(data3, aes_string(x=x)) + geom_histogram(aes(color=genomic_sex))
  ggsave(hist, file=paste0("hist_", x,".png"), width = 14, height = 10, units = "cm")
}


# now I'll do the species designations
for (i in contigs) {
  x <- print(i, quote = FALSE)
  hist <- ggplot(data2, aes_string(x=x)) + geom_histogram(aes(fill=Species)) + scale_fill_manual(values=species_palette)
  ggsave(hist, file=paste0("spp_hist_", x,".png"), width = 14, height = 10, units = "cm")
}

ggplot(data2, aes(x=tig00007935)) + geom_histogram(aes(fill=Species)) + scale_fill_manual(values=species_palette)


### intermediates I THINK THESE FILTERS ARE WRONG
IF02S01, QD21RABM01, QF16DMNS01 contig 1
GG13MVZ03, HF14AMN01, HG15MVZ01 contig 2
IF02S01, KE27S02, QF16DMNS01 contig 3

data2 %>% filter(contig_2 < 1)
