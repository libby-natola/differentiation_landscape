### make PCAs and manhattans
R

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#* IN R #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#install.packages("qqman")
library(qqman)

#windowed versions

YBSA_RNSA_25kb_fst <- read.table("summary_stats/YBSA-RNSA_25kb.windowed.weir.fst.qqman", header=TRUE)
YBSA_RNSA_25kb_fstsubset <- YBSA_RNSA_25kb_fst[complete.cases(YBSA_RNSA_25kb_fst),]
YBSA_RNSA_25kb_SNP<-c(1:(nrow(YBSA_RNSA_25kb_fstsubset)))
YBSA_RNSA_25kb_df<-data.frame(YBSA_RNSA_25kb_SNP,YBSA_RNSA_25kb_fstsubset)

png(file = "R_files/YBSA_RNSA_25kb_man.png")
manhattan(YBSA_RNSA_25kb_df,chr="CHR",bp="BP",p="P",logp=FALSE,ylab="Weir and Cockerham Fst", col = c("#2f5734", "gray40"))
dev.off()

YBSA_RNSA_50kb_fst <- read.table("summary_stats/YBSA-RNSA_50kb.windowed.weir.fst.qqman", header=TRUE)
YBSA_RNSA_50kb_fstsubset <- YBSA_RNSA_50kb_fst[complete.cases(YBSA_RNSA_50kb_fst),]
YBSA_RNSA_50kb_SNP<-c(1:(nrow(YBSA_RNSA_50kb_fstsubset)))
YBSA_RNSA_50kb_df<-data.frame(YBSA_RNSA_50kb_SNP,YBSA_RNSA_50kb_fstsubset)

png(file = "R_files/YBSA_RNSA_50kb_man.png")
manhattan(YBSA_RNSA_50kb_df,chr="CHR",bp="BP",p="P",logp=FALSE,ylab="Weir and Cockerham Fst", col = c("#2f5734", "gray40"))
dev.off()


### YBSA-RBSA

YBSA_RBSA_25kb_fst <- read.table("summary_stats/YBSA-RBSA_25kb.windowed.weir.fst.qqman", header=TRUE)
YBSA_RBSA_25kb_fstsubset <- YBSA_RBSA_25kb_fst[complete.cases(YBSA_RBSA_25kb_fst),]
YBSA_RBSA_25kb_SNP<-c(1:(nrow(YBSA_RBSA_25kb_fstsubset)))
YBSA_RBSA_25kb_df<-data.frame(YBSA_RBSA_25kb_SNP,YBSA_RBSA_25kb_fstsubset)

png(file = "R_files/YBSA_RBSA_25kb_man.png")
manhattan(YBSA_RBSA_25kb_df,chr="CHR",bp="BP",p="P",logp=FALSE,ylab="Weir and Cockerham Fst",col=c("#c6924d", "gray40"))
dev.off()

YBSA_RBSA_50kb_fst <- read.table("summary_stats/YBSA-RBSA_50kb.windowed.weir.fst.qqman", header=TRUE)
YBSA_RBSA_50kb_fstsubset <- YBSA_RBSA_50kb_fst[complete.cases(YBSA_RBSA_50kb_fst),]
YBSA_RBSA_50kb_SNP<-c(1:(nrow(YBSA_RBSA_50kb_fstsubset)))
YBSA_RBSA_50kb_df<-data.frame(YBSA_RBSA_50kb_SNP,YBSA_RBSA_50kb_fstsubset)

png(file = "R_files/YBSA_RBSA_50kb_man.png")
manhattan(YBSA_RBSA_50kb_df,chr="CHR",bp="BP",p="P",logp=FALSE,ylab="Weir and Cockerham Fst",col=c("#c6924d", "gray40"))
dev.off()


### RBSA-RNSA
RBSA_RNSA_25kb_fst <- read.table("summary_stats/RBSA-RNSA_25kb.windowed.weir.fst.qqman", header=TRUE)
RBSA_RNSA_25kb_fstsubset <- RBSA_RNSA_25kb_fst[complete.cases(RBSA_RNSA_25kb_fst),]
RBSA_RNSA_25kb_SNP<-c(1:(nrow(RBSA_RNSA_25kb_fstsubset)))
RBSA_RNSA_25kb_df<-data.frame(RBSA_RNSA_25kb_SNP,RBSA_RNSA_25kb_fstsubset)

png(file = "R_files/RBSA_RNSA_25kb_man.png")
manhattan(RBSA_RNSA_25kb_df,chr="CHR",bp="BP",p="P",logp=FALSE,ylab="Weir and Cockerham Fst",col=c("#b26ec4", "gray40"))
dev.off()

RBSA_RNSA_50kb_fst <- read.table("summary_stats/RBSA-RNSA_50kb.windowed.weir.fst.qqman", header=TRUE)
RBSA_RNSA_50kb_fstsubset <- RBSA_RNSA_50kb_fst[complete.cases(RBSA_RNSA_50kb_fst),]
RBSA_RNSA_50kb_SNP<-c(1:(nrow(RBSA_RNSA_50kb_fstsubset)))
RBSA_RNSA_50kb_df<-data.frame(RBSA_RNSA_50kb_SNP,RBSA_RNSA_50kb_fstsubset)

png(file = "R_files/RBSA_RNSA_50kb_man.png")
manhattan(RBSA_RNSA_50kb_df,chr="CHR",bp="BP",p="P",logp=FALSE,ylab="Weir and Cockerham Fst",col=c("#b26ec4", "gray40"))
dev.off()


############################


#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("SNPRelate")
#biocLite("SNPRelate", dependencies = TRUE)

#install.packages("tidyverse")
#install.packages("data.table")
#install.packages("grid")
#install.packages("gridExtra")


library(SNPRelate)
library(tidyverse)
library(data.table)
library(grid)
library(gridExtra)


## load SNP data in VCF format
snpgdsVCF2GDS(vcf.fn="snps/wgs_SNPs_filtered_missing80_mindepth3.recode.chroms_renamed.vcf", out.fn="snps/wgs_SNPs_filtered_missing80_mindepth3.vcf.recode.chroms_renamed.gds", method = c("copy.num.of.ref"),compress.annotation="ZIP.max", snpfirstdim=FALSE, verbose=TRUE)

snpgdsVCF2GDS(vcf.fn="snps/wgs_SNPs_filtered_missing80_mindepth3_nowisa.recode.chroms_renamed_num.vcf", out.fn="snps/wgs_SNPs_filtered_missing80_mindepth3_nowisa.recode.chroms_renamed_num.gds", method = c("copy.num.of.ref"),compress.annotation="ZIP.max", snpfirstdim=FALSE, verbose=TRUE)

## summarize input file
snpgdsSummary("snps/wgs_SNPs_filtered_missing80_mindepth3.vcf.recode.chroms_renamed.gds")

snpgdsSummary("snps/wgs_SNPs_filtered_missing80_mindepth3_nowisa.recode.chroms_renamed_num.gds")


## open file
genofile <- snpgdsOpen("snps/wgs_SNPs_filtered_missing80_mindepth3.vcf.recode.chroms_renamed.gds")
read.gdsn(index.gdsn(genofile, "sample.id"))
genofile_nowisa <- snpgdsOpen("snps/wgs_SNPs_filtered_missing80_mindepth3_nowisa.recode.chroms_renamed_num.gds")
read.gdsn(index.gdsn(genofile_nowisa, "sample.id"))

## add info on samples (sample ID, taxa ID, locality, pheno scores, etc.)
sample_info <- read.table("extras/sample_info.txt", sep="\t", header=TRUE)
sample_info_nowisa <- read.table("extras/sample_info_nowisa.txt", sep="\t", header=TRUE)

## assess missing data for each sample
miss <- snpgdsSampMissRate(genofile, sample.id=NULL, snp.id=NULL, with.id=TRUE)
miss <- as.data.frame(miss)
miss <- setDT(miss,keep.rownames=TRUE)[]
colnames(miss) <- c("ID", "missing")
miss_merge <- merge(miss, sample_info, by="ID")
miss_output <- select(miss_merge, samplenum, taxa, missing)

miss_nowisa <- snpgdsSampMissRate(genofile_nowisa, sample.id=NULL, snp.id=NULL, with.id=TRUE)
miss_nowisa <- as.data.frame(miss_nowisa)
miss_nowisa <- setDT(miss_nowisa,keep.rownames=TRUE)[]
colnames(miss_nowisa) <- c("ID", "missing")
miss_merge_nowisa <- merge(miss_nowisa, sample_info_nowisa, by="ID")
miss_output_nowisa <- select(miss_merge_nowisa, samplenum, taxa, missing)

## write missing data to file
write.table(miss_output,"snps/wgs_SNPs_filtered_missing80_mindepth3.vcf.recode.chroms_renamed.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(miss_output_nowisa,"snps/wgs_SNPs_filtered_missing80_mindepth3_nowisa.recode.chroms_renamed_num.txt", sep="\t", quote=FALSE, row.names=TRUE)



########### PCA ###########

## run PCA using SNPRelate
pca_wisa <- snpgdsPCA(gdsobj = genofile, autosome.only=FALSE, sample.id=NULL)
pca_nowisa <- snpgdsPCA(gdsobj = genofile_nowisa, autosome.only=FALSE, sample.id=NULL)

## get percent variation explained for each PC axis
pc.percent_wisa <- pca_wisa$varprop*100
head(round(pc.percent_wisa, 2))
# [1] 14.60  4.24  3.05  2.42  1.51  1.47
pc.percent_nowisa <- pca_nowisa$varprop*100
head(round(pc.percent_nowisa, 2))
# [1] 14.77  3.17  2.54  1.60  1.56  1.55

## pull sample ID + first four PC axes
pca_coords_wisa <- data.frame(ID = pca_wisa$sample.id,
                         pc1 = pca_wisa$eigenvect[,1],    # the first eigenvector
                         pc2 = pca_wisa$eigenvect[,2],    # the second eigenvector
                         pc3 = pca_wisa$eigenvect[,3],
                         pc4 = pca_wisa$eigenvect[,4],
                         stringsAsFactors = FALSE)
head(pca_coords_wisa)

pca_coords_nowisa <- data.frame(ID = pca_nowisa$sample.id,
                         pc1 = pca_nowisa$eigenvect[,1],    # the first eigenvector
                         pc2 = pca_nowisa$eigenvect[,2],    # the second eigenvector
                         pc3 = pca_nowisa$eigenvect[,3],
                         pc4 = pca_nowisa$eigenvect[,4],
                         stringsAsFactors = FALSE)
head(pca_coords_nowisa)


## merge PCA results with sample info by ID number
pca_coords_merged_wisa <- merge(pca_coords_wisa, sample_info, by.x="ID")
pca_coords_merged_nowisa <- merge(pca_coords_nowisa, sample_info_nowisa, by.x="ID")

## figure colors
fig_colors_wisa <- c("red4", "red", "purple", "blue", "green4", "black", "gray50", "gold", "darkorange")
fig_colors_nowisa <- c("#bd7777", "#bd5757", "#b26ec4", "#394ca2", "#2f5734", "grey25", "#e6e045", "#c6924d")
fig_shapes_nowisa <- c(16, 19, 5, 15, 6, 13, 17, 5)

## scatterplot of PC1 versus PC2
pdf(file = "R_files/wgs_pca_1_2_wisa.pdf")
ggplot() +
  geom_hline(aes(yintercept=0), color="gray") +
  geom_vline(aes(xintercept=0), color="gray") +
  geom_point(data=pca_coords_merged_wisa, aes(x=pc1, y=pc2, fill=taxa), size=4, alpha=0.75, shape=21, stroke=0.2) +
  labs(x="PC1 (14.60%)", y="PC2 (4.24%)") +
  scale_fill_manual(values=fig_colors_wisa) +
  theme_classic() +
  theme(legend.position="none", axis.line=element_line(color="black"), axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))
dev.off()

pdf(file = "R_files/wgs_pca_1_2_nowisa.pdf")
ggplot() +
  geom_hline(aes(yintercept=0), color="gray") +
  geom_vline(aes(xintercept=0), color="gray") +
  geom_point(data=pca_coords_merged_nowisa, aes(x=pc1, y=pc2, fill=taxa), size=4, alpha=0.75, shape=21, stroke=0.2) +
  labs(x="PC1 (14.77%)", y="PC2 (3.17%)") +
  scale_fill_manual(values=fig_colors_nowisa) +
  theme_classic() +
  theme(legend.position="none", axis.line=element_line(color="black"), axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))
dev.off()

pdf(file = "R_files/wgs_pca_1_2_nowisa_shapes.pdf", width = 9, height= 7)
ggplot() +
  geom_hline(aes(yintercept=0), color="gray") +
  geom_vline(aes(xintercept=0), color="gray") +
  geom_point(data=pca_coords_merged_nowisa, aes(x=pc1, y=pc2, color=taxa, shape = taxa), size=4, alpha=0.75, stroke=0.2) +
  labs(x="PC1 (14.77%)", y="PC2 (3.17%)") +
  scale_color_manual(name = "Species", labels = c("Red-breasted daggetti", "Red-breasted ruber", "Red-naped x Red-breasted ruber", "Red-naped", "Red-naped x Yellow-bellied", "Red-breasted x Red-naped x Yellow-bellied", "Yellow-bellied", "Yellow-bellied x Red-breasted ruber"), values=fig_colors_nowisa) +
  scale_shape_manual(name = "Species", labels = c("Red-breasted daggetti", "Red-breasted ruber", "Red-naped x Red-breasted ruber", "Red-naped", "Red-naped x Yellow-bellied", "Red-breasted x Red-naped x Yellow-bellied", "Yellow-bellied", "Yellow-bellied x Red-breasted ruber"), values=c(16, 19, 5, 15, 6, 13, 17, 0)) +
  theme_classic() +
  theme(legend.position="right", axis.line=element_line(color="black"), axis.title=element_text(size=12), axis.text=element_text(size=12), aspect.ratio=1)
dev.off()

## scatterplot of PC3 versus PC4
pdf(file = "R_files/wgs_pca_3_4_wisa.pdf")
ggplot() +
  geom_hline(aes(yintercept=0), color="gray") +
  geom_vline(aes(xintercept=0), color="gray") +
  geom_point(data=pca_coords_merged_wisa, aes(x=pc3, y=pc4, fill=taxa), size=4, alpha=0.75, shape=21, stroke=0.2) +
  labs(x="PC3 (3.05%)", y="PC4 (2.42%)") +
  scale_fill_manual(values=fig_colors_wisa) +
  theme_classic() +
  theme(legend.position="none", axis.line=element_line(color="black"), axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))
dev.off()

pdf(file = "R_files/wgs_pca_3_4_nowisa.pdf")
ggplot() +
  geom_hline(aes(yintercept=0), color="gray") +
  geom_vline(aes(xintercept=0), color="gray") +
  geom_point(data=pca_coords_merged_nowisa, aes(x=pc3, y=pc4, fill=taxa), size=4, alpha=0.75, shape=21, stroke=0.2) +
  labs(x="PC3 (2.54%)", y="PC4 (1.60%)") +
  scale_fill_manual(values=fig_colors_nowisa) +
  theme_classic() +
  theme(legend.position="none", axis.line=element_line(color="black"), axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))
dev.off()


