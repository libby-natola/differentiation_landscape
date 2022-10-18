## clear workspace
unlink(".RData")

#set huge memory
library(unix)
rlimit_as(Inf)

## set working directory
setwd("/scratch/st-darreni-1/libby/WGS/")
library(ggplot2)

## have it take args from command line execution
args <- commandArgs(trailingOnly = TRUE)

chr <- read.table(paste("chrs/wgs_SNPs_filtered_missing80_mindepth3_nowisa.recode.vcf.renamed_", args[1], ".recode.", args[2], ".vcf.ld", sep = ""), header = TRUE)
ggsave(paste("chrs/heatmaps/", args[1], "_", args[2], ".LD.png", sep = ""))
ggplot(chr, aes(BP_A,BP_B, color=R2))+geom_point()+
  scale_color_gradient(low="white", high="red") + labs(title=paste("LD ", args[1], args[2]))
dev.off()
