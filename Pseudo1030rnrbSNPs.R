setwd("/Users/libbynatola/Documents/UBC/Bioinformatics/wgs")

# Load functions
source ("plumage/genomics_R_functions_V2.R")

# install.packages("dplyr") 
library(dplyr) #DI ADDED THIS

# read in sample data
wgs_sample_info <- read.table("sample_info.txt", fill = T, sep="\t", header=T)

# read in genomic data
base.file.name <- "genotype.x.individual/wgs_SNPs_filtered_missing80_mindepth3_nowisa.chroms_renamed_num.1030_rnrb_fixedsnps.males"
pos <- read.table(paste0(base.file.name, ".012.pos"), col.names = c("chrom", "position"))
column_names <- c("null", paste("c", pos$chrom, pos$position, sep="."))
geno <- read.table(paste0(base.file.name, ".012NA"), colClasses = "integer", col.names = column_names)
SNPnum <- length(geno[1,]) -1   # because the first column is not a SNP (just a count from zero)
ind <- read.table(paste0(base.file.name, ".012.indv"))

#get location data from pg_sample_info and ind, so the pg_sample_info is in the right order
colnames(ind) <- "ID"
sample_info_ordered <- left_join(ind, wgs_sample_info) 
colnames(sample_info_ordered) <- c("ID", "MuseumID", "Institution", "group", "location", "date", "lat", "lon", "sex")

# indicate name of metadata file, a text file with these column headings, make the chin phenotypes the different groups:
# ID  location  group Fst_group plot_order
sample_info_ordered$Fst_group <- sample_info_ordered$group
locations <- as.data.frame(cbind(sample_info_ordered$ID, sample_info_ordered$location, sample_info_ordered$Fst_group ))
colnames(locations) <- c("ID", "location", "Fst_group")
# add plot order column

plot_order <- NULL
for (i in 1:nrow(locations)) {
  tmp <- if (locations[i,3] == "RBSAd") {
    print(7)
  } else if (locations[i,3] == "RBSAr") {
    print(6)
  } else if (locations[i,3] == "RBxRN") {
    print(5)
  } else if (locations[i,3] == "RNSA") {
    print(4)
  } else if (locations[i,3] == "RNxYB") {
    print(3)
  } else if (locations[i,3] == "YBxRB") {
    print(2)
  } else {
    print(1)
  }
  plot_order[i]<- tmp
}

locations$plot_order <- plot_order
locations <- locations[(order(locations$plot_order)),]
sample_info_ordered$plot_order <- plot_order

#order by plot order
sample_info_ordered <- sample_info_ordered[order(sample_info_ordered$plot_order),]

geno1 <- cbind(ind, geno)
geno2 <- geno1[order(match(geno1$ID, sample_info_ordered$ID)), ]
geno3 <- geno2[,-c(1)]
### DI: NOTE THAT geno3 STILL HAS A COLUMN CALLED "null" AT THE START, WHICH WERE ROW NUMBERS FROM 0 ON UP
# I THINK NEED TO REMOVE THIS, SO ADDING THIS LINE:
geno4 <- geno3[, -1]

start.pos <- min(pos$position)
end.pos <- max(pos$position)
num.inds <- nrow(geno4)
num_loc_cols <- length(locations[1,])

# NOT THE NUMBER OF COLUMNS IN geno. IN THIS CASE, num_loc_cols SHOULD BE ZERO, SO ADDING THIS NEXT LINE:

plot.group.colors <- c("#dddf5bff", "#c6924dff", "#2f5734ff", "#394ca2ff", "#b26ec4ff", "#bd5757ff","#bd7777ff") 

plot.groups <- unique(locations$Fst_group)
plot.groups2 <- unique(locations$plot_order)


SNP.positions_to_plot <- pos$position
num.inds <- nrow(geno4)

# Calculate allele freqs and sample sizes (use column Fst_group)
groups <- c("YBSA", "YBxRB", "RNxYB", "RNSA", "RBxRN", "RBSAr", "RBSAd")

combo <- cbind(locations, geno4)

temp.list <- getFreqsAndSampleSizes(combo, num_loc_cols, groups)
SNP.freqs <- temp.list$freqs
sample_size <- temp.list$sample_size
rm(temp.list)

group1 <- "YBSA"

#Get the genotypes for each bird at each snp
SNP.genotypes <- combo[,c(rep(TRUE, times=num_loc_cols))]

#narrow down to just the plot.groups you want
SNP.genotypes.subset <- SNP.genotypes[SNP.genotypes$Fst_group %in% plot.groups,]


alt.allele.hi.in.group1 <- which(SNP.freqs[rownames(SNP.freqs)==group1,] > 0.5) + num_loc_cols
SNP.genotypes.subset[,alt.allele.hi.in.group1] <- -1*SNP.genotypes.subset[,alt.allele.hi.in.group1] + 2

# plot everyone
pdf("genotype.x.individual/Pseudo1030rnrbSNPs.pdf", width=12, height=12)
chr.length <- max(pos$position)
genotype.colors <- c("#3f007d", "#807dba", "#dadaeb", "grey50")  # purple shades from colorbrewer

chr <- "Pseudo1030"
start.pos <- min(pos$position)
end.pos <- max(pos$position)

# choose which type to plot: 1=nucleotide positions; 2=spaced evenly; 3= both
plot.along.chromosome.type <- 2   

# plot along chromosome by nucleotide position:
if (plot.along.chromosome.type == 1 | plot.along.chromosome.type == 3) {
  plot(x=NULL, y=NULL, xlim=c(start.pos, end.pos+0.05*(end.pos-start.pos)), ylim=c(0, num.inds+1), main=paste0("Fixed FST SNPs on Pseudo1030 males"),xlab=paste0("Location along chromosome"), ylab="Individual")
  # generate my own plotting symbol (a rectangle)
  symbol.x <- c(-0.1, -0.1, 0.1, 0.1, -0.1)
  symbol.y <- c(1, -1, -1, 1, 1)
  plot.symbol <- cbind(symbol.x, symbol.y)
  symbol.size.x <- 5000  #width of box in nucleotides
  symbol.size.y <- 0.8  #height of box in units of individuals
  # cycle through individuals, graphing each type of genotype:
  for (i in 1:num.inds) {
    y <- (1*i)+num.inds+1  #reverses order of plot top-bottom  if (1*i) changed to (-1*i)
    lines(x = c(start.pos,end.pos), y = c(y,y), col = "grey")
    text(x = end.pos, y = y, labels=unlist(strsplit(as.character(sample_info_ordered$ID[i]), split='_', fixed=TRUE))[3], cex=0.3, pos=4)
    genotypes <- SNP.genotypes.subset[i, (num_loc_cols+1):length(SNP.genotypes.subset[1,])]
    hom.ref.locs <- SNP.positions_to_plot[genotypes == 0 & !is.na(genotypes)]
    my.symbols(hom.ref.locs, rep(y, times=length(hom.ref.locs)), plot.symbol, xsize=symbol.size.x, ysize=symbol.size.y, col="red")
    het.locs <- SNP.positions_to_plot[genotypes == 1 & !is.na(genotypes)]
    my.symbols(het.locs, rep(y, times=length(het.locs)), plot.symbol, xsize=symbol.size.x, ysize=symbol.size.y, col="orange")
    hom.alt.locs <- SNP.positions_to_plot[genotypes == 2 & !is.na(genotypes)]
    my.symbols(hom.alt.locs, rep(y, times=length(hom.alt.locs)), plot.symbol, xsize=symbol.size.x, ysize=symbol.size.y, col="yellow")
  }
}

# plot evenly spaced by SNP order along chromosome:
# make top part of fig (genotypes for individuals)
if (plot.along.chromosome.type == 2 | plot.along.chromosome.type == 3) {
  num.SNPs.to.plot <- length(SNP.positions_to_plot)
  
  plot(x=NULL, y=NULL, xlim=c(0.5-0.07*(num.SNPs.to.plot+0.5), 1.07*(num.SNPs.to.plot+0.5)), ylim=c(0.5-0.25*num.inds, num.inds+1), main=paste0("Fixed RNxRB FST SNPs on Pseudo1030 males"),xlab=paste0("Order along chromosome"), ylab="Individual")
  image.matrix <- t(as.matrix(SNP.genotypes.subset[, (num_loc_cols):length(SNP.genotypes.subset[1,])]))
  
  image.matrix[is.na(image.matrix)] <- 3
  
  group.color.box.loc.right <- 1.055*(num.SNPs.to.plot+0.5)
  group.color.box.loc.left <- 0.5-0.055*(num.SNPs.to.plot+0.5)
  box.width <- 0.005*num.SNPs.to.plot * 2
  group.color.box.x.right <- c(-box.width, -box.width, box.width, box.width, -box.width) + group.color.box.loc.right
  group.color.box.x.left <- c(-box.width, -box.width, box.width, box.width, -box.width) + group.color.box.loc.left
  group.color.box.y <- c(0.4, -0.4, -0.4, 0.4, 0.4)
  
  for (i in 1:num.inds) {
    y <- (-1*i)+num.inds+1  #reverses order of plot top-bottom
    name.text.bits <- unlist(strsplit(as.character(sample_info_ordered$ID[y]), split='_', fixed=TRUE))
    label.text <- name.text.bits[length(name.text.bits)]
    text(x = num.SNPs.to.plot+0.5, y = y, labels=label.text, cex=0.3, pos=4)
    text(x = 0.5, y = y, labels=label.text, cex=0.3, pos=2)
    polygon(group.color.box.x.right, y+group.color.box.y, border=NA, col=plot.group.colors[which(plot.groups2==locations$plot_order[y])])
    polygon(group.color.box.x.left, y+group.color.box.y, border=NA, col=plot.group.colors[which(plot.groups2==locations$plot_order[y])])
  }
  
  # generate my own plotting symbol (a rectangle)
  symbol.x <- c(-0.5, -0.5, 0.5, 0.5, -0.5)
  symbol.y <- c(0.4, -0.4, -0.4, 0.4, 0.4)
  # generate triangles for plotting heterozygotes
  triangle1.x <- c(-0.5, -0.5, 0.5, -0.5)
  triangle1.y <- c(0.4, -0.4, 0.4, 0.4)
  triangle2.x <- c(-0.5, 0.5, 0.5, -0.5)
  triangle2.y <- c(-0.4, -0.4, 0.4, -0.4)
  # cycle through individuals, graphing each type of genotype:
  for (i in 1:num.inds) {
    y <- (-1*i)+num.inds+1  #reverses order of plot top-bottom
    lines(x = c(0.5,num.SNPs.to.plot+0.5), y = c(y,y), col = "grey40")
    genotypes <- SNP.genotypes.subset[y, (num_loc_cols+1):length(SNP.genotypes.subset[1,])]
    
    hom.ref.locs <- which(genotypes == 0 & !is.na(genotypes))
    if (length(hom.ref.locs) > 0) {
      for (j in 1:length(hom.ref.locs)) {
        polygon(hom.ref.locs[j]+symbol.x, y+symbol.y, border=NA, col=genotype.colors[1])
      }
    }
    het.locs <- which(genotypes == 1 & !is.na(genotypes))
    if (length(het.locs) > 0) {
      for (j in 1:length(het.locs)) {
        #polygon(het.locs[j]+symbol.x, y+symbol.y, border=NA, col=genotype.colors[2])  # draws rectangle in hetero color
        polygon(het.locs[j]+triangle1.x, y+triangle1.y, border=NA, col=genotype.colors[1])  # draws triangle in hom ref color
        polygon(het.locs[j]+triangle2.x, y+triangle2.y, border=NA, col=genotype.colors[3])  # draws triangle in hom alt color
      }
    }
    hom.alt.locs <- which(genotypes == 2 & !is.na(genotypes))
    if (length(hom.alt.locs) > 0) {
      for (j in 1:length(hom.alt.locs)) {
        polygon(hom.alt.locs[j]+symbol.x, y+symbol.y, border=NA, col=genotype.colors[3])
      }
    }
  }
}

# make lower part of figure (indicating position along chromosome) #NOTE THAT CHROMOSOME LENGTH NOT QUITE THE TRUE LENGTH
chr.line.y <- 0.5-0.2*num.inds
top.hatch.line.y1 <- 0.5-0.005*num.inds
top.hatch.line.y2 <- 0.5-0.02*num.inds
low.hatch.line.y1 <- 0.5-0.18*num.inds
low.hatch.line.y2 <- 0.5-0.2*num.inds
lines(x = c(0.5,num.SNPs.to.plot+0.5), y = c(chr.line.y,chr.line.y), lwd=4, col = "black") #draws chromosome line
text(x=(0.5+(num.SNPs.to.plot+0.5)/2), y=chr.line.y-0.05*num.inds, paste0("Location along chromosome ",chr, sep=""))
text(x=0.5, y=chr.line.y-0.025*num.inds, start.pos)
text(x=num.SNPs.to.plot+0.5, y=chr.line.y-0.025*num.inds, end.pos)
chr.plot.ratio <- num.SNPs.to.plot/(end.pos-start.pos)
for (i in 1:length(SNP.positions_to_plot)) {
  lines(x=c(i,i), y=c(top.hatch.line.y1, top.hatch.line.y2), lwd=0.5, col="grey20")
  lines(x=c(i, 1+chr.plot.ratio*(SNP.positions_to_plot[i]-start.pos)), y=c(top.hatch.line.y2, low.hatch.line.y1), lwd=0.5, col="grey20")
  lines(x=c(1+chr.plot.ratio*(SNP.positions_to_plot[i]-start.pos), 1+chr.plot.ratio*(SNP.positions_to_plot[i]-start.pos)), y=c(low.hatch.line.y1, low.hatch.line.y2), lwd=0.5, col="grey20")
}
dev.off()


## rnrb only
plot.group.colors <- c("#394ca2ff", "#b26ec4ff", "#bd5757ff","#bd7777ff") 

plot.groups <- c("RNSA", "RBxRN", "RBSAr", "RBSAd")
plot.groups2 <- c(4:7)


SNP.positions_to_plot <- pos$position
num.inds <- nrow(geno4)

# Calculate allele freqs and sample sizes (use column Fst_group)
groups <- c("RNSA", "RBxRN", "RBSAr", "RBSAd")

combo <- cbind(locations, geno4)

temp.list <- getFreqsAndSampleSizes(combo, num_loc_cols, groups)
SNP.freqs <- temp.list$freqs
sample_size <- temp.list$sample_size
rm(temp.list)

group1 <- "RNSA"

#Get the genotypes for each bird at each snp
SNP.genotypes <- combo[,c(rep(TRUE, times=num_loc_cols))]

#narrow down to just the plot.groups you want
SNP.genotypes.subset <- SNP.genotypes[SNP.genotypes$Fst_group %in% plot.groups,]
num.inds <- nrow(SNP.genotypes.subset)

alt.allele.hi.in.group1 <- which(SNP.freqs[rownames(SNP.freqs)==group1,] > 0.5) + num_loc_cols
SNP.genotypes.subset[,alt.allele.hi.in.group1] <- -1*SNP.genotypes.subset[,alt.allele.hi.in.group1] + 2


#pdf("genotype.x.individual/Pseudo1030rnrbSNPs_rnrbonly.pdf", width=12, height=12)
chr.length <- max(pos$position)
genotype.colors <- c("#3f007d", "#807dba", "#dadaeb", "grey50")  # purple shades from colorbrewer

chr <- "Pseudo1030"
start.pos <- min(pos$position)
end.pos <- max(pos$position)

# choose which type to plot: 1=nucleotide positions; 2=spaced evenly; 3= both
plot.along.chromosome.type <- 2   

# plot along chromosome by nucleotide position:
if (plot.along.chromosome.type == 1 | plot.along.chromosome.type == 3) {
  plot(x=NULL, y=NULL, xlim=c(start.pos, end.pos+0.05*(end.pos-start.pos)), ylim=c(0, num.inds+1), main=paste0("Fixed RNxRB FST SNPs on Pseudo1030 males"),xlab=paste0("Location along chromosome"), ylab="Individual")
  # generate my own plotting symbol (a rectangle)
  symbol.x <- c(-0.1, -0.1, 0.1, 0.1, -0.1)
  symbol.y <- c(1, -1, -1, 1, 1)
  plot.symbol <- cbind(symbol.x, symbol.y)
  symbol.size.x <- 5000  #width of box in nucleotides
  symbol.size.y <- 0.8  #height of box in units of individuals
  # cycle through individuals, graphing each type of genotype:
  for (i in 1:num.inds) {
    y <- (1*i)+num.inds+1  #reverses order of plot top-bottom  if (1*i) changed to (-1*i)
    lines(x = c(start.pos,end.pos), y = c(y,y), col = "grey")
    text(x = end.pos, y = y, labels=unlist(strsplit(as.character(SNP.genotypes.subset$ID[i]), split='_', fixed=TRUE))[3], cex=0.3, pos=4)
    genotypes <- SNP.genotypes.subset[i, (num_loc_cols+1):length(SNP.genotypes.subset[1,])]
    hom.ref.locs <- SNP.positions_to_plot[genotypes == 0 & !is.na(genotypes)]
    my.symbols(hom.ref.locs, rep(y, times=length(hom.ref.locs)), plot.symbol, xsize=symbol.size.x, ysize=symbol.size.y, col="red")
    het.locs <- SNP.positions_to_plot[genotypes == 1 & !is.na(genotypes)]
    my.symbols(het.locs, rep(y, times=length(het.locs)), plot.symbol, xsize=symbol.size.x, ysize=symbol.size.y, col="orange")
    hom.alt.locs <- SNP.positions_to_plot[genotypes == 2 & !is.na(genotypes)]
    my.symbols(hom.alt.locs, rep(y, times=length(hom.alt.locs)), plot.symbol, xsize=symbol.size.x, ysize=symbol.size.y, col="yellow")
  }
}

# plot evenly spaced by SNP order along chromosome:
# make top part of fig (genotypes for individuals)
if (plot.along.chromosome.type == 2 | plot.along.chromosome.type == 3) {
  num.SNPs.to.plot <- length(SNP.positions_to_plot)
  
  plot(x=NULL, y=NULL, xlim=c(0.5-0.07*(num.SNPs.to.plot+0.5), 1.07*(num.SNPs.to.plot+0.5)), ylim=c(0.5-0.25*num.inds, num.inds+1), main=paste0("Fixed RNxRB FST SNPs on Pseudo1030 males"),xlab=paste0("Order along chromosome"), ylab="Individual")
  image.matrix <- t(as.matrix(SNP.genotypes.subset[, (num_loc_cols):length(SNP.genotypes.subset[1,])]))
  
  image.matrix[is.na(image.matrix)] <- 3
  
  group.color.box.loc.right <- 1.055*(num.SNPs.to.plot+0.5)
  group.color.box.loc.left <- 0.5-0.055*(num.SNPs.to.plot+0.5)
  box.width <- 0.005*num.SNPs.to.plot * 2
  group.color.box.x.right <- c(-box.width, -box.width, box.width, box.width, -box.width) + group.color.box.loc.right
  group.color.box.x.left <- c(-box.width, -box.width, box.width, box.width, -box.width) + group.color.box.loc.left
  group.color.box.y <- c(0.4, -0.4, -0.4, 0.4, 0.4)
  
  for (i in 1:num.inds) {
    y <- (-1*i)+num.inds+1  #reverses order of plot top-bottom
    name.text.bits <- unlist(strsplit(as.character(SNP.genotypes.subset$ID[y]), split='_', fixed=TRUE))
    label.text <- name.text.bits[length(name.text.bits)]
    text(x = num.SNPs.to.plot+0.5, y = y, labels=label.text, cex=0.4, pos=4)
    text(x = 0.5, y = y, labels=label.text, cex=0.4, pos=2)
    polygon(group.color.box.x.right, y+group.color.box.y, border=NA, col=plot.group.colors[which(plot.groups2==SNP.genotypes.subset$plot_order[y])])
    polygon(group.color.box.x.left, y+group.color.box.y, border=NA, col=plot.group.colors[which(plot.groups2==SNP.genotypes.subset$plot_order[y])])
  }
  
  # generate my own plotting symbol (a rectangle)
  symbol.x <- c(-0.5, -0.5, 0.5, 0.5, -0.5)
  symbol.y <- c(0.4, -0.4, -0.4, 0.4, 0.4)
  # generate triangles for plotting heterozygotes
  triangle1.x <- c(-0.5, -0.5, 0.5, -0.5)
  triangle1.y <- c(0.4, -0.4, 0.4, 0.4)
  triangle2.x <- c(-0.5, 0.5, 0.5, -0.5)
  triangle2.y <- c(-0.4, -0.4, 0.4, -0.4)
  # cycle through individuals, graphing each type of genotype:
  for (i in 1:num.inds) {
    y <- (-1*i)+num.inds+1  #reverses order of plot top-bottom
    lines(x = c(0.5,num.SNPs.to.plot+0.5), y = c(y,y), col = "grey40")
    genotypes <- SNP.genotypes.subset[y, (num_loc_cols+1):length(SNP.genotypes.subset[1,])]
    
    hom.ref.locs <- which(genotypes == 0 & !is.na(genotypes))
    if (length(hom.ref.locs) > 0) {
      for (j in 1:length(hom.ref.locs)) {
        polygon(hom.ref.locs[j]+symbol.x, y+symbol.y, border=NA, col=genotype.colors[1])
      }
    }
    het.locs <- which(genotypes == 1 & !is.na(genotypes))
    if (length(het.locs) > 0) {
      for (j in 1:length(het.locs)) {
        #polygon(het.locs[j]+symbol.x, y+symbol.y, border=NA, col=genotype.colors[2])  # draws rectangle in hetero color
        polygon(het.locs[j]+triangle1.x, y+triangle1.y, border=NA, col=genotype.colors[1])  # draws triangle in hom ref color
        polygon(het.locs[j]+triangle2.x, y+triangle2.y, border=NA, col=genotype.colors[3])  # draws triangle in hom alt color
      }
    }
    hom.alt.locs <- which(genotypes == 2 & !is.na(genotypes))
    if (length(hom.alt.locs) > 0) {
      for (j in 1:length(hom.alt.locs)) {
        polygon(hom.alt.locs[j]+symbol.x, y+symbol.y, border=NA, col=genotype.colors[3])
      }
    }
  }
}

# make lower part of figure (indicating position along chromosome) #NOTE THAT CHROMOSOME LENGTH NOT QUITE THE TRUE LENGTH
chr.line.y <- 0.5-0.2*num.inds
top.hatch.line.y1 <- 0.5-0.005*num.inds
top.hatch.line.y2 <- 0.5-0.02*num.inds
low.hatch.line.y1 <- 0.5-0.18*num.inds
low.hatch.line.y2 <- 0.5-0.2*num.inds
lines(x = c(0.5,num.SNPs.to.plot+0.5), y = c(chr.line.y,chr.line.y), lwd=4, col = "black") #draws chromosome line
text(x=(0.5+(num.SNPs.to.plot+0.5)/2), y=chr.line.y-0.05*num.inds, paste0("Location along chromosome ",chr, sep=""))
text(x=0.5, y=chr.line.y-0.025*num.inds, start.pos)
text(x=num.SNPs.to.plot+0.5, y=chr.line.y-0.025*num.inds, end.pos)
chr.plot.ratio <- num.SNPs.to.plot/(end.pos-start.pos)
for (i in 1:length(SNP.positions_to_plot)) {
  lines(x=c(i,i), y=c(top.hatch.line.y1, top.hatch.line.y2), lwd=0.5, col="grey20")
  lines(x=c(i, 1+chr.plot.ratio*(SNP.positions_to_plot[i]-start.pos)), y=c(top.hatch.line.y2, low.hatch.line.y1), lwd=0.5, col="grey20")
  lines(x=c(1+chr.plot.ratio*(SNP.positions_to_plot[i]-start.pos), 1+chr.plot.ratio*(SNP.positions_to_plot[i]-start.pos)), y=c(low.hatch.line.y1, low.hatch.line.y2), lwd=0.5, col="grey20")
}
#dev.off()



