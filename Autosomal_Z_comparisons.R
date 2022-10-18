# script to make plotGenomeFstDxyPi plots 
# adapted by Libby Natola from Darren Irwin's scripts from 2016 paper
# Started 12 August 2022

setwd("~/Documents/UBC/Bioinformatics/wgs")
#setwd("./")

# Load functions
source("genomics_R_functions_V3.R")
#source("genomics_R_functions_V5.R")

# install.packages("vroom")
library(vroom)   # for fread function for fast reading of data files

# choose the chromosomes to analyze in this run
#chromosomes.to.analyze <- c("wgs.genotypes.allSites.PseudoWWNC01000346.1_Melanerpes_aurifrons_Melanerpes_aurifrons_OMNH24340_contig_346_whole_genome_shotgun_sequence.filtered.missing80_mindepth3.012")
# to process other chromosomes, choose from among these:
#chromosomes.to.analyze <- c(88, 103, 106, 110, 111, 112, 114, 115, 120, 1073, 1076, 1080, 1082, 1095, 1097, 1099, 1100, 1144)
chromosomes.to.analyze <- c(1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 45, 46, 49, 51, 52, 54, 55, 60, 63, 65, 67, 70, 74, 75, 77, 78, 86, 87, 88, 92, 96, 103, 106, 110, 111, 112, 114, 115, 120, 125, 127, 131, 132, 134, 136, 139, 143, 145, 150, 153, 157, 159, 161, 163, 164, 172, 177, 181, 182, 185, 186, 196, 200,  207, 221, 239, 241, 257, 269, 298, 1003, 1030, 1044, 1048, 1059, 1073, 1076, 1080, 1082, 1095, 1097, 1099, 1100, 1144)


# choose path and filename for the 012NA files
base.name <- paste0("allSites_012NA_renamed/wgs.genotypes.allSites.",chromosomes.to.analyze,".filtered.missing80_mindepth3.wgs")
tag.name <- ".wgs_w50000."   # choose a tag name for this analysis
# indicate name of metadata file, a text file with these column headings:
# ID	location	group	Fst_group	plot_order
metadata.file <- "wgs.Fst_groups.txt"
# load metadata
locations <- read.table(paste0(metadata.file), sep = "\t", header=TRUE)
num_loc_cols <- length(locations[1,])
num.individuals <- 78  # specify number of individuals in file (good for error checking)
# specify window size (number of bp with info) and step size
window_size <- 50000
step_size <- window_size  # could change if wanting overlapping windows
# specify groups for calculation of statistics (these are in "Fst_group" column in metadata file)
groups <- c("YBSA", "YBxRB", "RNxYB", "RNSA", "RBxRN", "RBSA")
group.colors <- c("#dddf5bff", "#c6924dff", "#2f5734ff", "#394ca2ff", "#b26ec4ff", "#bd5757ff") 
group_count <- length(groups)
# specify groups for plotting, and their colors
groups.to.plot.WC84_Fst <- c("YBSA_RNSA",
                             "RNSA_RBSA",
                             "YBSA_RBSA"
)
group.colors.WC84_Fst <- c("#2f5734",
                           "#bda2c4",
                           "#c69f6b"
)
groups.to.plot.Dxy <- groups.to.plot.WC84_Fst   # or specify differences if necessary
group.colors.Dxy <- group.colors.WC84_Fst  # or specify differences if necessary
groups.to.plot.pi <- c("YBSA", "RNSA", "RBSA")
group.colors.pi <- c("#e6e045", "#394ca2", "#bd5757")
groups.to.plot.common.names <- c("Yellow-bellied and Red-naped", "Red-naped and Red-breasted", "Yellow-bellied and Red-breasted")

#region.text <- paste0("Chr",chromosomes.to.analyze,"_whole")
region.text <- paste0("Chr",chromosomes.to.analyze,"_whole") 
to_load <- paste0(base.name,tag.name,region.text,"_window",window_size,"_WindowStats_from_R.R")
lapply(to_load,load,.GlobalEnv)

# for(i in 1:length(chromosomes.to.analyze)){
#   chr <- chromosomes.to.analyze[i] # choose chromosome for the loop
#   region.text <- paste0("Chr",chr,"_whole") 
#   
#   #load(paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R")) # load the SiteStats
#   #print("Loaded saved summary stats")
#   to_load <- paste0(base.name[i],tag.name,region.text,"_window",window_size,"_WindowStats_from_R.R")
#   lapply(to_load_fun2,load,.GlobalEnv) # load the WindowStats JLD.rolling mean data
# }

# specify sex vs autosomes
#chromosomes.to.combine <- c('88', '103', '106', '110', '111', '112', '114', '115', '120')
chromosomes.to.combine <- c(1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 45, 46, 49, 51, 52, 54, 55, 60, 63, 65, 67, 70, 74, 75, 77, 78, 86, 87, 88, 92, 96, 103, 106, 110, 111, 112, 114, 115, 120, 125, 127, 131, 132, 134, 136, 139, 143, 145, 150, 153, 157, 159, 161, 163, 164, 172, 177, 181, 182, 185, 186, 196, 200,  207, 221, 239, 241, 257, 269, 298)

base.name <- paste0("allSites_012NA_renamed/wgs.genotypes.allSites.",chromosomes.to.combine,".filtered.missing80_mindepth3.wgs")
autosome.genome.rolling.stats <- compileWindowedStats(base.name, tag.name, chromosomes.to.combine, window_size)

# get number of autosome windows:
length(autosome.genome.rolling.stats$positions[1,])
#228
#[1] 16840

# compile Z-chromosome info on WC84_Fst, Dxy, pi (this part copied from above):
#chromosomes.to.combine <- c('1073', '1076', '1080', '1082', '1095', '1097', '1099', '1100', '1144')
chromosomes.to.combine <- c(1003, 1030, 1044, 1048, 1059, 1073, 1076, 1080, 1082, 1095, 1097, 1099, 1100, 1144)
base.name <- paste0("allSites_012NA_renamed/wgs.genotypes.allSites.",chromosomes.to.combine,".filtered.missing80_mindepth3.wgs")
Z.rolling.stats <- compileWindowedStats(base.name, tag.name, chromosomes.to.combine, window_size)

# get number of Z windows:
length(Z.rolling.stats$positions[1,])
#[1] 1710


# group1 <- "YBSA"
# group2 <- "RBSA"
# 
# 
# # MeanPi_Dxy autosome vs. Z plot ----
# # uses compiled data produced by function "compileWindowedStats";
# # for both autsome and Z
# plotMeanPi_Dxy.autosome_Z <- function(group1, group2,
#                                       autosome.genome.rolling.stats, Z.rolling.stats) {
#   groups.for.graph <- paste0(group1, "_", group2)
#   row.choice.1 <- which(rownames(autosome.genome.rolling.stats$pi) == group1)
#   row.choice.2 <- which(rownames(autosome.genome.rolling.stats$pi) == group2)
#   autosome.pi_1 <- autosome.genome.rolling.stats$pi[row.choice.1,]
#   autosome.pi_2 <- autosome.genome.rolling.stats$pi[row.choice.2,]
#   autosome.mean_pi <- (autosome.pi_1 + autosome.pi_2) / 2
#   Dxy.row.choice <- which(rownames(autosome.genome.rolling.stats$Dxy) == groups.for.graph)
#   autosome.Dxy <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice,])
#   # Z calcs:
#   row.choice.1 <- which(rownames(Z.rolling.stats$pi) == group1)
#   row.choice.2 <- which(rownames(Z.rolling.stats$pi) == group2)
#   Z.pi_1 <- Z.rolling.stats$pi[row.choice.1,]
#   Z.pi_2 <- Z.rolling.stats$pi[row.choice.2,]
#   Z.mean_pi <- (Z.pi_1 + Z.pi_2) / 2
#   Dxy.row.choice <- which(rownames(Z.rolling.stats$Dxy) == groups.for.graph)
#   Z.Dxy <- as.vector(Z.rolling.stats$Dxy[Dxy.row.choice,])
#   # make the figure:
#   quartz(title=paste0("Autosomes vs. Z: Scatterplot of windowed mean pi vs. Dxy between ", group1,"_", group2, sep=""), width=6, height=6)
#   par(oma=c(3,3,1,1))  # set outer margins
#   zones <- matrix(c(4,0,0,2,0,0,1,3,5), ncol=3, byrow=TRUE)  # numbers in matrix give order of plotting
#   layout(zones, widths=c(3/5,1/5,1/5), heights=c(1/5,1/5,3/5))
#   xlimits <- c(0, round(max(c(autosome.Dxy,Z.Dxy))*1.1, digits=3))  # rounds to 3 decimal places
#   ylimits <- c(0, round(max(c(autosome.mean_pi,Z.mean_pi))*1.1, digits=3))
#   limits <- c(min(c(xlimits,ylimits)), max(xlimits,ylimits)) # this to make axes have same limits
#   xhist <- hist(autosome.Dxy, breaks=seq(limits[1], limits[2], by=0.0005), plot=FALSE)
#   chrZ.xhist <- hist(Z.Dxy, breaks=seq(limits[1], limits[2], by=0.0005), plot=FALSE)
#   yhist <- hist(autosome.mean_pi, breaks=seq(limits[1], limits[2], by=0.0005), plot=FALSE)
#   chrZ.yhist <- hist(Z.mean_pi, breaks=seq(limits[1], limits[2], by=0.0005), plot=FALSE)
#   top <- max(c(xhist$counts, yhist$counts))
#   chrZ.top <- max(c(chrZ.xhist$counts, chrZ.yhist$counts))
#   par(mar=c(3,3,1,1))  # specifies number of lines around plot (bottom, left, top right)
#   plot(autosome.Dxy,autosome.mean_pi, col="grey70", xlim=limits, 
#        ylim=limits, pch=16, cex=0.3, asp=1)
#   lines(c(0,1), c(0,1))
#   points(Z.Dxy, Z.mean_pi, col='blue', pch=16, cex=0.5, asp=1)
#   par(mar=c(0,3,0,1))
#   barplot(xhist$counts, axes=FALSE, ylim=c(0, top*1.1), col="grey70", space=0)
#   axis(side=2, at=c(0,200))
#   par(mar=c(3,0,1,0.5))
#   barplot(yhist$counts, axes=FALSE, xlim=c(0, top*1.1), col="grey70", space=0, horiz=TRUE)
#   axis(side=1, at=c(0,200))
#   par(mar=c(0,3,1,1))
#   barplot(chrZ.xhist$counts, axes=FALSE, ylim=c(0, top*1.1), space=0, col="blue")
#   axis(side=2, at=c(0,200))
#   par(mar=c(3,0,1,1))
#   barplot(chrZ.yhist$counts, axes=FALSE, xlim=c(0, top*1.1), space=0, col="blue", horiz=TRUE)
#   axis(side=1, at=c(0,200))
#   par(oma=c(3,3,0,0))
#   label.size <- 1.2
#   mtext(expression("Between-group distance (" * italic(pi)[B] * ")"), side=1, line=0.5, outer=TRUE, 
#         at=1.5/5)
#   mtext(expression("Mean within-group variation (" * italic(pi)[W] * ")"), side=2, line=0, outer=TRUE, 
#         at=1.5/5)
#   mtext("Windows", side=1, line=-0.5, outer=TRUE, cex=0.8, 
#         at=3.15/5)
#   mtext("Windows", side=1, line=-0.5, outer=TRUE, cex=0.8, 
#         at=4.1/5)
#   mtext("Windows", side=2, line=-0.5, outer=TRUE, cex=0.8, 
#         at=3.15/5)
#   mtext("Windows", side=2, line=-0.5, outer=TRUE, cex=0.8, 
#         at=4.1/5)
#   # t-tests:
#   test.mean_pi <- t.test(autosome.mean_pi, Z.mean_pi)
#   test.Dxy <- t.test(autosome.Dxy, Z.Dxy)
#   return(list(test.mean_pi=test.mean_pi, test.Dxy=test.Dxy))
# }
# 
# plotMeanPi_Dxy.autosome_Z(group1, group2, autosome.genome.rolling.stats, Z.rolling.stats)
# 


# Make a plot of windowed Fst vs. Dxy for autosomes and Z
plotFst_Dxy_auto_z <- function(comp, autosome.genome.rolling.stats, Z.rolling.stats, comp_common, color) {
  pdf(file=(paste("Fst_PiB_autosomal_Z_", comp, ".pdf", sep="")), width = 4, height = 4)
  rowchoice.autosome <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == comp)
  rowchoice.Z <- which(rownames(Z.rolling.stats$WC84_Fst) == comp)
  plot(autosome.genome.rolling.stats$WC84_Fst[rowchoice.autosome,], autosome.genome.rolling.stats$Dxy[rowchoice.autosome,], 
       pch=16, cex=0.3, col = alpha('gray45', 0.6), xlab=paste0("Windowed FST for ", comp_common, sep = ""),
       ylab=paste0("Windowed piB for ", comp_common, sep=""), xlim = c(0, 1), ylim = c(0, 0.035))
  points(Z.rolling.stats$WC84_Fst[rowchoice.Z,], Z.rolling.stats$Dxy[rowchoice.Z,], 
       pch=16, cex=0.3, col = alpha(color, 0.6), xlim = c(0, 1), ylim = c(0, 0.035))
  dev.off()
}

comp <- groups.to.plot.WC84_Fst[1]
comp_common <- groups.to.plot.common.names[1]
color <- group.colors.WC84_Fst[1]
plotFst_Dxy_auto_z(comp, autosome.genome.rolling.stats, Z.rolling.stats, comp_common, color)


comp <- groups.to.plot.WC84_Fst[2]
comp_common <- groups.to.plot.common.names[2]
color <- group.colors.WC84_Fst[2]
plotFst_Dxy_auto_z(comp, autosome.genome.rolling.stats, Z.rolling.stats, comp_common, color)


comp <- groups.to.plot.WC84_Fst[3]
comp_common <- groups.to.plot.common.names[3]
color <- group.colors.WC84_Fst[3]
plotFst_Dxy_auto_z(comp, autosome.genome.rolling.stats, Z.rolling.stats, comp_common, color)

# # Fst_Fst plot
# ### but add z 
# # Make a bivariate plot of WC84_Fst using for two population pair comparisons:
# # return the correlation, according to cor.method
# plotFst_ZFst <- function(comparison, cor.method, 
#                         autosome.genome.rolling.stats, Z.rolling.stats) {
#   quartz(title=paste0("Bivariate plot of windowed WC84_Fst, based on all autosomes vs zchroms, ", comparison, sep=""), width=6, height=6)
#   rowchoice.autosome <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == comparison)
#   rowchoice.Z <- which(rownames(Z.rolling.stats$WC84_Fst) == comparison)
#   plot(autosome.genome.rolling.stats$WC84_Fst[rowchoice.autosome, ], Z.rolling.stats$WC84_Fst[rowchoice.Z, ], cex=0.5)
#   # lines(c(0,1), c(0,1))
#   test <- cor.test(autosome.genome.rolling.stats$WC84_Fst[rowchoice.autosome, ], Z.rolling.stats$WC84_Fst[rowchoice.Z, ], method=cor.method)
#   return(test)
# }
# 
# comparison <- 'RNSA_RBSA'
# cor.method <- 'pearson'
# plotFst_ZFst(comparison, cor.method, autosome.genome.rolling.stats, Z.rolling.stats)



# statistical test of WC84_Fst diff between autosomes and Z:
comparison <- groups.to.plot.WC84_Fst[1]
rowchoice.autosome <-which(rownames(autosome.genome.rolling.stats$WC84_Fst) == comparison)
autosome_data <- autosome.genome.rolling.stats$WC84_Fst[rowchoice.autosome, ]
rowchoice.Z <-which(rownames(Z.rolling.stats$WC84_Fst) == comparison)
Z_data <- Z.rolling.stats$WC84_Fst[rowchoice.Z, ]
print(comparison)
t.test(autosome_data, Z_data)

#[1] "YBSA_RNSA"
# Welch Two Sample t-test
# 
# data:  autosome_data and Z_data
# t = -62.818, df = 1795.5, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.3157181 -0.2966004
# sample estimates:
#   mean of x mean of y 
# 0.2671136 0.5732728 

comparison <- groups.to.plot.WC84_Fst[2]
rowchoice.autosome <-which(rownames(autosome.genome.rolling.stats$WC84_Fst) == comparison)
autosome_data <- autosome.genome.rolling.stats$WC84_Fst[rowchoice.autosome, ]
rowchoice.Z <-which(rownames(Z.rolling.stats$WC84_Fst) == comparison)
Z_data <- Z.rolling.stats$WC84_Fst[rowchoice.Z, ]
print(comparison)
t.test(autosome_data, Z_data)
# [1] "RNSA_RBSA"
# 
# Welch Two Sample t-test
# 
# data:  autosome_data and Z_data
# t = -64.229, df = 1758.3, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.3239673 -0.3047680
# sample estimates:
#   mean of x  mean of y 
# 0.09782254 0.41219017 

comparison <- groups.to.plot.WC84_Fst[3]
rowchoice.autosome <-which(rownames(autosome.genome.rolling.stats$WC84_Fst) == comparison)
autosome_data <- autosome.genome.rolling.stats$WC84_Fst[rowchoice.autosome, ]
rowchoice.Z <-which(rownames(Z.rolling.stats$WC84_Fst) == comparison)
Z_data <- Z.rolling.stats$WC84_Fst[rowchoice.Z, ]
print(comparison)
t.test(autosome_data, Z_data)
# [1] "YBSA_RBSA"
# 
# Welch Two Sample t-test
# 
# data:  autosome_data and Z_data
# t = -62.142, df = 1776.1, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.3481830 -0.3268769
# sample estimates:
#   mean of x mean of y 
# 0.2088815 0.5464115 

# statistical test of WC84_Fst diff between autosomes and Z:
comparison <- groups.to.plot.WC84_Fst[1]
rowchoice.autosome.dxy <- which(rownames(autosome.genome.rolling.stats$Dxy) == comparison)
autosome_data.dxy <- autosome.genome.rolling.stats$Dxy[rowchoice.autosome.dxy, ]
rowchoice.Z.dxy <-which(rownames(Z.rolling.stats$Dxy) == comparison)
Z_data.dxy <- Z.rolling.stats$Dxy[rowchoice.Z.dxy, ]
print(comparison)
t.test(autosome_data.dxy, Z_data.dxy)

# [1] "YBSA_RNSA"
# 
# Welch Two Sample t-test
# 
# data:  autosome_data.dxy and Z_data.dxy
# t = 16.447, df = 1893.9, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.001093282 0.001389313
# sample estimates:
#   mean of x   mean of y 
# 0.004826610 0.003585313 


comparison <- groups.to.plot.WC84_Fst[2]
rowchoice.autosome.dxy <-which(rownames(autosome.genome.rolling.stats$Dxy) == comparison)
autosome_data.dxy <- autosome.genome.rolling.stats$Dxy[rowchoice.autosome.dxy, ]
rowchoice.Z.dxy <-which(rownames(Z.rolling.stats$Dxy) == comparison)
Z_data.dxy <- Z.rolling.stats$Dxy[rowchoice.Z.dxy, ]
print(comparison)
t.test(autosome_data.dxy, Z_data.dxy)

# [1] "RNSA_RBSA"
# Welch Two Sample t-test
# 
# data:  autosome_data.dxy and Z_data.dxy
# t = 16.4, df = 1874.8, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.000945106 0.001201859
# sample estimates:
#   mean of x   mean of y 
# 0.003206930 0.002133447 


comparison <- groups.to.plot.WC84_Fst[3]
rowchoice.autosome.dxy <-which(rownames(autosome.genome.rolling.stats$Dxy) == comparison)
autosome_data.dxy <- autosome.genome.rolling.stats$Dxy[rowchoice.autosome.dxy, ]
rowchoice.Z.dxy <-which(rownames(Z.rolling.stats$Dxy) == comparison)
Z_data.dxy <- Z.rolling.stats$Dxy[rowchoice.Z.dxy, ]
print(comparison)
t.test(autosome_data.dxy, Z_data.dxy)

# [1] "YBSA_RBSA"
# Welch Two Sample t-test
# 
# data:  autosome_data.dxy and Z_data.dxy
# t = 14.219, df = 1855.4, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.001007245 0.001329568
# sample estimates:
#   mean of x   mean of y 
# 0.004779322 0.003610916 


species <- groups.to.plot.pi[1]
rowchoice.autosome.pi <-which(rownames(autosome.genome.rolling.stats$pi) == species)
autosome_data.pi <- autosome.genome.rolling.stats$pi[rowchoice.autosome.pi, ]
rowchoice.Z.pi <-which(rownames(Z.rolling.stats$pi) == species)
Z_data.pi <- Z.rolling.stats$pi[rowchoice.Z.pi, ]
print(species)
t.test(autosome_data.pi, Z_data.pi)

# [1] "YBSA"
# Welch Two Sample t-test
# 
# data:  autosome_data.pi and Z_data.pi
# t = 21.006, df = 1835.8, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.001815459 0.002189367
# sample estimates:
#   mean of x   mean of y 
# 0.004583915 0.002581502 

species <- groups.to.plot.pi[2]
rowchoice.autosome.pi <-which(rownames(autosome.genome.rolling.stats$pi) == species)
autosome_data.pi <- autosome.genome.rolling.stats$pi[rowchoice.autosome.pi, ]
rowchoice.Z.pi <-which(rownames(Z.rolling.stats$pi) == species)
Z_data.pi <- Z.rolling.stats$pi[rowchoice.Z.pi, ]
print(species)
t.test(autosome_data.pi, Z_data.pi)

# [1] "RNSA"
# Welch Two Sample t-test
# 
# data:  autosome_data.pi and Z_data.pi
# t = 24.803, df = 1955.5, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.001232276 0.001443883
# sample estimates:
#   mean of x   mean of y 
# 0.002696945 0.001358866 


species <- groups.to.plot.pi[3]
rowchoice.autosome.pi <-which(rownames(autosome.genome.rolling.stats$pi) == species)
autosome_data.pi <- autosome.genome.rolling.stats$pi[rowchoice.autosome.pi, ]
rowchoice.Z.pi <-which(rownames(Z.rolling.stats$pi) == species)
Z_data.pi <- Z.rolling.stats$pi[rowchoice.Z.pi, ]
print(species)
t.test(autosome_data.pi, Z_data.pi)

# [1] "RBSA"
# Welch Two Sample t-test
# 
# data:  autosome_data.pi and Z_data.pi
# t = 18.461, df = 1820.1, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.001283569 0.001588715
# sample estimates:
#   mean of x   mean of y 
# 0.003157573 0.001721431 

### and then take the autosomal loci, compare fst, dxy, and pi like in the 2018 paper to describe source of selection on these loci: loci causing RI, positive and/or BG selection, positive selection and selective sweeps
# make scatterplots of Fst vs. Dxy for three group comparisons.
# each of the "comp_" objects should be have two group names and a color name, e.g.:
# comp1 <- c("troch", "vir", "green3") 
# returns the statistical correlation tests, using the cor.method as specified
plots3Fst_Dxy <- function(comp1, comp2, comp3,
                          autosome.genome.rolling.stats, cor.method) 
  {
  pdf(file="Autosome_Fst_PiB.pdf", width=8, height=3)
  par(oma=c(1,3,1,3))  # set outer margins
  par(mfrow=c(1,3))
  label.size <- 1.5
  Fst.row.choice1 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == paste(comp1[1], comp1[2], sep="_"))
  Fst.row.choice2 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == paste(comp2[1], comp2[2], sep="_"))
  Fst.row.choice3 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == paste(comp3[1], comp3[2], sep="_"))
  Fst.vector1 <- autosome.genome.rolling.stats$WC84_Fst[Fst.row.choice1,]
  Fst.vector2 <- autosome.genome.rolling.stats$WC84_Fst[Fst.row.choice2,]
  Fst.vector3 <- autosome.genome.rolling.stats$WC84_Fst[Fst.row.choice3,]
  Dxy.row.choice1 <- which(rownames(autosome.genome.rolling.stats$Dxy) == paste(comp1[1], comp1[2], sep="_"))
  Dxy.row.choice2 <- which(rownames(autosome.genome.rolling.stats$Dxy) == paste(comp2[1], comp2[2], sep="_"))
  Dxy.row.choice3 <- which(rownames(autosome.genome.rolling.stats$Dxy) == paste(comp3[1], comp3[2], sep="_"))
  Dxy.vector1 <- autosome.genome.rolling.stats$Dxy[Dxy.row.choice1,]
  Dxy.vector2 <- autosome.genome.rolling.stats$Dxy[Dxy.row.choice2,]
  Dxy.vector3 <- autosome.genome.rolling.stats$Dxy[Dxy.row.choice3,]
  high.Dxy <- round(max(c(Dxy.vector1, Dxy.vector2, Dxy.vector3)), 3) + 0.001
  # first plot:
  plot(x=Fst.vector1, y=Dxy.vector1, xlim=c(0,1), xaxp=c(0,1,2), ylim=c(0,high.Dxy), yaxp=c(0,high.Dxy,2), 
       pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic('F')[ST] * ", " * phantom(italic(.(comp1[1])) * " to " * italic(.(comp1[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic('F')[ST] * ", ") * italic(.(comp1[1])) * " to " * italic(.(comp1[2]))), line=2.5, cex.lab=label.size, col.lab=comp1[3])
  title(ylab=bquote(italic(pi)[B] * ", " * phantom(italic(.(comp1[1])) * " to " * italic(.(comp1[2])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom(italic(pi)[B] * ", ") * italic(.(comp1[1])) * " to " * italic(.(comp1[2]))), line=2.5, cex.lab=label.size, col.lab=comp1[3])
  # Add cubic spline:
  lines(smooth.spline(x=Fst.vector1, y=Dxy.vector1, spar=1), lty = 1, col = alpha(comp1[3], 0.75), lwd=2)  # 0.8 is transparency
  # second plot:
  plot(x=Fst.vector2, y=Dxy.vector2, xlim=c(0,1), xaxp=c(0,1,2), ylim=c(0,high.Dxy), yaxp=c(0,high.Dxy,2), 
       pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic('F')[ST] * ", " * phantom(italic(.(comp2[1])) * " to " * italic(.(comp2[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic('F')[ST] * ", ") * italic(.(comp2[1])) * " to " * italic(.(comp2[2]))), line=2.5, cex.lab=label.size, col.lab=comp2[3])
  title(ylab=bquote(italic(pi)[B] * ", " * phantom(italic(.(comp2[1])) * " to " * italic(.(comp2[2])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom(italic(pi)[B] * ", ") * italic(.(comp2[1])) * " to " * italic(.(comp2[2]))), line=2.5, cex.lab=label.size, col.lab=comp2[3])
  lines(smooth.spline(x=Fst.vector2, y=Dxy.vector2, spar=1), lty = 1, col = alpha(comp2[3], 0.75), lwd=2)
  # third plot:
  plot(x=Fst.vector3, y=Dxy.vector3, xlim=c(0,1), xaxp=c(0,1,2), ylim=c(0,high.Dxy), yaxp=c(0,high.Dxy,2), 
       pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic('F')[ST] * ", " * phantom(italic(.(comp3[1])) * " to " * italic(.(comp3[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic('F')[ST] * ", ") * italic(.(comp3[1])) * " to " * italic(.(comp3[2]))), line=2.5, cex.lab=label.size, col.lab=comp3[3])
  title(ylab=bquote(italic(pi)[B] * ", " * phantom(italic(.(comp3[1])) * " to " * italic(.(comp3[2])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom(italic(pi)[B] * ", ") * italic(.(comp3[1])) * " to " * italic(.(comp3[2]))), line=2.5, cex.lab=label.size, col.lab=comp3[3])
  lines(smooth.spline(x=Fst.vector3, y=Dxy.vector3, spar=1), lty = 1, col = alpha(comp3[3], 0.75), lwd=2)
  dev.off()
  # do statistical tests of correlation:
  test1 <- cor.test(Fst.vector1, Dxy.vector1, method=cor.method)
  test2 <- cor.test(Fst.vector2, Dxy.vector2, method=cor.method)
  test3 <- cor.test(Fst.vector3, Dxy.vector3, method=cor.method)
  return(list(test1=test1, test2=test2, test3=test3))
}

comp1 <- c("RNSA", "RBSA", "#bda2c4")
comp2 <- c("YBSA", "RBSA", "#c69f6b")
comp3 <- c("YBSA", "RNSA", "#2f5734")
cor.method <- 'pearson'

plots3Fst_Dxy(comp1, comp2, comp3,
              autosome.genome.rolling.stats, cor.method)



# try to make scatterplots of Fst vs. piW for three group comparisons.
# each of the "comp_" objects should be have two group names and a color name, e.g.:
# comp1 <- c("troch", "vir", "green3") 
# returns the statistical correlation tests, using the cor.method as specified
plots3Fst_Dxy <- function(comp1, comp2, comp3,
                          autosome.genome.rolling.stats, cor.method) 
{
  pdf(file="Autosome_Fst_PiB.pdf", width=8, height=3)
  par(oma=c(1,3,1,3))  # set outer margins
  par(mfrow=c(1,3))
  label.size <- 1.5
  Fst.row.choice1 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == paste(comp1[1], comp1[2], sep="_"))
  Fst.row.choice2 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == paste(comp2[1], comp2[2], sep="_"))
  Fst.row.choice3 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == paste(comp3[1], comp3[2], sep="_"))
  Fst.vector1 <- autosome.genome.rolling.stats$WC84_Fst[Fst.row.choice1,]
  Fst.vector2 <- autosome.genome.rolling.stats$WC84_Fst[Fst.row.choice2,]
  Fst.vector3 <- autosome.genome.rolling.stats$WC84_Fst[Fst.row.choice3,]
  Dxy.row.choice1 <- which(rownames(autosome.genome.rolling.stats$Dxy) == paste(comp1[1], comp1[2], sep="_"))
  Dxy.row.choice2 <- which(rownames(autosome.genome.rolling.stats$Dxy) == paste(comp2[1], comp2[2], sep="_"))
  Dxy.row.choice3 <- which(rownames(autosome.genome.rolling.stats$Dxy) == paste(comp3[1], comp3[2], sep="_"))
  Dxy.vector1 <- autosome.genome.rolling.stats$Dxy[Dxy.row.choice1,]
  Dxy.vector2 <- autosome.genome.rolling.stats$Dxy[Dxy.row.choice2,]
  Dxy.vector3 <- autosome.genome.rolling.stats$Dxy[Dxy.row.choice3,]
  high.Dxy <- round(max(c(Dxy.vector1, Dxy.vector2, Dxy.vector3)), 3) + 0.001
  # first plot:
  plot(x=Fst.vector1, y=Dxy.vector1, xlim=c(0,1), xaxp=c(0,1,2), ylim=c(0,high.Dxy), yaxp=c(0,high.Dxy,2), 
       pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic('F')[ST] * ", " * phantom(italic(.(comp1[1])) * " to " * italic(.(comp1[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic('F')[ST] * ", ") * italic(.(comp1[1])) * " to " * italic(.(comp1[2]))), line=2.5, cex.lab=label.size, col.lab=comp1[3])
  title(ylab=bquote(italic(pi)[B] * ", " * phantom(italic(.(comp1[1])) * " to " * italic(.(comp1[2])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom(italic(pi)[B] * ", ") * italic(.(comp1[1])) * " to " * italic(.(comp1[2]))), line=2.5, cex.lab=label.size, col.lab=comp1[3])
  # Add cubic spline:
  lines(smooth.spline(x=Fst.vector1, y=Dxy.vector1, spar=1), lty = 1, col = alpha(comp1[3], 0.75), lwd=2)  # 0.8 is transparency
  # second plot:
  plot(x=Fst.vector2, y=Dxy.vector2, xlim=c(0,1), xaxp=c(0,1,2), ylim=c(0,high.Dxy), yaxp=c(0,high.Dxy,2), 
       pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic('F')[ST] * ", " * phantom(italic(.(comp2[1])) * " to " * italic(.(comp2[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic('F')[ST] * ", ") * italic(.(comp2[1])) * " to " * italic(.(comp2[2]))), line=2.5, cex.lab=label.size, col.lab=comp2[3])
  title(ylab=bquote(italic(pi)[B] * ", " * phantom(italic(.(comp2[1])) * " to " * italic(.(comp2[2])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom(italic(pi)[B] * ", ") * italic(.(comp2[1])) * " to " * italic(.(comp2[2]))), line=2.5, cex.lab=label.size, col.lab=comp2[3])
  lines(smooth.spline(x=Fst.vector2, y=Dxy.vector2, spar=1), lty = 1, col = alpha(comp2[3], 0.75), lwd=2)
  # third plot:
  plot(x=Fst.vector3, y=Dxy.vector3, xlim=c(0,1), xaxp=c(0,1,2), ylim=c(0,high.Dxy), yaxp=c(0,high.Dxy,2), 
       pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic('F')[ST] * ", " * phantom(italic(.(comp3[1])) * " to " * italic(.(comp3[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic('F')[ST] * ", ") * italic(.(comp3[1])) * " to " * italic(.(comp3[2]))), line=2.5, cex.lab=label.size, col.lab=comp3[3])
  title(ylab=bquote(italic(pi)[B] * ", " * phantom(italic(.(comp3[1])) * " to " * italic(.(comp3[2])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom(italic(pi)[B] * ", ") * italic(.(comp3[1])) * " to " * italic(.(comp3[2]))), line=2.5, cex.lab=label.size, col.lab=comp3[3])
  lines(smooth.spline(x=Fst.vector3, y=Dxy.vector3, spar=1), lty = 1, col = alpha(comp3[3], 0.75), lwd=2)
  dev.off()
  # do statistical tests of correlation:
  test1 <- cor.test(Fst.vector1, Dxy.vector1, method=cor.method)
  test2 <- cor.test(Fst.vector2, Dxy.vector2, method=cor.method)
  test3 <- cor.test(Fst.vector3, Dxy.vector3, method=cor.method)
  return(list(test1=test1, test2=test2, test3=test3))
}

comp1 <- c("RNSA", "RBSA", "#bda2c4")
comp2 <- c("YBSA", "RBSA", "#c69f6b")
comp3 <- c("YBSA", "RNSA", "#2f5734")
cor.method <- 'pearson'

plots3Fst_Dxy(comp1, comp2, comp3,
              autosome.genome.rolling.stats, cor.method)

