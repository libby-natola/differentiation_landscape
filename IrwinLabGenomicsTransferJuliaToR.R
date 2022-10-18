# script to import data from Julia output of IrwinLabGenomicsAnalysisScript.jl
# This will read in the data objects and give them names used in the R script.
# They can then be used to make graphs and summary statistics using the R code.
# By Darren Irwin, started 7August2022
# Adapted by Libby Natola 8Aug2022

# if rhdf5 not installed yet, run these two lines:
#install.packages("BiocManager") 
#BiocManager::install("rhdf5")

library("rhdf5")

# substitute in the correct working directory in this line"
setwd("./")

chromosomes.to.analyze <- commandArgs(trailingOnly = TRUE)
#chromosomes.to.analyze <- 4


# import windowed stats from Julia  ----

# put correct file name in next line:
#file_from_Julia <- "allSites_012NA_renamed/wgs.genotypes.allSites.239.filtered.missing80_mindepth3.wgs_w50000.Chr239_whole_window50000_WindowStats.jld2"
file_from_Julia <- paste0("allSites_012NA_renamed/wgs.genotypes.allSites.",chromosomes.to.analyze,".filtered.missing80_mindepth3.wgs_w50000.Chr",chromosomes.to.analyze,"_whole_window50000_WindowStats.jld2")




# The below will read the file above and put the saved Julia objects into the R objects. 
# IMPORTANT NOTE: In the below, I have preceded the R object name with "JLD." so that 
# running the below will not overwrite things produced by the R script. 
# I encourage you to compare the R and Julia results, and when confident you are ready to 
# have the below replace things in R memory, then just 
# do a search and replace on the below script, replacing "JLD." with nothing:

JLD.region.text <- h5read(file_from_Julia, "regionText")

# specify groups for calculation of statistics (these are in "Fst_group" column in metadata file)
JLD.groups <- c("YBSA", "YBxRB", "RNxYB", "RNSA", "RBxRN", "RBSA")
JLD.group.colors <- c("#dddf5bff", "#c6924dff", "#2f5734ff", "#394ca2ff", "#b26ec4ff", "#bd5757ff") 
JLD.group_count <- length(JLD.groups)

JLD.groups.to.plot.pi <- c("YBSA", "RNSA", "RBSA")
JLD.rolling.mean.pos.pi <- h5read(file_from_Julia, "windowedPi_pos") # looks just like R object
JLD.rolling.mean.pi <- h5read(file_from_Julia, "windowedPi") # good but need to add row and column names
colnames(JLD.rolling.mean.pi) <- round(JLD.rolling.mean.pos.pi)
JLD.groups.pi <- h5read(file_from_Julia, "groups") 
rownames(JLD.rolling.mean.pi) <- JLD.groups.pi # now this looks just like R object
JLD.group.colors.pi <- c("#e6e045", "#394ca2", "#bd5757")

JLD.groups.to.plot.Dxy <- c("YBSA_RNSA",
                            "RNSA_RBSA",
                            "YBSA_RBSA"
)
JLD.rolling.mean.pos.Dxy <- h5read(file_from_Julia, "windowedDxy_pos") # looks just like R object
JLD.rolling.mean.Dxy <- h5read(file_from_Julia, "windowedDxy") # good but need to add row and column names
colnames(JLD.rolling.mean.Dxy) <- round(JLD.rolling.mean.pos.Dxy)
JLD.groups.Dxy <- h5read(file_from_Julia, "pairwiseNamesDxy") # looks just like R object
rownames(JLD.rolling.mean.Dxy) <- JLD.groups.Dxy # now this looks just like R object
JLD.group.colors.Dxy <- c("#2f5734",
                           "#bda2c4",
                           "#c69f6b"
)

JLD.groups.to.plot.WC84_Fst <- c("YBSA_RNSA",
                                 "RNSA_RBSA",
                                 "YBSA_RBSA"
)
JLD.rolling.mean.pos.WC84_Fst <- h5read(file_from_Julia, "windowedFst_pos") # looks just like R object
JLD.rolling.mean.WC84_Fst <- h5read(file_from_Julia, "windowedFst") # good but need to add row and column names
colnames(JLD.rolling.mean.WC84_Fst) <- round(JLD.rolling.mean.pos.WC84_Fst)
JLD.groups.WC84_Fst <- h5read(file_from_Julia, "pairwiseNamesFst") # looks just like R object
rownames(JLD.rolling.mean.WC84_Fst) <- JLD.groups.WC84_Fst # now this looks just like R object
JLD.group.colors.WC84_Fst <- c("#2f5734",
                           "#bda2c4",
                           "#c69f6b"
)


# Save or load ----
base.name <- paste0("allSites_012NA_renamed/wgs.genotypes.allSites.",chromosomes.to.analyze,".filtered.missing80_mindepth3.wgs")
tag.name <- ".wgs_w50000."   # choose a tag name for this analysis
window_size <- 50000
# save the rolling mean data
save(JLD.rolling.mean.pos.pi, JLD.rolling.mean.pi, JLD.rolling.mean.pos.Dxy, JLD.rolling.mean.Dxy, JLD.rolling.mean.pos.WC84_Fst, JLD.rolling.mean.WC84_Fst, file=paste0(base.name,tag.name,JLD.region.text,"_window",window_size,"_WindowStats_from_R.R"))
print("Saved rolling mean stats")
