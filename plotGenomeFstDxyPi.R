# script to make plotGenomeFstDxyPi plots 
# adapted by Libby Natola from Darren Irwin's scripts from 2016 paper
# Started 12 August 2022

#setwd("~/Documents/UBC/Bioinformatics/wgs")
#setwd("./")

# Load functions
source("genomics_R_functions_V3.R")

# install.packages("vroom")
library(vroom)   # for fread function for fast reading of data files

# choose the chromosomes to analyze in this run
#chromosomes.to.analyze <- c("wgs.genotypes.allSites.PseudoWWNC01000346.1_Melanerpes_aurifrons_Melanerpes_aurifrons_OMNH24340_contig_346_whole_genome_shotgun_sequence.filtered.missing80_mindepth3.012")
# to process other chromosomes, choose from among these:
chromosomes.to.analyze <- c(1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 1003, 1030, 1044, 1048, 1059, 1073, 1076, 1080, 1082, 1095, 1097, 1099, 1100, 1144)

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

source("genomics_R_functions_V3.R")
chromosomes.to.plot <- chromosomes.to.analyze
max_Dxy_axis <-  0.015
max_pi_axis <- 0.015
transparency_Fst_Dxy <- 0
transparency_pi <- 0
line_transparency <- 0.5
# plotGenomeFstDxyPi(base.name, tag.name, window_size, chromosomes.to.plot, 
#                    max_Dxy_axis, max_pi_axis, transparency_Fst_Dxy, transparency_pi, line_transparency,
#                    groups.to.plot.WC84_Fst, group.colors.WC84_Fst,
#                    groups.to.plot.Dxy, group.colors.Dxy,
#                    groups.to.plot.pi, group.colors.pi)      


#################
plotGenomeFstDxyPi <- function(base.name, tag.name, window_size, chromosomes.to.plot, 
                               max_Dxy_axis, max_pi_axis, transparency_Fst_Dxy, transparency_pi, line_transparency,
                               groups.to.plot.WC84_Fst, group.colors.WC84_Fst,
                               groups.to.plot.Dxy, group.colors.Dxy,
                               groups.to.plot.pi, group.colors.pi) {
  #plot all chromosomes in one quartz window
  # plot dimensions in inches:
  window.width <- 11.97
  window.height <- 12
  # numbers below are in inches:
  left.margin <- 1 / window.width  # the left outer margin of the figure (number is inches)
  right.margin <- 1 / window.width
  top.margin <- 0.5 / window.height
  bottom.margin <- 0.5 / window.height
  row.height <- 1 / window.height  # height of each row of the figure (row contains Fst, Dxy, and pi plot)
  plot.height <- row.height / 3
  row.space <- 0.25 / window.height  # gap between rows 
  bp.per.inch <- 20000000  # number of bp of sequence per inch
  gap.between.chr <- 0.6 # gap in inches between chromosomes plots
  gap.between.chr.bp <- gap.between.chr*bp.per.inch
  
  # get sizes of chromosomes (note this actually isn't the true length, just the mean position of the rightmost window):
  chr.length <- NULL	
  for (i in 1:length(chromosomes.to.plot)) {
    chr.text <- paste0("Chr",chromosomes.to.plot[i], "_whole", sep="")
    #load(paste0("~/Dropbox/Darren's current work/genomics training 2014/GW_genomics_Darren/GW_SNP_summary_stats_by_chromosome_from_R/GW_SNP_RollingMean_stats_GW_Lane5_plus_Liz.45samples.troch_vir_plumb.",load.text,"_from_R.R"))
    load(paste0(base.name[i],tag.name,chr.text,"_window",window_size,"_WindowStats_from_R.R"))
    #chr.length[i] <- rolling.mean.pos.Fst[length(rolling.mean.pos.Fst)]
    chr.length[i] <- JLD.rolling.mean.pos.WC84_Fst[length(JLD.rolling.mean.pos.WC84_Fst)]
  }
  
  chr.plotted.already <- rep(FALSE, length(chromosomes.to.plot))
  
  # make plot order by row
  inch.per.row <- (window.width-(left.margin+right.margin)*window.width)
  bp.per.row <- bp.per.inch * inch.per.row
  row <- 1
  chr.in.row <- matrix(NA, 20, 20)
  chr.bp.start.in.row <- matrix(NA, 20, 20)
  bp.start <- 0
  remaining.row.length <- bp.per.row
  while (sum(chr.plotted.already==FALSE) > 0) {  # repeat until all chr have a row to be plotted in
    place.in.row <- 1
    # look through chromosomes
    for (i in 1:length(chromosomes.to.plot)) { ### LN ADDED GAP.BETWEEN.CHR.BP 
      # if not plotted and short enough, add to row:
      if (chr.plotted.already[i]==FALSE && chr.length[i]<=remaining.row.length) {
        chr.in.row[row, place.in.row] <- chromosomes.to.plot[i]
        chr.bp.start.in.row[row, place.in.row] <- bp.start
        remaining.row.length <- remaining.row.length - chr.length[i] - gap.between.chr.bp
        chr.plotted.already[i] <- TRUE
        place.in.row <- place.in.row+1
        bp.start <- bp.start + chr.length[i] + gap.between.chr.bp
      }		
    }
    row <- row+1
    bp.start <- 0
    remaining.row.length <- bp.per.row
}
  
  # Plot in the order defined above
  # png(file=("pi_Fst_Dxy_means_test.png"), width = 1920, height = 1920,
  #     units = "px", pointsize = 3, res = 600)

#pdf(file=("pi_Fst_Dxy_means_1-45_sexes.pdf"), width = 15, height = 15)
  
#quartz(title="WC84_Fst, Dxy, and pi for all chromosomes", width=window.width, height=window.height)  # width and height in inches
  # plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), xaxt='n', yaxt='n', xlab=NULL, ylab= NULL, bty='n', las=1)
  Sys.sleep(0.1)  # this is to prevent a bug
  num.rows <- sum(!is.na(chr.in.row[,1]))
  new.fig.parameter <- FALSE
  for (row.num in 1:num.rows) {
    num.chr.in.row <- sum(!is.na(chr.in.row[row.num,]))
    for (order.num in 1:num.chr.in.row) {
      chr.text <- paste0("Chr",chr.in.row[row.num,order.num], "_whole", sep="")
      #load(paste0("~/Dropbox/Darren's current work/genomics training 2014/GW_genomics_Darren/GW_SNP_summary_stats_by_chromosome_from_R/GW_SNP_RollingMean_stats_",load.text,"_from_R.R"))
      load(paste0("allSites_012NA_renamed/wgs.genotypes.allSites.",chr.in.row[row.num,order.num],".filtered.missing80_mindepth3.wgs",tag.name,chr.text,"_window",window_size,"_WindowStats_from_R.R"))
      print(paste0("Loaded saved rolling mean stats for Chr ", chr.in.row[row.num,order.num], sep=""))
      plot.bp.length <- chr.length[chromosomes.to.plot==chr.in.row[row.num,order.num]]
      # plot WC84_Fst:
      fig.left.loc <- left.margin + (chr.bp.start.in.row[row.num,order.num]/bp.per.row)*(inch.per.row/window.width)   
      fig.right.loc <- fig.left.loc + (plot.bp.length/bp.per.row)*(inch.per.row/window.width)
      par(fig=c(fig.left.loc, fig.right.loc, 1-(top.margin+plot.height+(row.num-1)*(row.height+row.space)), 1-(top.margin+(row.num-1)*(row.height+row.space))), new = new.fig.parameter, mai=c(0,0,0,0), cex=0.5)  # mai sets the margins of the plot to zero, matching
      new.fig.parameter <- TRUE
      #upper.xlim <- bp.per.inch*(window.width-(left.margin+right.margin)*window.width)
      plot(x=NULL, y=NULL, xlim=c(0, plot.bp.length), ylim=c(0, 1.2), yaxp=c(0, 1, n=1), ylab=expression(italic('F')[ST]*"  "), xlab=NA, xaxt='n', bty='n', las=1, xaxs="i", tcl=-0.25, xpd=NA, cex=1, cex.lab=1.55, cex.axis=0.75)
      text(-2000000, 1.4, chr.in.row[row.num,order.num], pos=4, cex=2, xpd=NA)
      for (j in 1:length(groups.to.plot.WC84_Fst)) {
        color.rgb <- col2rgb(group.colors.WC84_Fst[j]) /255   # divide to convert color scales to 0-1
        lines(JLD.rolling.mean.pos.WC84_Fst, JLD.rolling.mean.WC84_Fst[rownames(JLD.rolling.mean.WC84_Fst)==groups.to.plot.WC84_Fst[j],],
              col=rgb(color.rgb[1], color.rgb[2], color.rgb[3], alpha=line_transparency))
        xx <- c(JLD.rolling.mean.pos.WC84_Fst[1], JLD.rolling.mean.pos.WC84_Fst, JLD.rolling.mean.pos.WC84_Fst[length(JLD.rolling.mean.pos.WC84_Fst)])
        yy <- c(0, JLD.rolling.mean.WC84_Fst[rownames(JLD.rolling.mean.WC84_Fst)==groups.to.plot.WC84_Fst[j],], 0)
        polygon(xx, yy, col=rgb(color.rgb[1], color.rgb[2], color.rgb[3], alpha=transparency_Fst_Dxy), border=NA)
      }
      # plot Dxy:
      y_label <- expression(pi[B]*"    ")                     # y_label <- expression(italic('D')[xy]*"    ")
      par(fig=c(fig.left.loc, fig.right.loc, 1-(top.margin+2*plot.height+(row.num-1)*(row.height+row.space)), 1-(top.margin+plot.height+(row.num-1)*(row.height+row.space))), new = TRUE, mai=c(0,0,0,0))  # mai sets the margins of the plot to zero, matching
      plot(x=NULL, y=NULL, xlim=c(0, plot.bp.length), ylim=c(0, max_Dxy_axis), yaxp=c(0, max_Dxy_axis*5/8, n=1), ylab=y_label, xlab=NA, xaxt='n', bty='n', las=1, xaxs="i", tcl=-0.25, xpd=NA, cex.lab=1.5, cex.axis=0.75)  # 
      for (j in 1:length(groups.to.plot.Dxy)) {
        color.rgb <- col2rgb(group.colors.Dxy[j]) /255   # divide to convert color scales to 0-1
        lines(JLD.rolling.mean.pos.Dxy, JLD.rolling.mean.Dxy[rownames(JLD.rolling.mean.Dxy)==groups.to.plot.Dxy[j],], 
              col=rgb(color.rgb[1], color.rgb[2], color.rgb[3], alpha=line_transparency))
        xx <- c(JLD.rolling.mean.pos.Dxy[1], JLD.rolling.mean.pos.Dxy, JLD.rolling.mean.pos.Dxy[length(JLD.rolling.mean.pos.Dxy)])
        yy <- c(0, JLD.rolling.mean.Dxy[rownames(JLD.rolling.mean.Dxy)==groups.to.plot.Dxy[j],], 0)
        polygon(xx, yy, col=rgb(color.rgb[1], color.rgb[2], color.rgb[3], alpha=transparency_Fst_Dxy), border=NA)
      }
      # plot pi: 
      par(fig=c(fig.left.loc, fig.right.loc, 1-(top.margin+3*plot.height+(row.num-1)*(row.height+row.space)), 1-(top.margin+2*plot.height+(row.num-1)*(row.height+row.space))), new = TRUE, mai=c(0,0,0,0))  # mai sets the margins of the plot to zero, matching
      plot(x=NULL, y=NULL, xlim=c(0, plot.bp.length), ylim=c(0, max_pi_axis), yaxp=c(0, max_pi_axis*5/8, n=1), ylab=expression(pi[W]*"    "), xlab=NA, xaxt='n', bty='n', las=1, xaxs="i", tcl=-0.25, xpd=NA, cex.lab=1.5, cex.axis=0.75)
      title(xlab="Location", outer=TRUE)
      for (j in 1:length(groups.to.plot.pi)) {
        color.rgb <- col2rgb(group.colors.pi[j]) /255   # divide to convert color scales to 0-1
        lines(JLD.rolling.mean.pos.pi, JLD.rolling.mean.pi[rownames(JLD.rolling.mean.pi)==groups.to.plot.pi[j],], 
              col=rgb(color.rgb[1], color.rgb[2], color.rgb[3], alpha=line_transparency))
        xx <- c(JLD.rolling.mean.pos.pi[1], JLD.rolling.mean.pos.pi, JLD.rolling.mean.pos.pi[length(JLD.rolling.mean.pos.pi)])
        yy <- c(0, JLD.rolling.mean.pi[rownames(JLD.rolling.mean.pi)==groups.to.plot.pi[j],], 0)
        polygon(xx, yy, col=rgb(color.rgb[1], color.rgb[2], color.rgb[3], alpha=transparency_pi), border=NA)
      }
    }
  }
}                

# the main script for Fig. 5:
plotGenomeFstDxyPi(base.name, tag.name, window_size, chromosomes.to.plot, 
                   max_Dxy_axis, max_pi_axis, transparency_Fst_Dxy, transparency_pi, line_transparency,
                   groups.to.plot.WC84_Fst, group.colors.WC84_Fst,
                   groups.to.plot.Dxy, group.colors.Dxy,
                   groups.to.plot.pi, group.colors.pi)
#dev.off()
