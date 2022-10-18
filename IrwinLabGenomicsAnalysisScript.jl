# IrwinLabGenomicsAnalysisScript.jl

# Started by Darren Irwin on 2Aug2022 
# for purpose of speeding up our lab's genomics pipeline.
# Using the R code as a guideline, while making use of Julia's features
# This is the R file that this is based on: warbler_GBS_analysis_script_2021_reanalysis.R
# Based on logic and code used in these papers:
# Irwin, D.E., M. Alcaide, K.E. Delmore, J.H. Irwin, and G.L. Owens. 2016. Recurrent selection explains parallel evolution of genomic regions of high relative but low absolute differentiation in a ring species. Molecular Ecology 25: 4488-4507.
# Irwin, D.E., B. Milá, D.P.L. Toews, A. Brelsford, H.L. Kenyon, A.N. Porter, C. Grossen, K.E. Delmore, M. Alcaide, and J.H. Irwin. 2018. A comparison of genomic islands of differentiation across three young avian species pairs. Molecular Ecology 27: 4839-4855.

# Please cite one of the above paper if you use these scripts and/or functions.
#
# The starting point for this R script requires that the "012NA" files
# for each chromosome have already been produced, either by following the
# instructions in the "warbler_genomics_processing_scripts.txt",
# or by simply downloading the 012NA files from the Dryad package: https://doi.org/10.5061/dryad.4j2662g.

# I am happy to answer any questions: irwin@zoology.ubc.ca

# Set directories and functions ----

#import Pkg; Pkg.add("CSV")
#import Pkg; Pkg.add("DataFrames")
#import Pkg; Pkg.add("JLD2")
#import Pkg; Pkg.add("NaNStatistics")


using DelimitedFiles
using CSV
using DataFrames
using JLD2
using Statistics
using NaNStatistics

# set working directory (where files will be saved)
cd("/scratch/st-darreni-1/libby/WGS/allSites_012NA_renamed")

# choose directory of source data
source_data_dir = raw"/scratch/st-darreni-1/libby/WGS/allSites_012NA_renamed"

# Setup for main processing ----

# choose the chromosomes to analyze in this run
chromosomes_to_analyze = ARGS[1]


# to process other chromosomes, choose from among these:
# chromosomes_to_analyze = ["1","1A","1B","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28","Z"]

# Options to calculate the per-site and windowed stats, or to load already calculated stats:
calculate_or_load_stats = 1  # 1) calculate site stats;  
                              # 2) load previously calculated per-site stats; 
                              # 3) load per-site and windowed data from file
saveSiteInfo = true    # If TRUE, will save a file of per-site stats (set to FALSE if file already exists)
saveWindowedStats = true   # If TRUE, will save a file for per-window stats (set to FALSE if file already exists)

# choose path and filename for the 012NA files
baseName = "wgs.genotypes.allSites." * chromosomes_to_analyze * ".filtered.missing80_mindepth3"
tagName = ".wgs_w50000."   # choose a tag name for this analysis
# indicate name of metadata file, a text file with these column headings:
# ID	location	group	Fst_group	plot_order
metadataFile = "../wgs.Fst_groups.txt"
# load metadata
metadata = DataFrame(CSV.File(metadataFile)) # the CSV.File function interprets the correct delimiter
num_metadata_cols = ncol(metadata)
num_individuals = 78  # specify number of individuals in file (good for error checking)
# specify window size (number of bp with info) and step size
windowSize = 50000   
stepSize = windowSize  # could change if wanting overlapping windows (not built in yet to Julia version)
# specify groups for calculation of statistics (these are in "Fst_group" column in metadata file)
groups = ["YBSA", "YBxRB", "RNxYB", "RNSA", "RBxRN", "RBSA"]
#groups = ["MGWA", "MOWA", "TOWA"]
group_colors = ["#dddf5bff", "#c6924dff", "#2f5734ff", "#394ca2ff", "#b26ec4ff", "#bd5757ff"]  
group_count = length(groups)

# specify groups for plotting, and their colors
groups_to_plot_WC84_Fst = ["YBSA_RNSA",
                             "RNSA_RBSA",
                             "YBSA_RBSA"]
group_colors_WC84_Fst = ["#2f5734",
                         "#bda2c4",
                         "#c69f6b"]

groups_to_plot_Dxy = groups_to_plot_WC84_Fst   # or specify differences if necessary
group_colors_Dxy = group_colors_WC84_Fst  # or specify differences if necessary
groups_to_plot_pi = ["YBSA", "RNSA", "RBSA"]
group_colors_pi = ["#e6e045", "#394ca2", "#bd5757"]

# Option to focus on a region of chromosome ----
# (not presented in paper itself, but used for inspecting data in detail)
focus_region = false  # choose true for a subset of the chromosome, false for the whole thing)
if focus_region # meaning if focus_region == true
  position_min =  1500000
  position_max = 1750000
end

# note: leaving out the calculation of D-statistics, which is in the R code (maybe add here later)


# FUNCTIONS (load these before calling main loop far below) ----

# calculate sample size (of individuals) and frequency of alternate allele 
# (homozygote coded as 2; heterozygote as 1, ref homozygote as 0, missing as -1)
function getFreqsAndSampleSizes(genoData, indGroup, groupsToCalc) 
    groupCount = length(groupsToCalc)
    freqs = Array{Float32, 2}(undef, groupCount, size(genoData, 2))
    sampleSizes = Array{Int16, 2}(undef, groupCount, size(genoData, 2)) 
    for i in 1:groupCount
        selection = (indGroup .== groupsToCalc[i]) # gets the correct rows for individuals in the group 
        geno0counts = sum(genoData[selection,:] .== 0 , dims=1) # count by column the number of 0 genotypes (homozygous ref)
        geno1counts = sum(genoData[selection,:] .== 1 , dims=1) # same for 1 genotypes (heterozygous)
        geno2counts = sum(genoData[selection,:] .== 2 , dims=1) # same for 2 genotypes (homozygous alternate) 
        # genoMissingCounts = sum(genoData[selection,:] .== -1 , dims=1) # same for -1 genotypes (missing)
        sumGenoCounts = geno0counts .+ geno1counts .+ geno2counts 
        sampleSizes[i,:] = sumGenoCounts
        freqs[i,:] = ((2 .* geno2counts) .+ geno1counts) ./ (2 * sumGenoCounts)   
    end
    return freqs, sampleSizes
end 

# calculate non-biased nucleotide diversity (pi) at each site for each population,
# with correction for sample size, using allele freqs and sample sizes as input:
function getSitePi(freqs, sampleSizes)::Matrix{Float32}
    sitePi = @. 2freqs*(1 - freqs) * (2sampleSizes / (2sampleSizes - 1)) # the @. macro broadcasts the dot to each math operator in the formula
    # (technical note: there is never division by zero above, and division by -1 only occurs when freqs is already NaN)
end

# get row names for use in Dxy and Fst matrices, 
# or any others that compare populations pairwise:
function getPairwiseNames(groupsToCalc)::Vector{String}
    groupCount = length(groupsToCalc)
    pairwiseNames = []
    for i in 1:(groupCount-1)
        for j in (i+1):groupCount
            push!(pairwiseNames, string(groupsToCalc[i],"_",groupsToCalc[j]))
        end
    end
    return pairwiseNames
end 

#getPairwiseNames(groupsToCalc)

# calculate pairwise Dxy per site, using data in "freqs" and groups in "groups"
function getDxy(freqs, groupsToCalc)
    pairwiseNames = getPairwiseNames(groupsToCalc)
    groupCount = length(groupsToCalc)
    Dxy = Array{Float32, 2}(undef, length(pairwiseNames), size(freqs, 2))
    rowCount = 1
    for i in 1:(groupCount-1)
        for j in (i+1):groupCount
            Dxy[rowCount,:] = @. freqs[i,:]*(1-freqs[j,:]) + freqs[j,:]*(1-freqs[i,:])
            rowCount += 1 # adds one to rowCount
        end
    end
    return Dxy, pairwiseNames
end

# calculate Fst (and numerator and denominator) for each site (bp), 
# between pairs of groups (so pops (r) is 2), 
# using the Weir&Cockerham 1984 approach to correct for sample size and number of pops.
# Not yet ready: Set "among=true" if wanting it to include a line for Fst among all groups.
function getFst(freqs, sampleSizes, groupsToCalc; among=false)
    pairwiseNames = getPairwiseNames(groupsToCalc)
    groupCount = length(groupsToCalc)
    Fst = Array{Float32, 2}(undef, length(pairwiseNames), size(freqs, 2))
    FstNum = Array{Float32, 2}(undef, length(pairwiseNames), size(freqs, 2))
    FstDen = Array{Float32, 2}(undef, length(pairwiseNames), size(freqs, 2))
    rowCount = 1
    for i in 1:(groupCount-1)
        for j in (i+1):groupCount
            r = 2 # number of populations in comparison (from Weir & Cockerham 1984 p. 1363)
            n_bar = @. (sampleSizes[i,:] + sampleSizes[j,:]) / 2 # mean sample size 
            p_bar = @. (sampleSizes[i,:]*freqs[i,:] + sampleSizes[j,:]*freqs[j,:]) / (2*n_bar) # mean frequency
            C = @. sqrt(((sampleSizes[i,:] - n_bar)^2 + (sampleSizes[j,:] - n_bar)^2)/(r-1)) / n_bar  # The CV of sample sizes, based on the equation for sample SD--Mike & I solved the n_c equation in Weir&Cockerham1984 and confirmed this. (see the division by 1 in the numerator, prior to square root) 
            s2 = @. (sampleSizes[i,:]*((freqs[i,:] - p_bar)^2) + sampleSizes[j,:]*((freqs[j,:] - p_bar)^2)) / ((r-1) * n_bar)   # allele frequency variance, using sample variance, as per W&C84
            FstNum[rowCount,:] = @. s2 - ((p_bar*(1-p_bar) - ((r-1)/r) * s2)/(2 * n_bar - 1))
            FstDenTerm1 = @. (1 - ((2 * n_bar * C^2) / ((2 * n_bar - 1) *r ))) * p_bar * (1-p_bar)  
            FstDenTerm2 = @. (1 + ((2 * n_bar * (r-1) * C^2) / ((2 * n_bar - 1) * r))) * s2 / r
            FstDen[rowCount,:] = FstDenTerm1 .+ FstDenTerm2
            Fst[rowCount,:] = FstNum[rowCount,:] ./ FstDen[rowCount,:] 
            rowCount += 1 # adds one to rowCount
        end
    end
    return Fst, FstNum, FstDen, pairwiseNames
end 

# calculate windowed pi, in whole windows starting on left side of chromosome
function getWindowedPi(sitePi, pos, windowSize)
    windowCount = size(sitePi, 2) ÷ windowSize  # the division symbol ÷ produces the quotient (without the remainder)    
    windowedPi = Array{Float32, 2}(undef, size(sitePi, 1), windowCount)
    windowedPi_pos = Vector{Float32}(undef, windowCount) 
    for i in 1:windowCount
        windowedPi[:,i] = nanmean(sitePi[:, (i-1)*windowSize+1 : i*windowSize], dims=2)
        windowedPi_pos[i] = mean(pos.position[(i-1)*windowSize+1 : i*windowSize]) 
    end
    return windowedPi, windowedPi_pos     
end

# calculate windowed Dxy, in whole windows starting on left side of chromosome
function getWindowedDxy(Dxy, pos, windowSize)
    windowCount = size(Dxy, 2) ÷ windowSize  # the division symbol ÷ produces the quotient (without the remainder)    
    windowedDxy = Array{Float32, 2}(undef, size(Dxy, 1), windowCount)
    windowedDxy_pos = Vector{Float32}(undef, windowCount) 
    for i in 1:windowCount
        windowedDxy[:,i] = nanmean(Dxy[:, (i-1)*windowSize+1 : i*windowSize], dims=2)
        windowedDxy_pos[i] = mean(pos.position[(i-1)*windowSize+1 : i*windowSize]) 
    end
    return windowedDxy, windowedDxy_pos     
end 

# calculate windowed Fst according to according to Weir&Cockerham1984 
# (with sample size and pop number correction),
# calculated as windowed numerator over windowed denominator,
# in whole windows starting on left side of chromosome
function getWindowedFst(FstNum, FstDen, pos, windowSize)
    windowCount = size(FstNum, 2) ÷ windowSize  # the division symbol ÷ produces the quotient (without the remainder)    
    windowedFst = Array{Float32, 2}(undef, size(FstNum, 1), windowCount)
    windowedFst_pos = Vector{Float32}(undef, windowCount) 
    for i in 1:windowCount
        windowedFstNum = nansum(FstNum[:, (i-1)*windowSize+1 : i*windowSize], dims=2) 
        windowedFstDen = nansum(FstDen[:, (i-1)*windowSize+1 : i*windowSize], dims=2) 
        windowedFst[:,i] = windowedFstNum ./ windowedFstDen
        windowedFst_pos[i] = mean(pos.position[(i-1)*windowSize+1 : i*windowSize]) 
    end
    return windowedFst, windowedFst_pos     
end  


# MAIN LOOP NO LONGER LOOP----
chr = chromosomes_to_analyze
print(string("Starting chr ",chr,"\n"))

# Get chr data ----
# read in individual names for this chromosome dataset
individuals_file_name = string(baseName, ".012.indv")
ind = DataFrame(CSV.File(individuals_file_name; header=["ind"], types=[String])) 
indNum = size(ind, 1) # number of individuals
# read in position data for this chromosome
position_file_name = string(baseName, ".012.pos")
pos_whole_chr = DataFrame(CSV.File(position_file_name; header=["chrom", "position"], types=[String, Int]))
# read in genotype data for this chromosome
column_names = ["null"; string.("c.", pos_whole_chr.chrom, ".", pos_whole_chr.position)]    
genotype_file_name = string(baseName, ".012minus1") 
@time if 1 <= indNum <= 127   
    geno = readdlm(genotype_file_name, '\t', Int8, '\n'); # this has been sped up dramatically, by first coverting "NA" to -1
elseif 128 <= indNum <= 32767
    geno = readdlm(genotype_file_name, '\t', Int16, '\n'); # this needed for first column, which is number of individual; Int16 not much slower on import than Int8
else
    print("Error: Number of individuals in .indv appears outside of range from 1 to 32767")
end
loci_count = size(geno, 2) - 1   # because the first column is not a SNP (just a count from zero)
print(string("Read in genotypic data at ", loci_count," loci for ", indNum, " individuals. \n"))


# Add metadata ----
ind_with_metadata = hcat(ind, metadata) 
print(ind_with_metadata) 
print("\n")  # prints a line break 
print("Check first two columns to make sure the same \n")

# Note: in this Julia implementation, not combining the metadata and genotypes
# into one data structure, as I think will be faster as separate structures.
# Important to remember this throughout the rest of the script.

# Filter individuals ----
# If need to filter out individuals, based on low read number (or other reason):
filter = false
filter_out_inds = [20, 103] # if filtering out individuals, specify their row number here
if filter
    # Specify individuals to filter out:
    ind_with_metadata_indFiltered = ind_with_metadata[Not(filter_out_inds), :]
    geno_indFiltered = geno[Not(filter_out_inds), :]
else
    ind_with_metadata_indFiltered = ind_with_metadata
    geno_indFiltered = geno
end

# Get region text ----
if focus_region # if there is a specific region of focus on the chromosome (or scaffold)
    selection = position_min .<= pos_whole_chr.position .<= position_max
    pos = pos_whole_chr[selection,:]
    pushfirst!(selection, 1) # adds a 1 in first position, this needed because first column of geno is line number on import
    genoRegion = geno_indFiltered[:, selection]
    regionText = string("Chr",chr,"_from_",position_min,"_to_",position_max)
else # otherwise just grab the data for the whole chromosome 
    position_min = 1
    position_max = pos_whole_chr.position[length(pos_whole_chr.position)]
    pos = pos_whole_chr
    genoRegion = geno_indFiltered
    regionText = string("Chr",chr,"_whole")
end


# Make site stats ---- 
if calculate_or_load_stats == 1
    # Calculate allele freqs and sample sizes (use column Fst_group)
    @time freqs, sampleSizes = getFreqsAndSampleSizes(genoRegion[:,Not(1)], ind_with_metadata_indFiltered.Fst_group, groups) 
    print("Calculated population allele frequencies and sample sizes \n")

    # calculate unbiased nucleotide diversity (pi) at each site for each population
    @time sitePi = getSitePi(freqs, sampleSizes) 
    print("Calculated unbiased population pi values \n")

    # calculate Dxy (equivalent to pi_B) at each site, between pairs of groups
    @time Dxy, pairwiseNamesDxy = getDxy(freqs, groups)
    print("Calculated Dxy values \n")

    # calculate Fst (and numerator and denominator) for each site, 
    # between pairs of groups (so pops (r) is 2), 
    # using the Weir&Cockerham 1984 approach to correct for sample size and number of pops
    @time Fst, FstNumerator, FstDenominator, pairwiseNamesFst = getFst(freqs, sampleSizes, groups; among=true)  # set among to FALSE if no among Fst wanted (some things won't work without it) 
    print("Calculated Fst values \n")

    if saveSiteInfo == true   # save the per-site stats, if chosen to in Intro section
        filename = string(baseName,tagName,regionText,"_SiteStats.jld2")
        jldsave(filename; groups, freqs, sampleSizes, sitePi, Dxy, pairwiseNamesDxy, Fst, FstNumerator, FstDenominator, pairwiseNamesFst, regionText)
        print("Saved summary site stats \n")
        # WC84_Fst_among <- WC84_Fst[rownames(WC84_Fst)=="Fst_among",]
        # print(paste0(length(WC84_Fst_among)," markers in total (",sum(!is.na(WC84_Fst_among))," variable and ",sum(is.na(WC84_Fst_among)), " invariant)"))
    else print("Site stats not saved \n")
    end

elseif calculate_or_load_stats==2 || calculate_or_load_stats==3  # then load previously saved site stats 
    filename = string(baseName,tagName,regionText,"_SiteStats.jl")    
    load(filename)
    print("Loaded saved summary stats \n")
end

# Make windowed stats ---- 
if (calculate_or_load_stats == 1 || calculate_or_load_stats == 2) 
    
    # calculate windowed pi, in whole windows starting on left side of chromosome
    @time windowedPi, windowedPi_pos = getWindowedPi(sitePi, pos, windowSize)
    print("Calculated windowed Pi_within \n")
    
    # calculate windowed Dxy, in whole windows starting on left side of chromosome
    @time windowedDxy, windowedDxy_pos = getWindowedDxy(Dxy, pos, windowSize)
    print("Calculated windowed Pi_between \n") 
    
    # calculate windowed Fst according to according to Weir&Cockerham1984 
    # (with sample size and pop number correction),
    # calculated as windowed numerator over windowed denominator.
    @time windowedFst, windowedFst_pos = getWindowedFst(FstNumerator, FstDenominator, pos, windowSize) 
    print("Calculated windowed Fst \n")

    # save the rolling mean data, if chosen in Intro section
    if saveWindowedStats
        filename = string(baseName,tagName,regionText,"_window",windowSize,"_WindowStats.jld2") 
        jldsave(filename; groups, windowedPi, windowedPi_pos, pairwiseNamesDxy, windowedDxy, windowedDxy_pos, pairwiseNamesFst, windowedFst, windowedFst_pos, regionText)
        print("Saved windowed stats \n")
    end

elseif calculate_or_load_stats == 3 # load the window stats
    filename = string(baseName,tagName,regionText,"_window",windowSize,"_WindowStats.jld2")    
    load(filename)
    print("Loaded saved window stats \n")
end
