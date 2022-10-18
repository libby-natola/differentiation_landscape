# differentiation_landscape
scripts and data associated with dissertation chapter 4, future publications for sapsucker differentiation work
"Genomic divergence and speciation in red-naped, red-breasted, and yellow-bellied sapsuckers " readme

GENERAL INFORMATION

1. Title of Dataset: 
"Genomic divergence and speciation in red-naped, red-breasted, and yellow-bellied sapsuckers " dataset

2. Author Information:
Name: Libby Natola 
Institution: University of British Columbia
ORCID: 0000-0001-7576-2825
Email: libby.natola@gmail.com

3. Description of dataset
This readme details files and analyses used in "Genomic divergence and speciation in red-naped, red-breasted, and yellow-bellied sapsuckers " for Molecular Ecology, CITATION HERE

4. Date and geographic location of Data Collection:
Data were collected by Libby Natola from 2019-2022 at the University of British Columbia

SHARING/ACCESS INFORMATION

###1. Licenses/restrictions placed on data


FILE LIST
Genomic data filtering and processing:
1. WGS_processing_only_clean.txt scripts to process datasets. Analyses run in UBC's ARC Sockeye HPC cluster
2. sample_info.txt sample metadata
3. pseudochromosomes.fasta Red-breasted Sapsucker (RBSA, Sphyrapicus ruber) Pac-Bio assembly mapped to Golden-Fronted Woodpecker (GFWO, Melanerpes aurifrons) from Wiley and Miller 2020.
4. bamqc_output_conversion.txt reformat bamqc pdf data to txt for R
5. bamqc_fixed.csv bamqc data files
6. bamqc_explore.R R scripts to evaluate coverage, sex chromosomes
7. IrwinLabGenomicsAnalysisScript.jl Julia to R script
8. IrwinLabGenomicsTransferJuliaToR.R converts Julia output files into R readable files
9. plotGenomeFstDxyPi.processing.txt run and loop the Julia and Julia to R scripts
10. ld.heatmap.scripts.txt scripts to calculate LD statistics among and within species

Analyses/Figures:
1. WGS_sample_map.R R scripts for map, run on local disk, Figure 1 in text finalized in inkscape (wgs_sample_map.svg)
2. PCA_manhattan.R R scripts for pca and manhattan, run on HPC (in R with conda wgs environment, not by calling R script), Figures 2 and 3 in text made using wgs_pca_1_2_nowisa_shapes.pdf, RBSA_RNSA_25kb_man.png, YBSA_RBSA_25kb_man.png, YBSA_RNSA_25kb_man.png finalized in inkscape wgs_pca_1_2_nowisa_shapes.svg and all_spp_25kb.svg
3. plotGenomeFstDxyPi2.R R scripts for Fst, πB, πW plots, run on HPC, Figure 4 in text made using pi_Fst_Dxy_means_1-45_sexes.pdf finalized in inkscape Fst_pib_piw.svg
4. genomics_R_functions_V3.R R source scripts for plotGenomeFstDxyPi2.R
5. Autosomal_Z_comparisons.R R scripts for autosome x Z chromosome Fst x πB plots and Welch's tests, Table 1 in text and Figure 5 in text made using Fst_PiB_autosomal_Z_RNSA_RBSA.pdf, Fst_PiB_autosomal_Z_YBSA_RBSA.pdf, Fst_PiB_autosomal_Z_YBSA_RNSA.pdf finalized in inkscape Fst_Pib_autosomes_Z.svg 
6. genomics_R_functions_V5.R R source scripts for Autosomal_Z_comparisons.R
3. ld.heatmap.chr.R script to make heatmaps by chromosome or scripts/ld.heatmap.thin.chr.R for scaffolds > 50k SNPs, species pair, run on HPC (called using Rscript ld.heatmap.chr.R option option), Figure 6 in text made using *sa_Pseudo1.LD.png #sa_Pseudo1030.LD.png finalized in inkscape
7. Pseudo1030*SNPs.R, Pseudo1*SNPs.R R scripts for genotype by individual plots, run on HPC, Figure 7 in text made using Pseudo1*SNPs_*only.pdf Pseudo1030*SNPs_*only.pdf finalized in inkscape Pseudo1_GxI_plot.svg and Pseudo1030_GxI_plot.svg
8. genomics_R_functions_V2.R R source scripts for genotype by individual plots
