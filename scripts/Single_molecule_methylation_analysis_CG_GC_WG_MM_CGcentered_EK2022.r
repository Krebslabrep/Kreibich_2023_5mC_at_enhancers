#TITLE:   Single molecule CG GC methylation analysis with Cochran-Mantel-Haenszel
#AUTHOR:  Elisa Kreibich
#DATE:    27-10-2022
#VERSION: CG centric approach
#AIM:     Load CG-GC methylation calls in in/accessible fractions, load this for all replicates,
#           convert into 3-dimensional contingency table in an array format,
#           and calculate p-value and odds ratio with Cochran-Mantel-Haenszel test.
#           This is a proportionality test that accounts for the number of observations (in our case number of reads) by weighing the proportions
#           and computes an overall p value taking all replicates into account.
#           The odds ratio can be used instead of the methylation difference and tells how far away the ratios are from a random distribution between the fractions.


# Load arguments ---------------------------------------------------------------------
## Extract arguments when running on the cluster
args = commandArgs(trailingOnly=TRUE)
  if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
  } else if (length(args)==1) {
    # default output file
    args[2] = "out.txt"
  }
NAME        = args[1]  #celltype name (e.g. "ES", if sample name is "ES_R1_R2", or "ES_NO" if sample name is "ES_NO_R1_R2")
REP         = args[2]   #replicates (e.g. "R1_R2", if sample name is "ES_R1_R2")
CHR1        = args[3]  #which chromosome to start with (e.g. "1")
CHR2        = args[4]  #which chromosome to end with (e.g. "19")
PERFORM_CMH = args[5]  #logical argument whether to perform the statistical test ("TRUE"). If "FALSE" only coverage table will be made.
WDW_SIZE    = args[6]  #window size of the genomic bin (e.g. "101" for CA analysis)
DIR_WD      = args[7] #working directory (e.g. './Keibich_2022/')
DIR_QUASR   = args[8]  #input directory of QuasR files for qAlign (e.g. '/g/krebs/kreibich/HTS/SMF/MM/QuasR/QuasR_rmdup/')

## Define arguments when NOT running on the cluster
# NAME        = "ES_NO"
# REP         = "R1_R2_R5a6"
# CHR1        = "1"
# CHR2        = "19"
# PERFORM_CMH = TRUE
# WDW_SIZE    = 101
# DIR_WD      = './Kreibich_2023_5mC_at_enhancers/' #location of this repository
# DIR_QUASR   = 'data/example/' #location of QuasR files for qAlign

# Load libraries and dependencies -----------------------------------------------------
library(QuasR)
library("BSgenome.Mmusculus.UCSC.mm10")
library(tidyverse)
library(samplesizeCMH)
library(DescTools)

source('scripts/utilities/RFunctions_for_single_molecule_analysis.R')

INDIR   <- paste0(DIR_WD, 'data_results/single_molecule_call/')        # Input directory of Single Molecule methylation call CG-centered
OUTDIR  <- paste0(DIR_WD, 'data_results/single_molecule_call/')        # Output directory of Single Molecule methylation analysis CG-centered

# Load aligned data ---------------------------------------------------------------------
input_file <- paste0(DIR_QUASR, 'QuasR_aligned_files_', NAME, '_', REP, '_rmdup.txt')

Qproj <- qAlign(sampleFile = sample_file,
                genome ="BSgenome.Mmusculus.UCSC.mm10",
                paired = "fr",
                bisulfite = "undir")
Qproj@aligner <- "Rbowtie"
Qaln <- as.data.frame(alignments(Qproj)[[1]])
sampleNames <- Qaln[,2]     #e.g. SMF_MM_ES_NO_R1

# Run analysis ----------------------------------------------------------------------
## Arguments
cO            <- 30                                   # Coverage cutoff (>=)
wind_size     <- WDW_SIZE                             # Window size used for SM methylation call (101 bp for CA analysis)
SAMPLENUMBER  <- seq_along(sampleNames)
SAMPLENAMES   <- sampleNames[SAMPLENUMBER]
SAMPLEREPS    <- str_extract(SAMPLENAMES[-1], pattern = "R.{1,3}")
OUTPUTSUFFIX  <- paste0(SAMPLENAMES[1], '_', paste(SAMPLEREPS, collapse = '_'), '_', wind_size, 'bp')
minSample     <- 2                                    # Minimum number of samples needed for CMH (3 means there must be data for all three samples for all bins)
ADD           <- paste0("_min", minSample)  #""       # Additional info for output name, pasted before .rds
OUTPUTSUFFIX2 <- paste0('_Co', cO, '_pc', ADD)        # Complete info added after CMH_matrix before _chrXXX.rds

OUTPUTNAME    <- paste0('CMH_matrix_', OUTPUTSUFFIX, OUTPUTSUFFIX2)
OUTPUTNAME2   <- paste0('Coverage_matrix_', OUTPUTSUFFIX)
OUTPUTNAME3   <- paste0('CMH_matrix_BINNAMES_', OUTPUTSUFFIX)
OUTPUTNAME4   <- paste0('GpC_Coverage_matrix_', OUTPUTSUFFIX)


## Run loop over all chromosomes
for(CHR in seq(CHR1, CHR2)){
  # CHR         <- 11
  message(Sys.time())
  message("chr = ", CHR)
  message("SAMPLENAMES = ", paste(SAMPLENAMES, collapse = ", "))
  message("minSample = ", minSample)
  message("ADD = ", ADD)
  
  
  ## Load single molecule methylation data from samples per chromosome
    CG_GC_me_list <- lapply(SAMPLENUMBER, function(sp){
      print(SAMPLENAMES[sp])
      readRDS(paste0(INDIR, 'tmp_', SAMPLENAMES[sp], '_CG_binned_met_', wind_size, 'bp_chr', CHR,'.rds'))
    })
    message("Data loading complete")

  ## Extract CGmeth matrix from all samples
    CGme_list     <- lapply(seq_along(CG_GC_me_list), function(sample){CG_GC_me_list[[sample]][[2]][[1]]})
    names(CGme_list) <- SAMPLENAMES

  ## Make and save coverage tibble for CGmeth
    count_bin_tbl <- make.coverage.tibble(CGme_list)

    saveRDS(count_bin_tbl, paste0(OUTDIR, OUTPUTNAME2, '_chr', CHR, '.rds'))
    message(paste("Saving", paste0(OUTDIR, OUTPUTNAME2, '_chr', CHR, '.rds') ,"complete"))

if(PERFORM_CMH == TRUE){
  ## Extract bin names from all samples and filter for overlapping
    CGme_list1    <- lapply(seq_along(CGme_list), function(sample) {unique(CGme_list[[sample]])})
    BINNAMES_list <- lapply(seq_along(CGme_list), function(sample) {unique(names(CGme_list[[sample]]))})
    CGme_list1    <- lapply(seq_along(CGme_list), function(sample) {
                        names(CGme_list1[[sample]]) <- BINNAMES_list[[sample]]
                        return(CGme_list1[[sample]])
                      })
    BINNAMES      <- Reduce(intersect, BINNAMES_list) %>% unique()
  
  ## Filter for bins that are present in all samples
    CGme_list1.2  <- lapply(seq_along(CGme_list1), function(sample) {CGme_list1[[sample]][names(CGme_list1[[sample]]) %in% BINNAMES]})
    BINNAMES      <- names(CGme_list1.2[[1]])   # bug fix --> there are more bins than originally selected for
   
    saveRDS(BINNAMES, paste0(OUTDIR, OUTPUTNAME3, '_chr', CHR, '.rds'))
    message(paste("Saving", paste0(OUTDIR, OUTPUTNAME3, '_chr', CHR, '.rds'), "complete"))
    
  ## Make CMH input tibble - including filtering for coverage
    CGme_list2    <- filter.for.coverage(CGme_list1.2, cO, SAMPLENAMES)
    bin_tbl_list  <- make.CMH.input.tibble(CGme_list2, SAMPLENAMES, BINNAMES)
  
  ## Perform CMH test and make final tibble
    bin_arrays    <- make.3Dcontig.table(bin_tbl_list, BINNAMES)
    
    Final_tbl     <- perform.CMH.test(bin_arrays, BINNAMES)

    saveRDS(Final_tbl, paste0(OUTDIR, OUTPUTNAME, '_chr', CHR, '.rds'))
    message(paste("Saving", paste0(OUTDIR, OUTPUTNAME, '_chr', CHR, '.rds') ,"complete"))
}  
    
  ## Extract GCmeth matrix from all samples -  GpC coverage tibble
    GCme_list <- lapply(seq_along(CG_GC_me_list), function(sample){CG_GC_me_list[[sample]][[1]][[1]]})
    names(GCme_list) <- SAMPLENAMES
    
  ## Make and save coverage tibble for CGmeth
    GCcount_bin_tbl <- make.coverage.tibble(GCme_list) %>%
      mutate(fraction = str_extract(fraction, "\\d{1}$"))
    
    saveRDS(GCcount_bin_tbl, paste0(OUTDIR, OUTPUTNAME4, '_chr', CHR, '.rds'))
    message(paste("Saving", paste0(OUTDIR, OUTPUTNAME4, '_chr', CHR, '.rds') ,"complete"))   
    
    message(Sys.time())
}


# Combine all chromosomes to one big tibble -------------------------------------------------------------
## CMH tibble
  combine.chr.tibbles(OUTDIR, OUTPUTNAME)

## Coverage tibble CpGs
  combine.chr.tibbles(OUTDIR, OUTPUTNAME2)

## BINNAMES tibble
  combine.chr.tibbles(OUTDIR, OUTPUTNAME3)
  
## Coverage tibble GpCs
  combine.chr.tibbles(OUTDIR, OUTPUTNAME4)

message("Script finished")
