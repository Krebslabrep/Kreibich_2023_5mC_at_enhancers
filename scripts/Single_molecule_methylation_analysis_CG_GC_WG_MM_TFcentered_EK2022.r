#TITLE:   Single molecule CG GC methylation analysis with Cochran-Mantel-Haenszel - TF
#AUTHOR:  Elisa Kreibich
#DATE:    27-10-2022
#VERSION: TF-motif centric approach
#AIM:     Load CG-GC methylation calls in in/accessible fractions, load this for all replicates,
#           convert into 3-dimensional contingency table in an array format,
#           and calculate p-value and odds ratio with Cochran-Mantel-Haenszel test.
#           This is a proportionality test that accounts for the number of observations (in our case number of reads) by weighing the propotions
#           and computes a overall p value taking all replicates into account.
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
REP         = args[2]  #replicates (e.g. "R1_R2", if sample name is "ES_R1_R2")
CHR1        = args[3]  #which chromosome to start with (e.g. "1")
CHR2        = args[4]  #which chromosome to end with (e.g. "19")
PERFORM_CMH = args[5]  #logical argument whether to perform the statistical test ("TRUE"). If "FALSE" only coverage table will be made.
WDW_SIZE    = args[6]  #window size of the genomic bin (e.g. "101" for CA analysis, hard coded 30 bp for TFBS)
DIR_WD      = args[7]  #working directory (Github repository directory of Kreibich_2023_5mC_at_enhancers)
DIR_QUASR   = args[8]  #input directory of QuasR files for qAlign


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

INDIR   <- paste0(DIR_WD, 'data_results/single_molecule_call/')        # Input directory of Single Molecule methylation call TF-centered
OUTDIR  <- paste0(DIR_WD, 'data_results/single_molecule_call/')        # Output directory of Single Molecule methylation analysis TF-centered

# Load aligned data ---------------------------------------------------------------------
input_file <- paste0(DIR_QUASR, 'QuasR_aligned_files_', NAME, '_', REP, '_rmdup.txt')

Qproj <- qAlign(sampleFile = input_file,
                genome ="BSgenome.Mmusculus.UCSC.mm10",
                paired = "fr",
                bisulfite = "undir")
Qproj@aligner <- "Rbowtie"
Qaln <- as.data.frame(alignments(Qproj)[[1]])
sampleNames <- Qaln[,2]     #e.g. SMF_MM_ES_NO_R1

# # Position of the baits ---------------------------------------------------------------------
# baits_mm10=import('/g/krebs/krebs/DB/SureSelect/MouseMethyl_Bait_merged_mm10.bed')
# 
# ## TF bins - Jaspar2018
# MotDBf <- readRDS('/g/krebs/kreibich/DB/mapped_jaspar2018_ChIP_score10_inBaits_BANP.rds')
# # names(MotDBf) <- paste('TFBS_',seq_along(MotDBf),sep='')
# # MotDBf$TFBS <- names(MotDBf)
# MotDBf_tbl <- as_tibble(MotDBf)
# 
# col_wind <- MotDBf
# 
# ## AllC GRange object
# AllC <- readRDS('/g/krebs/kreibich/analysis/SMF/methCall/AllC_GRange_object.rds')
# AllC$cytosine <- c(1:length(AllC))
# CGi <- AllC$type == 'NWCGW'
# 
# ## Control bins for context
# col_wind_CGi <- countOverlaps(resize(col_wind, WDW_SIZE, fix='center'), resize(AllC[CGi], 1, fix = 'center'), type = 'equal', ignore.strand = TRUE) == 1
# 
# ## Map CGs bins back to AllC
# ov <- as.matrix(findOverlaps(resize(col_wind, WDW_SIZE, fix='center'), AllC, type = 'equal'))
# names(AllC)[ov[,2]] <- names(col_wind[ov[,1]])
# AllC2 <- AllC[!is.na(names(AllC))]
# AllC2$bin <- names(AllC2)
# 
# ## Create AllC tibble with bin name for every cytosine
# AllC_tbl <- as_tibble(AllC2)
# remove(CGi)
# remove(AllC)

# Compile methylation matrix for TFs as a function of states -----------------------------------
allPos <- expand.grid(c(0,1),c(0,1),c(0,1))
patternStrings <- names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))

states <- list(
  closed = patternStrings[c(1, seq_along(patternStrings)[!seq_along(patternStrings) %in% c(1, 6, 8)])],
  accessible = patternStrings[8],
  bound = patternStrings[6]
)                           # Using only 'pure' states

statesF <- as.factor(unlist(lapply(seq_along(states), function(i){rep(names(states[i]), length(states[[i]]))}))[order(unlist(states))])
statesF <- factor(statesF, levels = names(states))

states_tbl <- bind_rows(states) %>% 
  pivot_longer(cols = everything(), names_to = "state", values_to = "fraction") %>%
  unique()


# Run analysis ------------------------------------------------------------------
## Arguments
cO_TFBS       <- 30 # coverage cutoff < cO_TFBS
WDW_SIZE      <- WDW_SIZE   # Window size used for SM methylation call (30 bp for TFBS analysis)
SAMPLENUMBER  <- seq_along(sampleNames)
SAMPLENAMES   <- sampleNames[SAMPLENUMBER]
SAMPLEREPS    <- str_extract(SAMPLENAMES[-1], pattern = "R.{1,3}")
OUTPUTSUFFIX  <- paste0(SAMPLENAMES[1], '_', paste(SAMPLEREPS, collapse = '_'), '_TF_jaspar2018_BANP_inM', WDW_SIZE)
minSample     <- 2       # minimum number of samples needed for CMH (3 means there must be data for all three samples for all bins)
ADD           <- paste0("_min", minSample, "_closed2_NWCGW_2")  #""       # additional info for output name, pasted before .rds
OUTPUTSUFFIX2 <- paste0('_Co', cO_TFBS, '_pc', ADD) #Info added after CMH_matrix

OUTPUTNAME    <- paste0('TF_CMH_matrix_', OUTPUTSUFFIX, OUTPUTSUFFIX2)
OUTPUTNAME2   <- paste0('TF_Coverage_matrix_', OUTPUTSUFFIX)
WITH_COVERAGE_TABLE <- TRUE


## Run loop over all chromosomes
for(CHR in seq(CHR1, CHR2)){
  message(Sys.time())
  message("chr = ", CHR)
  message("SAMPLENAMES = ", paste(SAMPLENAMES, collapse = ", "))
  message("minSample = ", minSample)
  message("ADD = ", ADD)
  
  
  ## Load single molecule methylation data from samples per chromosome
    CG_GC_me_list <- lapply(SAMPLENUMBER, function(sp){
      print(SAMPLENAMES[sp])
      readRDS(paste0(INDIR, 'tmp_', SAMPLENAMES[sp], '_CG_TF_jaspar2018_ChIP_score10_inBaits_BANP_sorted_met_vect_inM', WDW_SIZE, '_chr', CHR,'_NWCGW.rds'))
    })
    message("Data loading complete")
  
  ## Extract CGmeth matrix from all samples
    CGme_list <- lapply(seq_along(CG_GC_me_list), function(sample) {CG_GC_me_list[[sample]][[2]][[2]]})
    names(CGme_list) <- SAMPLENAMES
    
  ## Make and save coverage tibble for CGmeth
    if(WITH_COVERAGE_TABLE == TRUE){
      count_bin_tbl <- make.coverage.tibble(CGme_list) %>% 
        mutate(fraction = str_remove(fraction, "TFBS_\\d{1,}_")) %>% 
        left_join(., states_tbl) %>% 
        select(-fraction) %>% 
        add_count(sample, bin, state, value, wt = n, name = "n") %>% 
        distinct() %>% 
        select(sample, bin, state, value, n)
      
      saveRDS(count_bin_tbl, paste0(OUTDIR, OUTPUTNAME2, '_chr', CHR, '.rds'))
      message(paste("Saving", paste0(OUTDIR, OUTPUTNAME2, '_chr', CHR, '.rds') ,"complete"))
    }

if(PERFORM_CMH == TRUE){
  ## Filter for coverage and remove bins with only one fraction - for every sample
    CGme_list2 <- filter.for.coverage(CGme_list, cO_TFBS, SAMPLENAMES, TF = TRUE)

  ## Remove NULL TFBS bins
    CGme_list2 <- lapply(seq_along(CGme_list2), function(sample){compact(CGme_list2[[sample]])})

  ## Convert to tibbles and filter for TFBS bin that contain bound fraction
    CGme_tbl_l <- make.CMH.input.tibble.TF(CGme_list2, SAMPLENAMES)

  ## Join tibbles for all samples to final tibble with state info (bound, accessible, closed)
    freq_tbl  <- bind_rows(CGme_tbl_l, .id = "sample") %>% 
      left_join(., states_tbl) %>% 
      select(-fraction)
  
  ## Compute 3-dimensional contingency tables in an array format for each bin
    BINNAMES <- pull(freq_tbl, TFBS) %>% unique()
    bin_arrays <- make.3Dcontig.table.TF(freq_tbl, BINNAMES)
  
  ## Perform Cochran-Mantel-Haenszel test 
    Final_tbl <- perform.CMH.test(bin_arrays, BINNAMES) %>% 
      dplyr::rename(TFBS = bin)
  
    saveRDS(Final_tbl, paste0(OUTDIR, OUTPUTNAME, '_chr', CHR, '.rds'))
    message(paste("Saving", paste0(OUTDIR, OUTPUTNAME, '_chr', CHR, '.rds') ,"complete"))
 }   
    message(Sys.time())
}


# Combine all chromosomes to one big tibble -------------------------------------------------------------
## CMH tibble
combine.chr.tibbles(OUTDIR, OUTPUTNAME)

## Coverage tibble
combine.chr.tibbles(OUTDIR, OUTPUTNAME2)


message("Script finished")
