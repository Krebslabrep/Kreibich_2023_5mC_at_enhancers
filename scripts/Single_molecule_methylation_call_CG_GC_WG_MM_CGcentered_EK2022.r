#TITLE:   WG MM CG single molecule methylation call for CpG and GpC
#AUTHOR:  Arnaud Krebs & Elisa Kreibich
#DATE:    08-06-2022
#VERSION: CA analysis strategy
#AIM:     Making the single molecule methylation call for CGs and GCs in a chromatin accessibility (CA) centered on a CpG fashion.
#         The output file is a multi-level list object with the following structure: [[C-context (1:GC|2:CG)]] [[bin(CG window)]] [[accessibility state(0|1)]] vector of methylation values of single molecule reads.

# Load libraries and dependencies ---------------------------------------------------
library(QuasR)
library("BSgenome.Mmusculus.UCSC.mm10")

source('scripts/utilities/RFunctions_for_single_molecule_call.R')

# Load arguments ---------------------------------------------------------------------
# Extract arguments when running on the cluster
args = commandArgs(trailingOnly=TRUE)
  if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
  } else if (length(args)==1) {
    # default output file
    args[2] = "out.txt"
  }
  INPUT_QUASR = args[1]    #QuasR input file for qAlign (tab-delimted file with FileName (location of BAM file) and SampleName) | best only one sample/ replicate per input file | see example file
  WDW_SIZE    = args[2]    #101 bp for CA analysis used in Kreibich et al, 2023
  CHR1        = args[3]    #chromosome to start with (e.g. "1")
  CHR2        = args[4]    #chromosome to end with (e.g. "19")
  DIR_WD      = args[5]    #working directory (Github repository directory of Kreibich_2023_5mC_at_enhancers)
  INPUT_BAITS = args[6]    #input bed file for baits regions

## Define arguments when NOT running on the cluster
  # INPUT_QUASR = 'data/example/QuasR_aligned_files_ES_NO_R1_R2_R5a6_rmdup.txt'
  # WDW_SIZE    = 101
  # CHR1        = 1
  # CHR2        = 19
  # DIR_WD      = './Kreibich_2023_5mC_at_enhancers/' #location of this repository
  # INPUT_BAITS = 'data/MouseMethyl_Bait_merged_mm10.bed'

# Load environment ---------------------------------------------------------------------
dir.create(path = paste0(DIR_WD, 'data_results/')
dir.create(path = paste0(DIR_WD, 'data_results/single_molecule_call/')
OUTDIR   <- paste0(DIR_WD, 'data_results/single_molecule_call/')
OUT_BINS <- paste0(DIR_WD, 'data/CGs_over_baits_', WDW_SIZE,'bp.rds')


# Load aligned data ---------------------------------------------------------------------
Qproj = qAlign(sampleFile = INPUT_QUASR, 
                genome = "BSgenome.Mmusculus.UCSC.mm10", 
                paired = "fr", 
                bisulfite = "undir")
Qproj@aligner = "Rbowtie"
sampleNames = unique(alignments(Qproj)[[1]][,2])
message("Aligned data loaded")

# Step 1: Create CG centered collection windows ---------------------------------------------------------------------

if(!file.exists(OUT_BINS)){
  ## Load baits regions -> extending
  baits_mm10 = import(INPUT_BAITS)
  baits_mm10_ext = GRanges(seqnames(baits_mm10), IRanges(start(baits_mm10) - 500, end(baits_mm10) + 500))
  
  ## Get location of all CGs genome-wide
  musmus_chr <- paste0('chr', c(seq(19), 'X', 'Y'))                                                    # Length of all chromosome of Mus Musculus
  mm10_CGs <- mclapply(seq_along(musmus_chr), function(i){
    reg <- musmus_chr[i]                                                                              # Define which chromosome to look at
    CGs_chr_XSV <- matchPattern(DNAString("CG"), Mmusculus[[reg]], fixed = "subject")                 # Find positions of all CGs in the chromosome
    CGs_chr <- GRanges(seqnames = Rle(as.character(reg)), ranges = IRanges(start(CGs_chr_XSV), end = end(CGs_chr_XSV)), strand = "+")  # Create GRange object of CGs of the chromosome
    return(CGs_chr)  
  }, mc.cores = 10)
  mm10_CGs <- unlist(GRangesList(mm10_CGs))                                                           # Get CG position info from all chromosomes into on GRange object
  mm10_CGs <- resize(mm10_CGs, 1, fix = "center")                                                     # Extract position of only the C
  
  
  ## Get location of all CGs genome-wide (NWCGW only)
  musmus_chr <- paste0('chr', c(seq(19), 'X', 'Y'))														                      # Length of all chromosome of Mus Musculus
  mm10_NWCGWs <- mclapply(seq_along(musmus_chr), function(i){
    reg <- musmus_chr[i]
    NWCGWs_chr_XSV <- matchPattern(DNAString("NWCGW"), Mmusculus[[reg]], fixed = "subject")
    NWCGWs_chr <- GRanges(seqnames = Rle(as.character(reg)), ranges = IRanges(start(NWCGWs_chr_XSV), end = end(NWCGWs_chr_XSV)), strand = "+")
    return(NWCGWs_chr)  
  }, mc.cores = 10)
  mm10_NWCGWs <- unlist(GRangesList(mm10_NWCGWs))
  mm10_NWCGWs <- resize(mm10_NWCGWs, 1, fix = "center")
  
  
  ## Subset to the ones in baits
  regs <-  subsetByOverlaps(mm10_NWCGWs, baits_mm10_ext, ignore.strand = T)
      #total CGs 21.8M
      #6.1M CGs considered  (30%)
      #1M within baits (5%)
  
  ## Create X bp collection windows around the CGs
  message('CG window size is ', WDW_SIZE)
  col_wind <- resize(regs, WDW_SIZE, fix='center')
  names(col_wind) <- paste0('bin_', seq_along(col_wind))
  
  saveRDS(col_wind, OUT_BINS)
  message('CG centered collection windows created')
} else {
  message('CG window size is ', WDW_SIZE)
  col_wind <- readRDS(OUT_BINS)
  message('CG centered collection windows loaded')
}


# Step 2: Sort CG windows ---------------------------------------------------------------------------------------------

wmMapped_SOu <- col_wind
wms <- (wmMapped_SOu[order(start(wmMapped_SOu))])
wms_c <- split(wms, seqnames(wms))    # Split by chr

regs <- unlist(GRangesList(lapply(seq_along(wms_c), function(i){
  x = wms_c[[i]]
  l = length(x)
  if(l>0){
    #cut by CG window position
    seq.i <- seq(1, l, round(l/25))
    st.seq <- start(x[seq.i])
    lowB <- st.seq[1:(length(st.seq)-1)]-1000
    lowB[lowB<0] <- 1
    highB <- st.seq[2:(length(st.seq))]+1000
    highB[highB>seqlengths(Mmusculus)[as.character(seqnames(wms_c[[i]])[1])]] <- seqlengths(Mmusculus)[as.character(seqnames(wms_c[[i]])[1])]
    gr <- GRanges(seqnames(x[1]), IRanges(lowB, highB))
    gr
  } else(GRanges())
})))
message('CG centered collection windows sorted')


# Step 3: Perform vectorial read extraction ------------------------------------------------------------------------------

## Define the sample you want to analyse
sampleNames <- unique(alignments(Qproj)[[1]][,2])
sp_sbs <- seq_along(sampleNames)
# sp_sbs <- seq_along(sampleNames)[seq_along(sampleNames) %in% c(grep(cellLine, sampleNames))]
message(sampleNames[sp_sbs])

## Read extraction
for(CHR in seq(CHR1, CHR2)){           # de-activate for loop if chromosome is defined by argument
  print(CHR)
  for(sp in sp_sbs){
    print(sampleNames[sp])
    
    regsD <- cbind(as.data.frame(regs)[,1:3], sample = sampleNames[sp])
    regsD <- regsD[regsD[,1] %in% paste0('chr', CHR), ]
    regsL <- split(regsD, 1:nrow(regsD))
    
    print(as.character(sampleNames[sp]))
    sRs <- mclapply(seq_along(regsL), function(i){
      sR <- extract_CG_GC_overs_CG_bins(regsL[[i]], Qproj, col_wind)  	
      sR
    }, mc.cores = 10)
  
  
    ### Merge the genomic bins
    sRs_t <- lapply(seq(2), function(context){
      lapply(seq(1), function(bin){
        sR_g <- lapply(seq(sRs), function(gbin){
          sRs[[gbin]][[context]][[bin]]
        })
        do.call(c, sR_g)
      })
    })
    
    saveRDS(sRs_t, paste0(OUTDIR, 'tmp_', sampleNames[sp], '_CG_binned_met_', WDW_SIZE, 'bp_chr', CHR,'.rds'))
    message(paste0(OUTDIR, 'tmp_', sampleNames[sp], '_CG_binned_met_', WDW_SIZE, 'bp_chr', CHR,'.rds'), ' saved')
  }
}

message(Sys.time(), ': Script finished. Input arguments were:\nfile: ', INPUT_QUASR, '\nWDW_SIZE: ', WDW_SIZE, '\nchromosomes: ', CHR1, ' to ', CHR2)


