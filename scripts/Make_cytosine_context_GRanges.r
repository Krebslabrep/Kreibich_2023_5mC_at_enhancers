#TITLE:   Creating GRange objects for different cytosine contexts
#AUTHOR:  Elisa Kreibich
#DATE:    15-04-2021
#AIM:     Creating GRange objects for different cytosine contexts (NWCGW for CpG; and DGCHN for GpC) 
#         needed for cytosine annotations in subsequent analyses and scripts.


# Load libraries -----------------------------------------------------
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
library(parallel)
library(tidyverse)

WD <- '/g/krebs/kreibich/Kreibich_2023_5mC_at_enhancers/'
source(paste0(WD, 'scripts/utilities/Input_arguments_CA.R')

# Make GRange objects -----------------------------------------------------
# Potential cytosine contexts: "DGCHN" "CGH"   "NWCGW" "GCH"   "GCG"  
INDIR       <- 'data/'
OUTDIR      <- 'data/'
Cpatterns   <- c("NWCGW", "DGCHN")
musmus_chr  <- paste0('chr', c(seq(19), 'X', 'Y'))

for(pattern in Cpatterns) {
  mm10_Cs <- mclapply(seq_along(musmus_chr), function(i){
    reg <- musmus_chr[i]
    Cs_chr_XSV <- matchPattern(DNAString(pattern), Mmusculus[[reg]], fixed = "subject")
    Cs_chr <- GRanges(seqnames = Rle(as.character(reg)), ranges = IRanges(start(Cs_chr_XSV), end = end(Cs_chr_XSV)), strand = "+")
    return(Cs_chr)  
  }, mc.cores = 4)
  mm10_Cs <- unlist(GRangesList(mm10_Cs))
  mm10_Cs <- resize(mm10_Cs, 1, fix = "center")
  mm10_Cs$type <- pattern
  mm10_Cs$cytosine <- paste0(pattern, "_", seq_along(mm10_Cs))
  names(elementMetadata(mm10_Cs)) <- c("type", pattern)
  
  saveRDS(mm10_Cs, paste0(OUTDIR, pattern,'_GRange.rds'))
}

# NWCGW_gr <- readRDS(paste0(OUTDIR, Cpatterns[1],'_GRange.rds'))
# DGCHN_gr <- readRDS(paste0(OUTDIR, Cpatterns[2],'_GRange.rds'))



# Combine NWCGW GRange with CG bin information ---------------------------------------------------------
## Load bin information
bin_gr <- readRDS(paste0(INDIR, 'CGs_over_baits_', WDW_SIZE,'bp.rds'))
bin_gr <- resize(bin_gr, width = 1, fix = "center")
bin_gr$bin <- names(bin_gr)
  
  
## Reload NWCGW info
NWCGW_gr <- readRDS(paste0(OUTDIR, Cpatterns[1],'_GRange.rds'))
names(NWCGW_gr) <- NWCGW_gr$NWCGW

## Get all NWCGW within bins and add NWCGW information
hits <- findOverlaps(bin_gr, NWCGW_gr)
NWCGW <- CharacterList(split(names(NWCGW_gr)[subjectHits(hits)],
                           queryHits(hits)))
NWCGW_bins <- subsetByOverlaps(bin_gr, NWCGW_gr)
mcols(NWCGW_bins) <- DataFrame(mcols(NWCGW_bins), NWCGW)

NWCGW_bins_tbl <- as_tibble(NWCGW_bins) %>%
  mutate_at(vars(NWCGW), as.character)

## Check for multiple NWCGW annotation per bin
bin_ov_count <- NWCGW_bins_tbl %>% 
  mutate(count = str_count(NWCGW, 'NWCGW_\\d+')) %>% 
  pull(count) %>% 
  max()
if(bin_ov_count > 1) {message("Multiple NWCGWs annotated per bin. >>> Check bin size.")}


## Create final tibble & GRange object
# NWCGW_bin_gr <- as_granges(NWCGW_bins_tbl)

saveRDS(NWCGW_bins_tbl, paste0(OUTDIR, 'CGs_over_baits_', WDW_SIZE,'bp_', Cpatterns[1], '.rds'))

