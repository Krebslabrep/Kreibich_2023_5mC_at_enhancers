#TITLE:   Creating chromHMM annotation tibble for CG bins
#AUTHOR:  Elisa Kreibich
#DATE:    15-04-2022
#AIM:     Creating tibble with chromHMM annotations for each CG bin, including chromHMM clusters.
#OUTPUT:  Tibble of bin GRanges with columns: seqnames, start, end, width, strand, bin, chromHMM, chromHMMcl


# Load libraries and environment-----------------------------------------------------
library(tidyverse)
library(GenomicRanges)
library(plyranges)

WD <- '/g/krebs/kreibich/Kreibich_2023_5mC_at_enhancers/'
source(paste0(WD, 'scripts/utilities/Input_arguments_CA.R')


INDIR   <- 'data/'
OUTDIR  <- 'data/'
OUTNAME <- 'chromHMM_annotations_bins.rds'

# Load data -------------------------------------------------------------------------
## Load bin information
new_WDW_SIZE <- 15
bin_gr <- readRDS(paste0(WD, INDIR, 'CGs_over_baits_', WDW_SIZE,'bp.rds'))
bin_gr <- resize(bin_gr, width = new_WDW_SIZE, fix = "center")
bin_gr$bin <- names(bin_gr)


## Load ChromHMM annotation
### ChromHMM input file gives chromatin annotations for the whole genome segmented into **200bp** windows.
### Data source: previously defined chromHMM states for mouse embryonic stem cells for USCS mm10 (https://github.com/guifengwei/ChromHMM_mESC_mm10) 
### Pintacuda et al., 2017 (doi: 10.1016/j.molcel.2017.11.013)

file <- paste0(WD, INDIR, 'mESC_E14_12_dense.annotated.bed')
chromHMM_anno <- read_bed(file)
chromHMM_anno <- chromHMM_anno[,1]
colnames(elementMetadata(chromHMM_anno)) <- "chromHMM"
names(chromHMM_anno) <- chromHMM_anno$chromHMM

chromHMM_anno_names <- sapply(c(12,2,1,3,5,6,7,9,10,4,8,11), function(x){unique(chromHMM_anno$chromHMM)[grepl(paste0("^", x,"_"), unique(chromHMM_anno$chromHMM))]})
chromHMM_anno_names <- str_replace(chromHMM_anno_names, "6_BivalentChromatin", "6_BivalentPromoter")

# Perform bin annotation -------------------------------------------------------------------------
## Get annotations for all bins
hits <- findOverlaps(bin_gr, chromHMM_anno)
anno <- CharacterList(split(names(chromHMM_anno)[subjectHits(hits)],
                            queryHits(hits)))
bin_gr_anno <- subsetByOverlaps(bin_gr, chromHMM_anno)
mcols(bin_gr_anno) <- DataFrame(mcols(bin_gr_anno), anno)

bin_gr_anno_tbl <- as_tibble(bin_gr_anno) %>%
  mutate_at(vars(anno), as.character)

## Resolve multiple annotations for one bin
bin_ov_count <- bin_gr_anno_tbl %>% 
  mutate(count = str_count(anno, '_')) %>% 
  pull(count) %>% 
  max()

bin_gr_anno_tbl_ov_l <- lapply(seq(bin_ov_count), function(x){
  bin_gr_anno_tbl %>% 
    mutate(anno = str_extract_all(anno, "\\d+_\\w+"),
           anno = map_chr(anno, nth, x)) %>% 
    filter(!is.na(anno))
})

bin_gr_anno_tbl <- bind_rows(bin_gr_anno_tbl_ov_l) %>% 
  dplyr::rename(chromHMM = anno) %>% 
  mutate(chromHMM = str_replace(chromHMM, "6_BivalentChromatin", "6_BivalentPromoter"))

## Assign duplicates to one annotation
anno_newOrder <- c(8, 11, 4, 7, 9, 10, 6, 5, 3, 2, 1, 12)
chromHMM_anno_names2 <- sapply(anno_newOrder, function(x){chromHMM_anno_names[grepl(paste0("^", x,"_"), chromHMM_anno_names)]})

bin_gr_anno_tbl2 <- bin_gr_anno_tbl %>% 
  mutate(chromHMM = factor(chromHMM, level = chromHMM_anno_names2)) %>% 
  arrange(bin, chromHMM) %>% 
  distinct(bin, .keep_all = T)

## Collapse annotations into larger groups
### "1_Insulator"                     --> C1_CTCF
### "2_Intergenic"                    --> C2_Inactive
### "3_Heterochromatin"               --> C2_Inactive
### "4_Enhancer"                      --> C5_CREs_enhancer
### "5_RepressedChromatin"            --> C2_Inactive
### "6_BivalentPromoter"              --> C2_Inactive
### "7_ActivePromoter"                --> C4_CREs_promoter
### "8_StrongEnhancer"                --> C5_CREs_enhancer
### "9_TranscriptionTransition"       --> C3_Transcription
### "10_TranscriptionElongation"      --> C3_Transcription
### "11_WeakEnhancer"                 --> C5_CREs_enhancer
### "12_LowSignal/RepetitiveElements" --> C2_Inactive

chromHMMcl_anno_names <- c("C1_CTCF", "C2_Inactive", "C3_Transcription", "C4_CREs_promoter", "C5_CREs_enhancer")

bin_gr_anno_tbl2 <- bin_gr_anno_tbl2 %>% 
  mutate(chromHMMcl = case_when(
    str_detect(chromHMM, "^1_") ~ chromHMMcl_anno_names[1],
    
    str_detect(chromHMM, "^2_") ~ chromHMMcl_anno_names[2],
    str_detect(chromHMM, "^3_") ~ chromHMMcl_anno_names[2],
    str_detect(chromHMM, "^5_") ~ chromHMMcl_anno_names[2],
    str_detect(chromHMM, "^6_") ~ chromHMMcl_anno_names[2],
    str_detect(chromHMM, "^12_") ~ chromHMMcl_anno_names[2],
    
    str_detect(chromHMM, "^9_") ~ chromHMMcl_anno_names[3],
    str_detect(chromHMM, "^10_") ~ chromHMMcl_anno_names[3],
    
    str_detect(chromHMM, "^7_") ~ chromHMMcl_anno_names[4],
    
    str_detect(chromHMM, "^4_") ~ chromHMMcl_anno_names[5],
    str_detect(chromHMM, "^8_") ~ chromHMMcl_anno_names[5],
    str_detect(chromHMM, "^11_") ~ chromHMMcl_anno_names[5],
  ))

## Check final numbers
bin_gr_anno_tbl2 %>% 
  count(chromHMMcl, chromHMM) 

## Save final tibble
saveRDS(bin_gr_anno_tbl2, paste0(WD, OUTDIR, OUTNAME))

        