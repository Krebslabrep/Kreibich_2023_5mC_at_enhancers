#TITLE:   ChromHMM annotation of TFBS
#AUTHOR:  Elisa Kreibich
#DATE:    07/12/2021
#AIM:     Creating tibble with chromHMM annotations for each TFBS, including chromHMM clusters.
#OUTPUT:  Tibble of bin GRanges with columns: seqnames, start, end, width, strand, TFBS, chromHMM, chromHMMcl

# Load libraries and environment-----------------------------------------------------
library(tidyverse)
library(genomation)
library(GenomicRanges)

WD <- '/g/krebs/kreibich/Kreibich_2023_5mC_at_enhancers/'
source(paste0(WD, 'scripts/utilities/Input_arguments_TFBS.R')

INDIR   <- 'data/'
OUTDIR  <- 'data/'
OUTNAME <- 'chromHMM_annotations_bins.rds'


#GRanges of TFBS ------------------
MotDBf = readRDS(paste0(WD, 'data/GRange_mapped_jaspar2018_ChIP_score10_inBaits_BANP.rds'))
bin_gr <- resize(MotDBf, width = WDW_SIZE, fix = "center")

#ChromHMM annotation -----------------------
input = paste0(WD, INDIR, 'mESC_E14_12_dense.annotated.bed')
genomation::readBed(input, remove.unusual = T) #doesn't work for some reason
skip = detectUCSCheader(input)
df = read_delim(input, skip = skip, n_max = 2, col_names = FALSE, 
                delim = "\t")
numcol = ncol(df)
df = readGeneric(input, skip = skip, keep.all.metadata = TRUE)
chromHMM_anno = df[,"V4"]
colnames(elementMetadata(chromHMM_anno)) <- "ChromHMM_E14_mm10"
names(chromHMM_anno) <- chromHMM_anno$ChromHMM_E14_mm10
saveRDS(chromHMM_anno, paste0(WD, OUTDIR, 'chromHMM_anno.rds')

chromHMM_anno_names <- sapply(seq(12), function(x){unique(chromHMM_anno$ChromHMM_E14_mm10)[grepl(paste0("^", x,"_"), unique(chromHMM_anno$ChromHMM_E14_mm10))]})
chromHMM_anno_names <- str_replace(chromHMM_anno_names, "6_BivalentChromatin", "6_BivalentPromoter")
saveRDS(chromHMM_anno_names, paste0(WD, OUTDIR, 'chromHMM_anno_names.rds'))

chromHMM_anno <- readRDS(paste0(WD, OUTDIR, 'chromHMM_anno.rds'))
chromHMM_anno_names <- readRDS(paste0(WD, OUTDIR, 'chromHMM_anno_names.rds'))

#Get annotations for all bins -----------------
hits <- findOverlaps(bin_gr, chromHMM_anno)
anno <- CharacterList(split(names(chromHMM_anno)[subjectHits(hits)],
                            queryHits(hits)))
bin_gr_anno <- subsetByOverlaps(bin_gr, chromHMM_anno)
mcols(bin_gr_anno) <- DataFrame(mcols(bin_gr_anno), anno)

bin_gr_anno_tbl <- as_tibble(bin_gr_anno) %>%
  mutate_at(vars(anno), as.character)

bin_anno_tbl <- bin_gr_anno_tbl %>% dplyr::select(TFBS, anno)

#Resolve multiple annotations for one bin
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
  mutate(chromHMM = str_replace(chromHMM, "6_BivalentChromatin", "6_BivalentPromoter")) %>% 
  mutate(chromHMM = factor(chromHMM, levels = (chromHMM_anno_names)))

bin_gr_anno_tbl %>% count(TFBS) %>% count(n) #163111 

saveRDS(bin_gr_anno_tbl, paste0(WD, OUTDIR, 'chromHMM_', WDW_SIZE,'bp_TFBS.rds'))


# Assign duplicates to one annotation -----------------------
anno_newOrder <- c(8, 11, 4, 7, 9, 10, 6, 5, 3, 2, 1, 12)
chromHMM_anno_names <- sapply(anno_newOrder, function(x){chromHMM_anno_names[grepl(paste0("^", x,"_"), chromHMM_anno_names)]})

bin_gr_anno_tbl_2 <- bin_gr_anno_tbl %>% 
  mutate(chromHMM = factor(chromHMM, level = chromHMM_anno_names)) %>% 
  arrange(TFBS, chromHMM) %>% 
  distinct(TFBS, .keep_all = T)

#Collapse annotations into larger groups ----------------------------------
# "1_Insulator"                     --> CTCF  
# "2_Intergenic"                    --> Inactive  
# "3_Heterochromatin"               --> Inactive  
# "4_Enhancer"                      --> CREs_enhancer  
# "5_RepressedChromatin"            --> Inactive  
# "6_BivalentPromoter"              --> Inactive  
# "7_ActivePromoter"                --> CREs_promoter  
# "8_StrongEnhancer"                --> CREs_enhancer  
# "9_TranscriptionTransition"       --> Transcription  
# "10_TranscriptionElongation"      --> Transcription  
# "11_WeakEnhancer"                 --> CREs_enhancer  
# "12_LowSignal/RepetitiveElements" --> Inactive  


chromHMMcl_anno_names <- c("C1_CTCF", "C2_Inactive", "C3_Transcription", "C4_CREs_promoter", "C5_CREs_enhancer")

bin_gr_anno_tbl_3 <- bin_gr_anno_tbl_2 %>% 
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

bin_gr_anno_tbl_3 %>% count(chromHMMcl, chromHMM) 

saveRDS(bin_gr_anno_tbl_3, paste0(WD, OUTDIR, 'chromHMM_', WDW_SIZE,'bp_TFBS_with_clusters.rds'))
