#TITLE:   Chromatin Modifications Heatmap
#AUTHOR:  Elisa Kreibich
#DATE:    15-09-2021
#AIM:     Creating Heatmap of different chromatin modifications for antagonist sites

#Set envorinment -------------
WD <- '/g/krebs/kreibich/Kreibich_2023_5mC_at_enhancers/'

#Load libraries --------------
library(tidyverse)
library(GenomicRanges)
library(scales)
library(ggpubr)
library(patchwork)
library(plyranges)
library(rtracklayer)
library(circlize)
library(EnrichedHeatmap)

source(paste0(WD, 'data/COLORS.R'))
DATE <- Sys.Date()

# Region of interest --------------
SMF <- readRDS() #GRange Final data tibble ES_NO (e.g.'./data_results/Final_data_tibble_CA_2022-10-21_SMF_MM_ES_NO_R1_R2_R5a6_curatedICR_cObin10_cOCMH30_states_gr.rds')
ROI <- SMF %>% filter(state2 == "antagonist")

# Histone Modifications --------------
#Load in ChIP-seq data (bigwig file)
HM <- c(  '/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K27Ac/NT/PMID_2763441/signal/H3K27Ac_mESC_NT_SRX112928_rpgc.bigwig',
          '/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K27me3/NT/PMID_28212747/signal/H3K27me3_mESC_NT_SRX1342331_rpgc.bigwig',
          '/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K4me1/NT/PMID_22170606/signal/H3K4me1_mESC_NT_SRX080175_rpgc.bigwig',
          '/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K4me3/NT/PMID_2763441/signal/H3K4me3_mESC_NT_SRX062993_rpgc.bigwig'
        )
HM_names <- lapply(HM, function(x){
  name <- str_remove(x, ".{1,}/signal/") %>% 
    str_remove("_rpgc.bigwig") %>% 
    str_remove("_mESC.{1,}")
})

HM_data <- lapply(HM, import.bw)
names(HM_data) <- HM_names

mat1 = EnrichedHeatmap::normalizeToMatrix(HM_data$H3K27Ac, ROI, extend = 2000, value_column = "score", w = 50)
EnrichedHeatmap(mat1, name = "H3K27Ac")

mat1_l = lapply(HM_data, function(gr){
  EnrichedHeatmap::normalizeToMatrix(gr, ROI, extend = 2000, value_column = "score", w = 50, smooth = TRUE)
})
EnrichedHeatmap(mat1_l[[1]], name = HM_names[[1]], 
                top_annotation = HeatmapAnnotation(enrich = anno_enriched(axis_param = list(side = "left")))) + 
  
  EnrichedHeatmap(mat1_l[[2]], name = HM_names[[2]]) +
  EnrichedHeatmap(mat1_l[[3]], name = HM_names[[3]]) +
  EnrichedHeatmap(mat1_l[[4]], name = HM_names[[4]])
  


# DHS -------------------------
#Load in DNase-seq data (bigwig file)
DHS <- c(
  '/g/krebs/kreibich/HTS/GEO/MM/wigs/DNAse-seq_mESC_WT_rep1_GSM1657364.bw'
)
DHS_names <- "DHS"

DHS_data <- lapply(DHS, import.bw)
names(DHS_data) <- DHS_names

mat1_DHS = EnrichedHeatmap::normalizeToMatrix(DHS_data[[1]], ROI, extend = 2000, value_column = "score", w = 50, smooth = TRUE)
EnrichedHeatmap(mat1_DHS, name = DHS_names)

# WGBS ------------------------
#Load in WGBS data (bigwig file)
WGBS <- c(
  '/g/krebs/kreibich/analysis/methylomes/MM/Context_methylation_call_WGBS_MM_ESC_R1_3.bw'
)
WGBS_names <- "WGBS"

WGBS_data <- lapply(WGBS, import.bw)
names(WGBS_data) <- WGBS_names

mat1_WGBS = EnrichedHeatmap::normalizeToMatrix(WGBS_data[[1]], ROI, extend = 2000, value_column = "score", w = 50, smooth = TRUE)
EnrichedHeatmap(mat1_WGBS, name = WGBS_names)


# CTCF ChIP ------------------------------------------------
#Load in CTCF ChIP-seq data (bigwig file)
CTCF <- c(
  '/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/CTCF/NT/PMID_22170606/signal/CTCF_mESC_NT_SRX080167_rpgc.bigwig'
)
CTCF_names <- lapply(CTCF, function(x){
  str_remove(x, ".{1,}/signal/") %>% 
    str_remove("_rpgc.bigwig") %>% 
  str_remove("_mESC.{1,}")
})

CTCF_data <- lapply(CTCF, import.bw)
names(CTCF_data) <- CTCF_names

mat1_CTCF = EnrichedHeatmap::normalizeToMatrix(CTCF_data[[1]], ROI, extend = 2000, value_column = "score", w = 50, smooth = TRUE)
EnrichedHeatmap(mat1_CTCF, name = CTCF_names[[1]])

# Combine plots ----------------------------
col_fun1 = colorRamp2(quantile(mat1_l[[1]], c(0.0, 0.75, 0.99)), c("blue", "white", "red"))
col_fun2 = colorRamp2(quantile(mat1_l[[2]], c(0.0, 0.75, 0.99)), c("blue", "white", "red"))
col_fun3 = colorRamp2(quantile(mat1_l[[3]], c(0.0, 0.75, 0.99)), c("blue", "white", "red"))
col_fun4 = colorRamp2(quantile(mat1_l[[4]], c(0.0, 0.75, 0.99)), c("blue", "white", "red"))
col_fun5 = colorRamp2(quantile(mat1_DHS, c(0.0, 0.75, 0.99)), c("blue", "white", "red"))
col_fun6 = colorRamp2(quantile(mat1_CTCF, c(0.0, 0.75, 0.99)), c("blue", "white", "red"))

png(file= paste0(WD, "data_results/plots/Heatmap_histoneMods_antagonist_sites", DATE, "_.png"), width=680, height=750)
# pdf(file= paste0(WD, "data_results/plots/Heatmap_histoneMods_antagonist_sites", DATE, "_.pdf")

EnrichedHeatmap(mat1_l[[1]], col = col_fun1, name = HM_names[[1]], column_title =  HM_names[[1]],
                top_annotation = HeatmapAnnotation(enrich = anno_enriched(axis_param = list(side = "left")))) + 
  
  EnrichedHeatmap(mat1_l[[2]], col = col_fun2, name = HM_names[[2]], column_title =  HM_names[[2]]) +
  EnrichedHeatmap(mat1_l[[3]], col = col_fun3, name = HM_names[[3]], column_title =  HM_names[[3]]) +
  EnrichedHeatmap(mat1_l[[4]], col = col_fun4, name = HM_names[[4]], column_title =  HM_names[[4]]) +
  EnrichedHeatmap(mat1_DHS, col = col_fun5, name = DHS_names, column_title =  DHS_names) +
  EnrichedHeatmap(mat1_CTCF, col = col_fun6, name = CTCF_names[[1]], column_title =  CTCF_names[[1]], top_annotation = HeatmapAnnotation(enrich = anno_enriched(axis_param = list(side = "left"))))
dev.off()
