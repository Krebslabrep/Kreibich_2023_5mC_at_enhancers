library(tidyverse)
library(GenomicRanges)
library(scales)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(plyranges)
library(rtracklayer)
library(parallel)
library(ggbio)
library(Mus.musculus)
library(QuasR)
library(SingleMoleculeFootprinting)


source('/g/krebs/kreibich/analysis/SMF/manuscript/COLORS.R')
PLOTOUT <- paste('/g/krebs/kreibich/analysis/SMF/manuscript/')
DATE <- Sys.Date()

#Plotting functions --------------
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}
plot_enhancer <- function(tbl, COLOR){
  tbl %>% 
    # select(start, end) %>% 
    # mutate(ID = paste0("enhancer_", seq(1:nrow(.)))) %>% 
    # pivot_longer(cols = -ID, names_to = "position", values_to = "seq") %>% 
    # mutate(target = "enhancer") %>% 
    ggplot(aes(seq, 1, fill = state2)) +
    geom_area(show.legend = F) +
    scale_fill_manual(values = COLOR) +
    xlim(start(ROI), end(ROI)) +
    facet_wrap(~target, ncol = 1, strip.position = "right") +
    theme_EK() +
    theme(strip.background = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          plot.background = element_blank(), 
          panel.background = element_blank()
    )
}

plot_signal_area <- function(tbl, TARGET, MAX, COLOR){
  tbl %>% 
    filter(target == TARGET) %>% 
    ggplot(aes(start, score, fill = target)) +
    geom_area() +
    scale_fill_manual(values = COLOR) +
    scale_y_continuous(limits = c(0, MAX), breaks = c(0, MAX)) +
    scale_x_continuous(expand = expansion(mult = 0, add = 0), limits = c(start(ROI), end(ROI))) +
    labs(x = seqnames(ROI)) +
    facet_wrap(~target, ncol = 1, strip.position = "right") +
    theme_EK() +
    theme(strip.background = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          plot.background = element_blank(), 
          panel.background = element_blank()
          ) 
}

plot_signal_area2 <- function(tbl, TARGET, MAX, COLOR){
  tbl %>% 
    filter(target == TARGET) %>% 
    ggplot(aes(start, score2, fill = target)) +
    geom_area() +
    scale_fill_manual(values = COLOR) +
    scale_y_continuous(limits = c(0, MAX), breaks = c(0, MAX)) +
    scale_x_continuous(expand = expansion(mult = 0, add = 0), limits = c(start(ROI), end(ROI))) +
    labs(x = seqnames(ROI)) +
    # facet_wrap(~target, ncol = 1, strip.position = "right") +
    theme_EK() +
    theme(strip.background = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          plot.background = element_blank(), 
          panel.background = element_blank()) 
}

plot_signal_area3 <- function(tbl, TARGET, MAX, COLOR){
  tbl %>% 
    filter(target == TARGET) %>% 
    ggplot(aes(start, score2, fill = target)) +
    geom_area() +
    scale_fill_manual(values = COLOR) +
    scale_y_continuous(limits = c(0, MAX), breaks = c(0, MAX), expand = expansion(mult = 0, add = 0)) +
    scale_x_continuous(expand = expansion(mult = 0, add = 0), limits = c(start(ROI), end(ROI))) +
    labs(x = seqnames(ROI)) +
    facet_wrap(~target, ncol = 1, strip.position = "right") +
    theme_EK() +
    theme(strip.background = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          panel.border = element_blank(),
          axis.line.y = element_line(),
          plot.background = element_blank(), 
          panel.background = element_blank()) 
}

plot_signal_point <- function(tbl, TARGET, COLOR){
  tbl %>% 
    filter(target == TARGET) %>% 
    ggplot(aes(start, score, color = target)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = COLOR) +
    scale_y_continuous(limits = c(0,1), expand = expansion(add = 0.1), breaks = c(0,1)) +
    scale_x_continuous(expand = expansion(mult = 0, add = 0), limits = c(start(ROI), end(ROI))) +
    labs(x = seqnames(ROI)) +
    facet_wrap(~target, ncol = 1, strip.position = "right") +
    theme_EK() +
    theme(strip.background = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          plot.background = element_blank(), 
          panel.background = element_blank()) 
}

plot_signal_point2 <- function(tbl, TARGET, COLOR){
  tbl %>% 
    filter(target == TARGET) %>% 
    ggplot(aes(start, score, color = target)) +
    geom_point(size = 1) +
    scale_color_manual(values = COLOR) +
    scale_y_continuous(limits = c(0,1), expand = expansion(mult = 0, add = 0.1), breaks = c(0,1)) +
    scale_x_continuous(expand = expansion(mult = 0, add = 0), limits = c(start(ROI), end(ROI))) +
    labs(x = seqnames(ROI)) +
    facet_wrap(~target, ncol = 1, strip.position = "right") +
    theme_EK() +
    theme(strip.background = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          plot.background = element_blank(), 
          panel.border = element_blank(),
          axis.line.y = element_line(),
          axis.line.x = element_line(),
          panel.background = element_blank()) 
}

plot_signal_point3 <- function(tbl, TARGET, COLOR){
  tbl %>% 
    filter(target == TARGET) %>%
    ggplot(aes(start, score, color = target, size = 1)) +
    geom_line() +
    geom_point(data = (. %>% filter(dplyr::between(score, 0.1, 0.6))), size = 1) +
    scale_color_manual(values = COLOR) +
    scale_y_continuous(limits = c(0,1), expand = expansion(mult = 0, add = 0.1), breaks = c(0,1)) +
    scale_x_continuous(expand = expansion(mult = 0, add = 0), limits = c(start(ROI), end(ROI))) +
    labs(x = seqnames(ROI)) +
    facet_wrap(~target, ncol = 1, strip.position = "right") +
    theme_EK() +
    theme(strip.background = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          plot.background = element_blank(), 
          panel.border = element_blank(),
          axis.line.y = element_line(),
          panel.background = element_blank()) 
}

#Plot signal point with enhancer points colored
plot_signal_point4 <- function(tbl, TARGET, COLOR, COLOR_EOI){
  COLOR2 <- c(COLOR, COLOR_EOI)
  names(COLOR2) <- c(FALSE, TRUE)
  SIZE <- c(1, 1)
  names(SIZE) <- c(FALSE, TRUE)
  
  tbl %>% 
    filter(target == TARGET) %>% 
    ggplot(aes(start, score, color = EOI, size = EOI)) +
    geom_point() +
    scale_color_manual(values = COLOR2) +
    scale_size_manual(values = SIZE) +
    scale_y_continuous(limits = c(0,1), expand = expansion(mult = 0, add = 0.15), breaks = c(0,1)) +
    scale_x_continuous(expand = expansion(mult = 0, add = 0), limits = c(start(ROI), end(ROI))) +
    labs(x = seqnames(ROI)) +
    facet_wrap(~target, ncol = 1, strip.position = "right") +
    theme_EK() +
    theme(strip.background = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          plot.background = element_blank(), 
          panel.border = element_blank(),
          axis.line.y = element_line(),
          panel.background = element_blank()) 
}

# Plotting version  --------------
# VERSION <- "h5mC"
VERSION <- "F1"

# Load SMF Data --------------------------------------------------------
SAMPLES   <- c("ES_NO", "TKO_DE", "TETTKO_NO")
TKO_TREAT <- "DE"
TET_TREAT <- "NO"

Qinput    <- '/g/krebs/kreibich/HTS/SMF/MM/QuasR/QuasR_rmdup/QuasR_aligned_files_ALL_rmdup.txt'

MySamples_1 <- c("ES_NO_R1", "ES_NO_R2", "ES_NO_R5a6")
print(MySamples_1)

if(TKO_TREAT == "DE"){
  MySamples_2 <- c("TKO_DE_R1a2", "TKO_DE_R3a4", "TKO_DE_R5a6")
} else if (TKO_TREAT == "NO"){
  MySamples_2 <- c("TKO_NO_R3", "TKO_NO_R4", "TKO_NO_R3")
}
print(MySamples_2)

MySamples_3 <- c("TETTKO_NO_R1", "TETTKO_NO_R2", "TETTKO_NO_R1")

## BAMs  
QuasRprj <- qAlign(sampleFile = Qinput, 
                   genome = "BSgenome.Mmusculus.UCSC.mm10",
                   paired = "fr",
                   bisulfite = "undir")
QuasRprj@aligner <- "Rbowtie"
Samples <- QuasR::alignments(QuasRprj)[[1]]$SampleName

### ES
Samples_1 <- Samples[grepl(paste(MySamples_1, collapse = "|"), Samples)]
print(Samples_1)
### TKO
Samples_2 <- Samples[grepl(paste(MySamples_2, collapse = "|"), Samples)]
print(Samples_2)
### TETTKO
Samples_3 <- Samples[grepl(paste(MySamples_3, collapse = "|"), Samples)]
print(Samples_3)


# Region of interest --------------
ICRs_curated <- readRDS('/g/krebs/kreibich/analysis/SMF/revision/R1_3_ICRs/Data_tbl_ICRs_curated_3_2022-10-21.rds')
ICRs_curated_names <- ICRs_curated %>% pull(names) %>% unique() %>% as.character()


ICR_gr <- readRDS('/g/krebs/kreibich/GitHub_Kreibich_et_al/rds/Impr_regions_NWCGW_all_info_gr.rds') %>% as_granges()
# ICR_gr <- ICR_gr[!grepl("meth", ICR_gr$names)]
ICR_gr <- ICR_gr[ICR_gr$names %in% ICRs_curated_names]

# C_ICR_gr <- readRDS('/g/krebs/kreibich/GitHub_Kreibich_et_al/rds/Impr_regions_NWCGW_gr.rds')
# ROIs <- ICR_gr[ICR_gr$names == "Dlk1-Gtl2 IG"]
# ROIs <- resize(ROIs, width = 6000, fix = "center")
# ROI <- resize(ROI, width = WIDTH, fix = "start")

for (x in seq_along(ICR_gr)){
  # x=1
  # print(x)
# WIDTH_l <- seq(8000, 9000, 500)
# for (x in seq_along(WIDTH_l)){
  # WIDTH = WIDTH_l[x]
# print(x)
# x=1

ROI_ori <- ICR_gr[x]
WIDTH_ORI <- width(ROI_ori)
WIDTH_ORI_top <- round_any(WIDTH_ORI, 500)
WIDTH_l <- seq(WIDTH_ORI_top, WIDTH_ORI_top+3000, 1000)

for (WIDTH in WIDTH_l){
ROI <- resize(ROI_ori, width = WIDTH, fix = "center")
print(ROI)

# NWCGWs in region of interest --------------
if(exists("data_gr") == FALSE){
  data_gr <- readRDS(paste0('/g/krebs/kreibich/GitHub_Kreibich_et_al/single_molecule_call/rds/Final_data_tibble_CA_SMF_MM_ES_NO_R1_R2_R5a6_cO30_wCMH_2022-06-23_states_gr.rds'))
}

EOI_all <- data_gr[data_gr$ICR == TRUE]
EOI     <- subsetByOverlaps(EOI_all, ROI) 
EOI     <- resize(EOI, width = 1, fix = "center")
EOI_c   <- resize(EOI, width = 1, fix = "center")
if(length(EOI) == 0) {next}

# Histone Modifications ---------------------------------------
COLORS_HM <- c("#000000",'#d7b5d8','#dd1c77', '#980043')
names(COLORS_HM) <- c("H3K27me3", "H3K4me3","H3K27Ac", "H3K4me1")


HM <- c(
          # '/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K27Ac/DNMT_WT/PMID_26675734/signal/H3K27Ac_mESC_DNMT_WT_SRX1280444_rpgc.bigwig',
           # '/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K27Ac/DNMT_WT/PMID_26675734/signal/H3K27Ac_mESC_DNMT_WT_SRX1280445_rpgc.bigwig',
           '/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K27Ac/NT/PMID_2763441/signal/H3K27Ac_mESC_NT_SRX112928_rpgc.bigwig',
           '/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K27me3/NT/PMID_28212747/signal/H3K27me3_mESC_NT_SRX1342331_rpgc.bigwig',
           # '/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K27me3/NT/PMID_22179133/signal/H3K27me3_mESC_NT_SRX032466_rpgc.bigwig',
           # '/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K4me1/NT/PMID_2763441/signal/H3K4me1_mESC_NT_SRX062992_rpgc.bigwig',
           '/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K4me1/NT/PMID_22170606/signal/H3K4me1_mESC_NT_SRX080175_rpgc.bigwig',
           '/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K4me3/NT/PMID_2763441/signal/H3K4me3_mESC_NT_SRX062993_rpgc.bigwig'
           # '/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K4me3/NT/PMID_28212747/signal/H3K4me3_mESC_NT_SRX1342333_rpgc.bigwig'
        )
HM_names <- lapply(HM, function(x){
 str_remove(x, ".{1,}/signal/") %>% 
    str_remove("_rpgc.bigwig") %>% 
    str_remove("_mESC.{1,}")
})

if(exists("HM_data") == FALSE){
  HM_data <- lapply(HM, import.bw)
  names(HM_data) <- HM_names
}

HM_data_ROI <- lapply(HM_data, subsetByOverlaps, ROI)
HM_data_ROI <- lapply(HM_data_ROI, as_tibble)
HM_data_ROI <- lapply(HM_data_ROI, function(x){
  x %>% 
    mutate(score2 = caTools::runmean(score, k = 2)) 
})

HM_data_ROI_tbl <- HM_data_ROI %>% 
  bind_rows(., .id = "target") %>% 
  # mutate(target = str_remove(name, "_mESC.{1,}")) %>% 
  mutate(target = factor(target, levels = c("H3K27me3", "H3K4me3","H3K27Ac", "H3K4me1")))


p_H1 <- plot_signal_area3(HM_data_ROI_tbl, "H3K27me3", 180, COLORS_HM)
p_H2 <- plot_signal_area3(HM_data_ROI_tbl, "H3K4me3", 180, COLORS_HM)
p_H3 <- plot_signal_area3(HM_data_ROI_tbl, "H3K27Ac", 180, COLORS_HM)
p_H4 <- plot_signal_area3(HM_data_ROI_tbl, "H3K4me1", 50, COLORS_HM)
# p_H1 + p_H2 + p_H3 + p_H4 + plot_layout(ncol = 1, heights = c(rep(1,4)))
message("Done: Histone Mods")

# DHS --------------------------------------------------
COLOR_DHS <- "grey30"
names(COLOR_DHS) <- "DHS"

DHS <- c(
  '/g/krebs/kreibich/HTS/GEO/MM/wigs/DNAse-seq_mESC_WT_rep1_GSM1657364.bw'
)
DHS_name <- "DHS"

if(exists("DHS_data") == FALSE){
  DHS_data <- lapply(DHS, import.bw)
  names(DHS_data) <- DHS_name
}

DHS_data_ROI <- lapply(DHS_data, subsetByOverlaps, ROI)
DHS_data_ROI <- lapply(DHS_data_ROI, as_tibble)
DHS_data_ROI <- lapply(DHS_data_ROI, function(x){
  x %>% 
    mutate(score2 = caTools::runmean(score, k = 2)) 
})
DHS_data_ROI_tbl <- DHS_data_ROI %>% 
  bind_rows(., .id = "target")

p_D <- plot_signal_area3(DHS_data_ROI_tbl, TARGET="DHS", MAX=2, COLOR=COLOR_DHS)
message("Done: DHS")

# RNAseq --------------------------------------------------
# COLOR_RNA <- rep(c("black", "goldenrod4"),1)
# # names(COLOR_RNA) <- c("RNAseq WT 1", "RNAseq TET TKO 1", "RNAseq WT 2", "RNAseq TET TKO 2", "RNAseq WT 3", "RNAseq TET TKO 3")
# names(COLOR_RNA) <- c("RNAseq WT 1", "RNAseq TET TKO 1")
# 
# samples_index <- read_csv('/g/krebs/barzaghi/analyses/nf-core_runs/181021_Elisa_rnaseq/sampleSheet.csv')
# 
# RNA <- c(
#   '/g/krebs/barzaghi/analyses/nf-core_runs/181021_Elisa_rnaseq_stranded/results/star_salmon/bigwig/RNA_WT1_rep1.reverse.bigWig',
#   '/g/krebs/barzaghi/analyses/nf-core_runs/181021_Elisa_rnaseq_stranded/results/star_salmon/bigwig/RNA_TETTKO_rep1.reverse.bigWig'
#   # '/g/krebs/barzaghi/analyses/nf-core_runs/181021_Elisa_rnaseq_stranded/results/star_salmon/bigwig/RNA_WT2_rep1.reverse.bigWig',
#   # '/g/krebs/barzaghi/analyses/nf-core_runs/181021_Elisa_rnaseq_stranded/results/star_salmon/bigwig/RNA_TETTKO_rep2.reverse.bigWig',
#   # '/g/krebs/barzaghi/analyses/nf-core_runs/181021_Elisa_rnaseq_stranded/results/star_salmon/bigwig/RNA_WT3_rep1.reverse.bigWig',
#   # '/g/krebs/barzaghi/analyses/nf-core_runs/181021_Elisa_rnaseq_stranded/results/star_salmon/bigwig/RNA_TETTKO_rep3.reverse.bigWig'
# )
# RNA_names <- lapply(RNA, function(x){
#   name <- str_remove(x, ".{1,}/bigwig/") %>%
#     str_remove(".bigWig") 
# })
# # RNA_names <- "RNAseq"
# 
# # RNA_data <- lapply(RNA, import.bw)
# # names(RNA_data) <- RNA_names
# 
# RNA_data_ROI <- lapply(RNA_data, subsetByOverlaps, ROI.1)
# RNA_data_ROI <- lapply(RNA_data_ROI, as_tibble)
# RNA_data_ROI <- lapply(RNA_data_ROI, function(x){
#   x %>% 
#     mutate(score2 = caTools::runmean(score, k = 15)) 
# })
# RNA_data_ROI_tbl <- RNA_data_ROI %>% 
#   bind_rows(., .id = "target") %>% 
#   mutate(target = factor(target, levels = names(RNA_data), labels = names(COLOR_RNA)))
# 
# 
# height.i <- 100
# p_R1 <- plot_signal_area3(RNA_data_ROI_tbl, "RNAseq WT 1", height.i, COLOR_RNA)
# p_R2 <- plot_signal_area3(RNA_data_ROI_tbl, "RNAseq TET TKO 1", height.i, COLOR_RNA)
# # p_R1 + p_R2 + plot_layout(ncol = 1, heights = c(rep(1,2)))
# 
# # p_R3 <- plot_signal_area3(RNA_data_ROI_tbl, "RNAseq WT 2", height.i, COLOR_RNA)
# # p_R4 <- plot_signal_area3(RNA_data_ROI_tbl, "RNAseq TET TKO 2", height.i, COLOR_RNA)
# # p_R5 <- plot_signal_area3(RNA_data_ROI_tbl, "RNAseq WT 3", height.i, COLOR_RNA)
# # p_R6 <- plot_signal_area3(RNA_data_ROI_tbl, "RNAseq TET TKO 3", height.i, COLOR_RNA)
# 
# # p_R1 + p_R2 + p_R3 + p_R4 + p_R5 + p_R6 + plot_layout(ncol = 1, heights = c(rep(1,6)))

# WGBS ------------------------------------------------
COLOR_WGBS <- COLORS_METH[1]
names(COLOR_WGBS) <- "WGBS"

WGBS <- c(
  '/g/krebs/kreibich/analysis/methylomes/MM/Context_methylation_call_WGBS_MM_ESC_R1_3.bw'
)
WGBS_names <- "WGBS"

if(exists("WGBS_data") == FALSE){
  WGBS_data <- lapply(WGBS, import.bw)
  names(WGBS_data) <- WGBS_names
}


WGBS_data_ROI <- lapply(WGBS_data, subsetByOverlaps, ROI)
WGBS_data_ROI <- lapply(WGBS_data_ROI, as_tibble)

WGBS_data_ROI_tbl <- WGBS_data_ROI %>%
  bind_rows(., .id = "target")

p_W <- plot_signal_point2(WGBS_data_ROI_tbl, "WGBS", COLOR_WGBS)
# p_W <- plot_signal_point3(WGBS_data_ROI_tbl, "WGBS", COLOR_WGBS)
message("Done: WGBS")

# WGBS - enhancer ------------------------------------------------
COLOR_WGBS_EOI <- "darkred"
names(COLOR_WGBS_EOI) <- "WGBS|EOI"

# WGBS <- c(
#   '/g/krebs/kreibich/analysis/methylomes/MM/Context_methylation_call_WGBS_MM_ESC_R1_3.bw'
# )
# WGBS_names <- "WGBS"
# 
# WGBS_data <- lapply(WGBS, import.bw)
# names(WGBS_data) <- WGBS_names

WGBS_data_EOI <- lapply(WGBS_data, subsetByOverlaps, EOI)
WGBS_data_EOI <- lapply(WGBS_data_EOI, as_tibble)

WGBS_data_EOI_tbl <- WGBS_data_EOI %>% 
  bind_rows(., .id = "target") %>% 
  mutate(EOI = TRUE)

WGBS_data_tbl2 <- WGBS_data_ROI_tbl %>% 
  left_join(., WGBS_data_EOI_tbl) %>% 
  mutate(EOI = case_when(EOI == TRUE ~ EOI, TRUE~FALSE))

p_W_EOI <- plot_signal_point4(WGBS_data_tbl2, "WGBS", COLOR_WGBS, COLOR_WGBS_EOI)
# p_W_EOI
# p_W <- plot_signal_point3(WGBS_data_ROI_tbl, "WGBS", COLOR_WGBS)

message("Done: WGBS-covered")


# WGBS - F1 ------------------------------------------------
COLOR_WGBS_F1_R <- COLORS_METH[1]
names(COLOR_WGBS_F1_R) <- "5mC_F1_R"

COLOR_WGBS_F1_A <- COLORS_METH[1]
names(COLOR_WGBS_F1_A) <- "5mC_F1_A"

WGBS_F1 <- c(
  '/g/krebs/kreibich/analysis/SMF/methCall/new/ES_F1_CAST_129_NO_ALL_av_meth_context_NArm_Co1_NWCGW_gr.rds'
)
WGBS_F1_names <- "5mC_F1"

if(exists("WGBS_F1_data") == FALSE){
  WGBS_F1_data <- lapply(WGBS_F1, readRDS)
  names(WGBS_F1_data) <- WGBS_F1_names
}

WGBS_F1_data_EOI <- lapply(WGBS_F1_data, subsetByOverlaps, ROI)
WGBS_F1_data_EOI <- lapply(WGBS_F1_data_EOI, as_tibble)

WGBS_F1_data_EOI_tbl <- WGBS_F1_data_EOI %>% 
  bind_rows(., .id = "target") %>% 
  # dplyr::rename(score = av_me) %>% 
  group_by(cytosine, allel) %>% 
  mutate(score = mean(av_me)) %>% 
  ungroup() %>% 
  dplyr::select(-c(replicate, av_me)) %>% 
  distinct(cytosine, allel, .keep_all = T) %>% 
  mutate(target = paste0(target, "_", allel))

p_W_F1_R <- WGBS_F1_data_EOI_tbl %>% 
  filter(allel == "R") %>% 
  plot_signal_point2(., "5mC_F1_R", COLOR_WGBS_F1_R) +
  theme(strip.text = element_text(angle = -90, size = 7, vjust = 1))

p_W_F1_A <- WGBS_F1_data_EOI_tbl %>% 
  filter(allel == "A") %>% 
  plot_signal_point2(., "5mC_F1_A", COLOR_WGBS_F1_A) +
  theme(strip.text = element_text(angle = -90, size = 7, vjust = 1))

message("Done: WGBS SMF F1")
# p_W_F1_R / p_W_F1_A

# Accessibility - F1 ------------------------------------------------
COLOR_SMF_F1_R <- "black"
names(COLOR_SMF_F1_R) <- "SMF_F1_R"

COLOR_SMF_F1_A <- "black"
names(COLOR_SMF_F1_A) <- "SMF_F1_A"

SMF_F1 <- c(
  '/g/krebs/kreibich/analysis/SMF/methCall/new/ES_F1_CAST_129_NO_ALL_av_meth_context_NArm_Co1_gr.rds'
)
SMF_F1_names <- "SMF_F1"

if(exists("SMF_F1_data") == FALSE){
  SMF_F1_data <- lapply(SMF_F1, readRDS)
  names(SMF_F1_data) <- SMF_F1_names
}

SMF_F1_data_EOI <- lapply(SMF_F1_data, subsetByOverlaps, ROI)
SMF_F1_data_EOI <- lapply(SMF_F1_data_EOI, as_tibble)

SMF_F1_data_EOI_tbl <- SMF_F1_data_EOI %>% 
  bind_rows(., .id = "target") %>% 
  filter(context == "DGCHN") %>% 
  # dplyr::rename(score = av_me) %>% 
  group_by(cytosine, allel) %>% 
  mutate(score = 1-mean(av_me)) %>% 
  ungroup() %>% 
  dplyr::select(-c(replicate, av_me)) %>% 
  distinct(cytosine, allel, .keep_all = T) %>% 
  mutate(target = paste0(target, "_", allel))

p_SMF_F1_R <- SMF_F1_data_EOI_tbl %>% 
  filter(allel == "R") %>% 
  plot_signal_point2(., "SMF_F1_R", COLOR_SMF_F1_R) +
  theme(strip.text = element_text(angle = -90, size = 7, vjust = 1))

p_SMF_F1_A <- SMF_F1_data_EOI_tbl %>% 
  filter(allel == "A") %>% 
  plot_signal_point2(., "SMF_F1_A", COLOR_SMF_F1_A) +
  theme(strip.text = element_text(angle = -90, size = 7, vjust = 1))


# p_SMF_F1_R / p_SMF_F1_A
message("Done: Accessibility SMF F1")

# 5hmC ------------------------------------------------
COLOR_hmC <- COLORS_METH[3]
names(COLOR_hmC) <- "5hmC"

hmC <- c(
  '/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/ESC/5hmC/IP/PMID_21514197/signal/5hmC__ESC_SRX057744_rpgc.bigwig'
)
hmC_names <- "5hmC"

if(exists("hmC_data") == FALSE){
  hmC_data <- lapply(hmC, import.bw)
  names(hmC_data) <- hmC_names
}

# resROI <- resize(ROI, width = width(ROI), fix = "center")
hmC_data_ROI <- lapply(hmC_data, subsetByOverlaps, ROI)
hmC_data_ROI <- lapply(hmC_data_ROI, as_tibble)
hmC_data_ROI <- lapply(hmC_data_ROI, function(x){
  x %>% 
    mutate(score2 = caTools::runmean(score, k = 1)) 
})

hmC_data_ROI_tbl <- hmC_data_ROI %>% 
  bind_rows(., .id = "target")

p_5hmC <- plot_signal_area3(hmC_data_ROI_tbl, "5hmC", 30, COLOR_hmC)
  # scale_x_continuous(expand = expansion(mult = 0, add = 0), limits = c(start(resROI), end(resROI)))

p_5hmC

message("Done: 5hmC")

# SMF WT ------------------------------------------------
COLOR_SMF_WT <- c("black", COLORS_METH[3])
names(COLOR_SMF_WT) <- c("SMF_WT_CA","SMF_WT_5mC")
COV           <- 10
SAMPLE        <- Samples_1

MethGR_1_l <- lapply(seq_along(SAMPLE), function(x){
  # x=2
  MethGR  <- try(CallContextMethylation(sampleSheet = Qinput, sample = SAMPLE[x], genome = BSgenome.Mmusculus.UCSC.mm10, coverage = COV, range = ROI))
  if(!("try-error" %in% class(MethGR))) {
    MethGR  <- MethGR[[1]]
    MethGR2 <- lapply(MethGR, as_tibble)
    MethGR2 <- bind_rows(MethGR2, .id = "context") %>% 
      dplyr::rename(meth = 7)
  } else {
    MethGR2 <- NULL
  }
  return(MethGR2)
})
names(MethGR_1_l) <- SAMPLE

SMF_WT_data_ROI_tbl <- bind_rows(MethGR_1_l, .id = "sample") %>% 
  group_by(seqnames, start, end) %>% 
  mutate(score = mean(meth, na.rm = TRUE)) %>% 
  ungroup() %>% 
  dplyr::select(-c(sample, meth)) %>% 
  distinct(seqnames, start, .keep_all = TRUE) %>% 
  mutate(target = "SMF_WT", 
         target = case_when(
           context == "DGCHN" ~ paste0(target, "_CA"),
           context == "NWCGW" ~ paste0(target, "_5mC")
         ))

p_SMF_WT_CA <- SMF_WT_data_ROI_tbl %>% 
  filter(target == "SMF_WT_CA") %>% 
  mutate(score = 1-score) %>% 
  plot_signal_point2(., "SMF_WT_CA", COLOR_SMF_WT) +
  geom_line()+
  theme(strip.text = element_text(angle = -90, size = 7, vjust = 1))

p_SMF_WT_5mC <- SMF_WT_data_ROI_tbl %>% 
  filter(target == "SMF_WT_5mC") %>% 
  plot_signal_point2(., "SMF_WT_5mC", COLOR_SMF_WT) +
  theme(strip.text = element_text(angle = -90, size = 7, vjust = 1))

p_SMF_WT_CA + p_SMF_WT_5mC + plot_layout(ncol=1)

message("Done: SMF_WT")

# SMF TET ------------------------------------------------
COLOR_SMF_TET <- c("goldenrod", "turquoise3")
names(COLOR_SMF_TET) <- c("SMF_TET_CA","SMF_TET_5mC")
COV           <- 10
SAMPLE        <- Samples_3

MethGR_3_l <- lapply(seq_along(SAMPLE), function(x){
  # x=2
  MethGR  <- try(CallContextMethylation(sampleSheet = Qinput, sample = SAMPLE[x], genome = BSgenome.Mmusculus.UCSC.mm10, coverage = COV, range = ROI))
  if(!("try-error" %in% class(MethGR))) {
    MethGR  <- MethGR[[1]]
    MethGR2 <- lapply(MethGR, as_tibble)
    MethGR2 <- bind_rows(MethGR2, .id = "context") %>% 
      dplyr::rename(meth = 7)
  } else {
    MethGR2 <- NULL
  }
  return(MethGR2)
})
names(MethGR_3_l) <- SAMPLE

SMF_TET_data_ROI_tbl <- bind_rows(MethGR_3_l, .id = "sample") %>% 
  group_by(seqnames, start, end) %>% 
  mutate(score = mean(meth, na.rm = TRUE)) %>% 
  ungroup() %>% 
  dplyr::select(-c(sample, meth)) %>% 
  distinct(seqnames, start, .keep_all = TRUE) %>% 
  mutate(target = "SMF_TET", 
         target = case_when(
           context == "DGCHN" ~ paste0(target, "_CA"),
           context == "NWCGW" ~ paste0(target, "_5mC")
         ))

p_SMF_TET_CA <- SMF_TET_data_ROI_tbl %>% 
  filter(target == "SMF_TET_CA") %>% 
  mutate(score = 1-score) %>% 
  plot_signal_point2(., "SMF_TET_CA", COLOR_SMF_TET) +
  geom_line()+
  theme(strip.text = element_text(angle = -90, size = 7, vjust = 1))

p_SMF_TET_5mC <- SMF_TET_data_ROI_tbl %>% 
  filter(target == "SMF_TET_5mC") %>% 
  plot_signal_point2(., "SMF_TET_5mC", COLOR_SMF_TET) +
  theme(strip.text = element_text(angle = -90, size = 7, vjust = 1))

p_SMF_TET_CA + p_SMF_TET_5mC + plot_layout(ncol=1)

message("Done: SMF_TET")

# cCREs ------------------------------------------------
#candidate CREs data from ENCODE

COLOR_cCREs <- "orange"
names(COLOR_cCREs) <- "cCREs"

cCREs <- c(
  '/g/krebs/kreibich/DB/ENCFF693RIQ.bed'
)
cCREs_names <- "cCREs"

if(exists("cCREs_data") == FALSE){
  cCREs_data <- read.delim(file = cCREs, header = F, 
                           col.names = c("seqnames", "start", "end", "name", "x0", "x1", "x2", "x3", "x4", "info", "info2")) %>% 
    as_tibble() %>% 
    dplyr::select(-starts_with("x")) %>% 
    as_granges() %>% 
    list()
  names(cCREs_data) <- cCREs_names
}

cCREs_data_ROI <- lapply(cCREs_data, subsetByOverlaps, ROI)
cCREs_data_ROI <- lapply(cCREs_data_ROI, as_tibble)

cCREs_data_ROI_tbl <- cCREs_data_ROI %>% 
  bind_rows(., .id = "target")

p_cCRES <- cCREs_data_ROI_tbl %>%
  ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = target)) +
  geom_rect(show.legend = F) +
  # geom_label_repel(aes(x = (start+end)/2, y = 0.5, label = info), min.segment.length = 0, color = "black", fill = "white", size = 2) +
  geom_label(aes(x = (start+end)/2, y = 0.5, label = info), color = "black", fill = "white", size = 2) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 1), expand = expansion(mult = 0, add = 0)) +
  scale_x_continuous(expand = expansion(mult = 0, add = 0), limits = c(start(ROI), end(ROI)), breaks = c(start(ROI), end(ROI))) +
  scale_fill_manual(values = COLOR_cCREs) +
  theme_void() 

p_cCRES <- p_cCRES +
  geom_vline(xintercept = start(ROI_ori), color = "grey30", linetype = 3, lwd = 1) +
  geom_vline(xintercept = end(ROI_ori), color = "grey30", linetype = 3, lwd = 1) 


message("Done: candidate CREs")

# Gene annotation ----------------------------------------------------------
#load gene symbol : GRanges, one gene/row
data(genesymbol, package = "biovizBase")
#Plot the different transcripts  for our region of interest
p.txdb <- autoplot(Mus.musculus, which = ROI) +
  theme_void() +
  scale_x_continuous(expand = expansion(mult = 0, add = 0), limits = c(start(ROI), end(ROI)))

message("Done: Gene annotation")

# CpGs covered ----------------------------------------------------------
# EOI_c_tbl <- EOI_c %>% 
#   as_tibble() %>% 
#   mutate(target = "CpGs intMe")
# 
# p_C <- COI %>% 
#   as_tibble() %>% 
#   mutate(target = "CpGs intMe") %>% 
#   ggplot(aes(start, ES_NO_NWCGW_me_REPMEAN)) +
#   geom_point(color = COLORS_METH[1], size = 1) +
#   geom_point(data = EOI_c_tbl, color = "darkred", size = 2) +
#   scale_y_continuous(limits = c(0,1), expand = expansion(mult = 0, add = 0.0), breaks = c(0,1)) +
#   scale_x_continuous(expand = expansion(mult = 0, add = 0), limits = c(start(ROI), end(ROI))) +
#   labs(x = seqnames(ROI)) +
#   facet_wrap(~target, ncol = 1, strip.position = "right") +
#   theme_EK() +
#   theme(strip.background = element_blank(),
#         axis.text.x = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.position = "none",
#         plot.background = element_blank(), 
#         panel.border = element_blank(),
#         axis.line.y = element_line(),
#         panel.background = element_blank()) 



# Genomic location --------------------------------------------------------
# chr <- as.character(seqnames(ROI))
# COLOR_ROI <- "white"
# names(COLOR_ROI) <- chr
# 
# EOIs <- resize(EOI, width = 1, fix = "center")
# 
# ROI_tbl <- tibble(target = chr,
#          start = c(start(ROI), end(ROI)),
#          score2 = 1)
# 
# p_ROI <- ROI_tbl %>% 
#   ggplot(aes(start, score2)) +
#   geom_area(fill = "white") +
#   scale_y_continuous(limits = c(0, 1), breaks = c(0, 1), expand = expansion(mult = 0, add = 0)) +
#   scale_x_continuous(expand = expansion(mult = 0, add = 0), limits = c(start(ROI), end(ROI)), breaks = c(start(ROI), start(EOIs), end(ROI))) +
#   labs(x = as.character(seqnames(ROI))) +
#   facet_wrap(~target, ncol = 1, strip.position = "right") +
#   theme_EK() +
#   theme(strip.background = element_blank(),
#         axis.title = element_blank(),
#         legend.position = "none",
#         panel.border = element_blank(),
#         axis.line.x = element_line(),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank(),
#         plot.background = element_blank(), 
#         panel.background = element_blank()) 

# CpG antagonist/neutral locations --------------------------------------------------------
chr <- as.character(seqnames(ROI))
COLOR_ROI_CpGs <- COLORS_STATE_v6[c(1,3)]

EOIs <- resize(EOI, width = 1, fix = "center")
EOIs <- EOIs %>% 
  as_tibble() %>% 
  mutate(score2 = 1)
  
ROI_tbl <- tibble(target = chr,
                  start = c(start(ROI), end(ROI)),
                  score2 = 1)
p_ROI_Cs <- EOIs %>% 
  ggplot(aes(xmin = start-10, xmax = start+10, ymin = 0, ymax = 1, fill = state0)) +
  geom_rect(show.legend = F) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 1), expand = expansion(mult = 0, add = 0)) +
  scale_x_continuous(expand = expansion(mult = 0, add = 0), limits = c(start(ROI), end(ROI)), breaks = c(start(ROI), end(ROI))) +
  scale_fill_manual(values = COLOR_ROI_CpGs) +
  # labs(x = as.character(seqnames(ROI))) +
  theme_EK() +
  theme(strip.background = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank()) 

p_ROI_Cs <- p_ROI_Cs +
  geom_vline(xintercept = start(ROI_ori), color = "grey30", linetype = 3, lwd = 1) +
  geom_vline(xintercept = end(ROI_ori), color = "grey30", linetype = 3, lwd = 1) 

message("Done: state CpGs covered")

# Combined plot -----------------------------------------------------------

#Version 1 - with F1 data
if(VERSION == "F1"){
p_final <- p_H1 + p_H2 + p_H3 + p_H4 + 
  p_D + 
  p_5hmC +
  p_W_EOI + 
  p_W_F1_R + p_W_F1_A + 
  p_SMF_F1_R + p_SMF_F1_A + 
  (if(length(p.txdb@ggplot$layers) != 0){ p.txdb@ggplot } else { plot_spacer() }) +
  p_cCRES + 
  (p_ROI_Cs + xlab(ROI$names) + theme(axis.title.x = element_text())) + 
  plot_layout(ncol = 1, heights = c(rep(1,11),1.5,0.5,0.2))

p_final
}

#Version 2 - with SMF TET data
if(VERSION == "h5mC"){
p_final <- p_H1 + p_H2 + p_H3 + p_H4 + 
  p_D + 
  p_5hmC +
  p_W_EOI + 
  p_SMF_WT_CA + p_SMF_WT_5mC + 
  p_SMF_TET_CA + p_SMF_TET_5mC + 
  (if(length(p.txdb@ggplot$layers) != 0){ p.txdb@ggplot } else { plot_spacer() }) +
  p_cCRES + 
  (p_ROI_Cs + xlab(ROI$names) + theme(axis.title.x = element_text())) + 
  plot_layout(ncol = 1, heights = c(rep(1,11),1.5,0.5,0.2))

p_final
}
message("Done: Final plot")

ROI_name <- ROI$names %>% stringr::str_replace(., "/", "-")
ggplot2::ggsave(plot = p_final, filename = paste0(PLOTOUT, paste("IGV_ICR", ROI_name, "all", sep="_"), "_v_", VERSION, "_", DATE, "_v2_width", WIDTH/1000,"kb.pdf"), width = 8, height = 12)
message("Done: Final plot saved")
}
}


