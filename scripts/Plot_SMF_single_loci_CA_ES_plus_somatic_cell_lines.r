#TITLE:   Plotting SMF single locus plots - for CA analysis - ES and somatic cell line
#AUTHOR:  Elisa Kreibich
#DATE:    15-11-2022
#AIM:     Create single locus SMF plots with average and single molecules plots. 
#         For ES and somatic cell lines.

#Set environment -------------
WD <- '/g/krebs/kreibich/Kreibich_2023_5mC_at_enhancers/'

#Load libraries --------------
library(tidyverse)
library(GenomicRanges)
library(plyranges)
library(QuasR)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(SingleMoleculeFootprinting)

source(paste0(WD, 'scripts/utilities/RFunctions_for_SingleMoleculeFootprinting_plotting.r'))
source(paste0(WD, 'scripts/utilities/COLORS.R'))

# Import Arguments --------------------------------------------------------

  
  #Define file path of QuasR input file with all SMF data sets.
  Qinput    <- '/g/krebs/kreibich/HTS/SMF/MM/QuasR/QuasR_rmdup/QuasR_aligned_files_ALL_rmdup.txt'
  
  SAMPLES           <- c("ES_NO") #short ES sample name
  REPOI             <- 1          #replicate of interest to plot
  MySamples_1       <- c("ES_NO_R1", "ES_NO_R2", "ES_NO_R5a6")
  canidateSamples_1 <- c(MySamples_1[REPOI])
  print(canidateSamples_1)
  
  SAMPLENAMES       <- c("C2C12_NO", "MEL_NO", "NP_NO") #short somatic cell sample name
  SAMPLEREPS        <- list(c("R1", "R2"), c("R1", "R2"),  c("R1", "R2"))
  SAMPLES_l         <- lapply(seq_along(SAMPLENAMES), function(x){paste("SMF_MM", SAMPLENAMES[[x]], SAMPLEREPS[[x]], sep = "_")})
  SOI               <- 3          #sample of interest to plot (from SAMPLES_l on top)
  REPOI             <- 1          #replicate of interest to plot
  canidateSamples_3 <- SAMPLES_l[[SOI]][REPOI]
  CL_TREAT          <- str_extract(canidateSamples_3, "\\w{2}(?=_R.{1,}$)")
  print(canidateSamples_3)

  canidateSamples = c(canidateSamples_1, NA, canidateSamples_3)
# Load Data ---------------------------------------------------------------
## BAMs  
  QuasRprj <- qAlign(sampleFile= Qinput, 
                  genome="BSgenome.Mmusculus.UCSC.mm10",
                  paired= "fr",
                  bisulfite="undir")
  QuasRprj@aligner <- "Rbowtie"
  Samples <- QuasR::alignments(QuasRprj)[[1]]$SampleName
  
  ### ES
  Samples_1 <- Samples[grepl(canidateSamples[1], Samples)]
  print(Samples_1)
  ### Cell line
  Samples_3 <- Samples[grepl(canidateSamples[3], Samples)]
  print(Samples_3)
  
## TFBS - load if you want to plot TFBS motifs in the region
  # MotDBf <- readRDS(paste0(WD, 'data/mapped_jaspar2018_ChIP_score10_inBaits_BANP.rds'))
  # names(MotDBf) <- paste('TFBS_', seq_along(MotDBf), sep = '')
  
## Candidate sites 
  #Load in GRange object including states annotation made in scripts/Make_Final_CA_CMH_data_tibbles_cell_line.Rmd
  INPUT_FILE <- '/g/krebs/kreibich/analysis/SMF/revision/Final_data_tibble_CA_2022-10-21_SMF_MM_ES_NO_R1_R2_R5a6_curatedICR_cObin10_cOCMH30_states_gr.rds'
  candidates <- readRDS(INPUT_FILE)
  names(candidates) <- paste0(candidates$bin)
  
  #If a file is loaded in without GRange annotation, run the following code to annotate the genomic information per bin.
  # CpGs in bins
  # bin_gr <- readRDS('paste0(WD, 'data/CGs_over_baits_101bp.rds'))
  # bin_gr$bin <- names(bin_gr)
  # bin_gr <- as_tibble(bin_gr)
  # candidates <- candidates %>% left_join(., bin_gr) %>% as_granges()
  # names(candidates) <- paste0(candidates$bin)
  
# Choose regions of interest (ROIs)  ----------------------------------------
  ROIs <- candidates[grepl("antagonist", candidates$state2)]
  # ROIs <- candidates[grepl("bin_34655$",candidates$bin)]

  #Resize ROI for plotting (401 is the usual plotting width for SMF data)
  ROIs_rsz <- resize(ROIs, 401, fix = "center")

  # Plotting settings ------------------------------------------------------
  #Define which plots you want to create
  WITH_WT        <- TRUE     #Include wildtype ES data
  WITH_CELL_LINE <- TRUE     #Include somatinc cell line data
  WITH_TFBS      <- FALSE    #Sort for TFBS of interest
  PLOT_MOTIFS    <- FALSE    #Plot TF motifs within the genomic window (MotDBf has to be loaded)
  MOTIFS_OI       <- ""       #If PLOT_MOTIFS==TRUE, do you want to plot only a specific TF of interest (e.g. "CTCF")
  
  PLOT_AVG_ONLY  <- TRUE    #Plot only average plots
  PLOT_INDIVIDU  <- TRUE    #Plot individual plots
  PLOT_COMBI     <- TRUE    #Plot WT plus the Cell line
  
  OUTDIR        <- paste0(WD, 'data_results/plots_single_loci/')
  OUTDIR        <- paste0(WD, 'test/')
  dir.create(OUTDIR)
  
  COV           <- 10     #Coverage cutoff for plotting (10x is normally used, as in the single molecule methylation analysis)
  WDW_SIZE      <- 101    #Genomic window of chromatin accessibility (CA) analysis. 
  START         <- "SM_plot" #Start of final filename of all the plots
  ADD           <- ""     #Additional info added to the end of the file name (FILENAME, ADD, ".png")
  WHAT          <- "CA"   #Type of analysis "CA" for chromatin accessibility analysis | "TF" for TFBS analysis
  COLOR_CL_GpC  <- "grey40"
  COLOR_CL_CpG  <- "turquoise3"
  
# SM plotting -------------------------------------------------------------
# Before you start:
#   - define which features you want to include in the plotname (PLOT.OUT)
#   - define which features you want to include in the plot title (TITLE)
#   - define which features you want to include in the plot caption (CAPTION)

i = 3
for(i in seq_along(ROIs)){
  print(i)
  ICR_name  <- ROIs[i]$names %>% str_replace("/", "-")
  if(WHAT == "CA"){
    PLOT.OUT  <- paste0(OUTDIR, paste(START, ROIs[i]$state2, ROIs[i]$bin, sep = "_"))
  } else if (WHAT == "TF"){
    PLOT.OUT  <- paste0(OUTDIR, paste(START, ROIs[i]$state2, ROIs[i]$representative_motif, ROIs[i]$TFBS, sep = "_"))
  }
  
  TITLE     <- paste(ROIs[i]$bin, ROIs[i]$state2, ROIs[i]$chromHMM, sep = " | ")
  CAPTION   <- paste("COR =", round(ROIs[i]$COR, 2), "| -log10(p) = ", round(-log10(ROIs[i]$pval)))
  TFBS      <- ROIs[i]
  
  if(WHAT == "CA"){
    TFBS_center     <- start(TFBS) + (end(TFBS)-start(TFBS))/2
    BinsCoordinates <- IRanges(start = c(TFBS_center-(WDW_SIZE-1)/2),
                               end = c(TFBS_center+(WDW_SIZE-1)/2))
    range_rs        <- resize(ROIs_rsz[i], width = 3, fix = "center")
  }
  #TFBS window
  if(WITH_TFBS == TRUE){
    TFBS_MOTIF.OI   <- str_replace(ROIs_rsz[i]$TFBS_motif, "_", "::")
    TFBSs           <- subsetByOverlaps(MotDBf, ROIs_rsz[i])
    TFBSs           <- subsetByOverlaps(TFBSs, resize(ROIs_rsz[i], width = 30, fix = "center"))
    TFBS            <- TFBSs[grepl(paste0("^", TFBS_MOTIF.OI, "$"), TFBSs$name)]
    if(length(TFBS) != 1){
      TFBS            <- TFBSs[1]
    }
    bins            <- list(c(-35,-25), c(-15,+15), c(+25,+35))
    TFBS_center     <- start(TFBS) + (end(TFBS)-start(TFBS))/2
    BinsCoordinates_TFBS <- IRanges(start = c(TFBS_center+bins[[1]][1], TFBS_center+bins[[2]][1], TFBS_center+bins[[3]][1]),
                                    end = c(TFBS_center+bins[[1]][2], TFBS_center+bins[[2]][2], TFBS_center+bins[[3]][2]))
  } else {
    if(PLOT_MOTIFS == TRUE){
      TFBSs           <- subsetByOverlaps(MotDBf, ROIs_rsz[i])
      TFBSs           <- TFBSs[grepl(paste(MOTIFS_OI, collapse = "|"), TFBSs$name)]
    } else {
      TFBSs            <- NULL
    }
  }
  
## ES ----------
    MethGR_1 <- try(CallContextMethylation(sampleSheet = Qinput, sample = Samples_1, genome = BSgenome.Mmusculus.UCSC.mm10, coverage = COV, RegionOfInterest = ROIs_rsz[i]))
    if("try-error" %in% class(MethGR_1)) {next}
    if(length(subsetByOverlaps(MethGR_1[[1]][[2]], range_rs, ignore.strand = TRUE)) == 0)  {next}
    Sample_name   <- str_remove(Samples_1, "SMF_MM_")
    
    Sort          <- SortReads.CA.EK(MethSM = MethGR_1, BinsCoordinates, range = ROIs_rsz[i], DE = FALSE, WHAT = WHAT) #Simple sorting for CA
    Sort_me       <- SortReads.CA.CGme.EK(MethSM = MethGR_1, BinsCoordinates, range = ROIs_rsz[i], DE = FALSE, WHAT = WHAT) #Double sorting for CA and then for 5mC -- THIS IS THE PREFERED SORTING!
    SM_count      <- length(Sort[[1]]) + length(Sort[[2]])
    
    TFBSs_range   <- GRanges() #empty GRange object for CA anlaysis
    
  ### Plot ES only 
    pAvg_1        <- PlotAvgSMF.EK(MethGR = MethGR_1[[1]], TFBSs = TFBSs_range, range = ROIs_rsz[i]) + 
                        scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.25), name = "SMF", sec.axis = sec_axis(~ ., name = "CG methylation", breaks = seq(0, 1, 0.25))) +
                        labs(title = TITLE, caption = CAPTION) +
                        theme(panel.grid.major = element_blank())
    
    #Additional plot to see CA window
    range_CpG     <- resize(ROIs_rsz[i], width = 3, fix = "center")
    CpG_loc       <- subsetByOverlaps(MethGR_1[[1]][[2]], range_CpG)
    pAvg_1_window <- pAvg_1 + 
      geom_vline(xintercept = c(start(CpG_loc) + ((WDW_SIZE-1)/2), start(CpG_loc) - ((WDW_SIZE-1)/2)), linetype = 3) 


    ### Plot sorting according to CA binding + 5mC
    pSM_CG_1.CA      <- PlotSM.CpG.EK(MethSM = MethGR_1, range = ROIs_rsz[i], SortedReads = Sort_me, accessibility = TRUE, CGme_sort = TRUE, WHAT = WHAT, title = Sample_name)
    
    pSM_GC_1.CA      <- PlotSM.GpC.EK(MethSM = MethGR_1, range = ROIs_rsz[i], SortedReads = Sort_me, accessibility = TRUE, CGme_sort = TRUE, WHAT = WHAT, title = Sample_name)
    
    pStates_1.CA     <- StateQuantificationPlot.GpC.EK(SortedReads = Sort_me, accessibility = TRUE) + labs(subtitle = "CA")
    
    pStatesMeth_1.CA <- StateQuantificationPlot.CpG.EK(MethSM = MethGR_1, SortedReads = Sort_me, range = ROIs_rsz[i], range_width = 101, accessibility = TRUE, WHAT = "TF")
    
    ### Combine and plot - with double sorting (accessibility & 5mC)
    pAvg_1.2        <- pAvg_1 + plot_spacer() + plot_layout(ncol = 2, widths = c(1, 0.75))
    
    pSMs_1.CA       <- (pSM_GC_1.CA + pStates_1.CA + pSM_CG_1.CA + pStatesMeth_1.CA & theme(plot.background = element_blank())) + plot_layout(ncol = 4, widths = c(1,0.2,1,0.2))
    pSMs_1.CA.l     <- (pSM_GC_1.CA + pStates_1.CA + pSM_CG_1.CA + pStatesMeth_1.CA & theme(plot.background = element_blank(), legend.position = "none")) + plot_layout(ncol = 4, widths = c(1,0.2,1,0.2)) #Version without legend when combined with TKO plots

    final_plot_1.CA <- pAvg_1.2 / pSMs_1.CA + plot_layout(nrow = 2, heights = c(1, 0.75))
    # final_plot_1.CA
  
  ### Save plots
    PLOTNAME <- paste0(PLOT.OUT, "_", Samples_1, ADD)
    ggsave(plot = final_plot_1.CA, filename = paste0(PLOTNAME, ".pdf"), height = 8, width = 8)

## CELL LINE -------
    if(WITH_CELL_LINE == TRUE){
      MethGR_3 <- CallContextMethylation(sampleSheet = Qinput, sample = Samples_3, genome = BSgenome.Mmusculus.UCSC.mm10, coverage = COV, RegionOfInterest = ROIs_rsz[i])
      
      if (CL_TREAT == "DE"){DE_INPUT_3 <- TRUE} else if (CL_TREAT == "NO"){DE_INPUT_3 <- FALSE}
      Sample_name3      <- str_remove(Samples_3, "SMF_MM_")
      Sample_name_short <- str_split(Sample_name3, "_", simplify = T)[1]
      
      ### Sorting SM  
      Sort_3        <- SortReads.CA.EK(MethSM = MethGR_3, BinsCoordinates, range = ROIs_rsz[i], DE = DE_INPUT_3, WHAT = WHAT)
      SM_count_3    <- length(Sort_3[[1]]) + length(Sort_3[[2]])
      
      if (CL_TREAT != "DE"){
        Sort_me_3   <- SortReads.CA.CGme.EK(MethSM = MethGR_3, BinsCoordinates, range = ROIs_rsz[i], DE = DE_INPUT_3, WHAT = WHAT)
      }
      
      ### Plot CELL LINE - only CELL LINE
      pAvg_3    <- PlotAvgSMF.EK(MethGR = MethGR_1[[1]], TFBSs = TFBSs_range, range = ROIs_rsz[i],
                                 TETTKO = TRUE, DE = DE_INPUT_3, MethGR_TETTKO = MethGR_3[[1]], TETonly = TRUE, 
                                 CONTEXTS_3 = c(paste(c("GC", "CG"), Sample_name_short)), 
                                 COLORS_3 = c(COLOR_CL_GpC, COLOR_CL_CpG)) +
        scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.25), name = "SMF", sec.axis = sec_axis(~ ., name = "CG methylation", breaks = seq(0, 1, 0.25))) +
        labs(title = TITLE, caption = CAPTION) +
        theme(panel.grid.major = element_blank())
      
      ### Plot CELL LINE - ES + CELL LINE
      pAvg_3.2    <- PlotAvgSMF.EK(MethGR = MethGR_1[[1]], TFBSs = TFBSs_range, range = ROIs_rsz[i],
                                   TETTKO = TRUE, DE = DE_INPUT_3, MethGR_TETTKO = MethGR_3[[1]], 
                                   CONTEXTS_3 = c(paste(c("GC", "CG"), Sample_name_short)), 
                                   COLORS_3 = c(COLOR_CL_GpC, COLOR_CL_CpG)) +
        scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.25), name = "SMF", sec.axis = sec_axis(~ ., name = "CG methylation", breaks = seq(0, 1, 0.25))) +
        labs(title = TITLE, caption = CAPTION) +
        theme(panel.grid.major = element_blank())
      
      ### Plot CELL LINE SM
      if (CL_TREAT != "DE"){
        ### Plot  CELL LINE SM - with 5mC
        pSM_CG_3.CA        <- PlotSM.CpG.EK(MethSM = MethGR_3, range = ROIs_rsz[i], SortedReads = Sort_me_3, accessibility = TRUE, WHAT = WHAT, title = Sample_name3)
        pStatesMeth_3.CA   <- StateQuantificationPlot.CpG.EK(MethSM = MethGR_3, SortedReads = Sort_me_3, range = ROIs_rsz[i], range_width = 201, accessibility = TRUE, WHAT = WHAT)
        
        pSM_GC_3.CA       <- PlotSM.GpC.EK(MethSM = MethGR_3, range = ROIs_rsz[i], SortedReads = Sort_me_3, accessibility = TRUE, WHAT = WHAT, title = Sample_name3)
        pStates_3.CA      <- StateQuantificationPlot.GpC.EK(Sort_me_3, accessibility = TRUE) + labs(subtitle = "CA")
        
      } else if (CL_TREAT == "DE"){
        ### Plot CELL LINE SM - withOUT 5mC
        pSM_GC_3.CA        <- PlotSM.GpC.EK(MethSM = MethGR_3, range = ROIs_rsz[i], SortedReads = Sort_3, DE = DE_INPUT_3, accessibility = TRUE, title = Sample_name3)
        pStates_3.CA       <- StateQuantificationPlot.GpC.EK(Sort_3, accessibility = TRUE) + labs(subtitle = "CA")
      }
      
      ### Plot CELL LINE 5mC Box
      pBox_3            <- Plot.5mC.Box.EK(MethGR1 = MethGR_1[[1]], n_samples = 2, range = ROIs_rsz[i], MethGR2 = MethGR_3[[1]], DE2 = DE_INPUT_3) +
        labs(title = TITLE)
      
      ### Combine and plot - AVERAGE PLOT
      pAvg_3.1          <- pAvg_3 + plot_spacer() + plot_layout(ncol = 2, widths = c(1,0.75))
      pAvg_3.3          <- pAvg_3.2 + plot_spacer() + plot_layout(ncol = 2, widths = c(1,0.75))
      pAvg_box_3        <- pBox_3 + (pAvg_3.2 + theme(title = element_blank())) + plot_layout(ncol = 1, height = c(0.1, 1))
      
      ### Combine and plot
      if (CL_TREAT != "DE"){
        ###  SMF PLOT with double sorting
        pSMs_3.CA       <- (pSM_GC_3.CA + pStates_3.CA + pSM_CG_3.CA + pStatesMeth_3.CA & theme(plot.background = element_blank())) + plot_layout(ncol = 4, widths = c(1,0.2,1,0.2))
        pSMs_3.CA.l     <- (pSM_GC_3.CA + pStates_3.CA + pSM_CG_3.CA + pStatesMeth_3.CA & theme(plot.background = element_blank(), legend.position = "none")) + plot_layout(ncol = 4, widths = c(1,0.2,1,0.2))
        final_plot_3.   <- pAvg_3.1 / pSMs_3.CA + plot_layout(nrow = 2, heights = c(1,0.75))
        
      } else if (CL_TREAT == "DE"){
        ###  SMF PLOT withOUT double sorting & withOUT 5mC
        pSMs_3.CA       <- (pSM_GC_3.CA + pStates_3.CA + plot_spacer() + plot_spacer() & theme(plot.background = element_blank())) + plot_layout(ncol = 4, widths = c(1,0.2,1,0.2))
        pSMs_3.CA.l     <- (pSM_GC_3.CA + pStates_3.CA + plot_spacer() + plot_spacer() & theme(plot.background = element_blank(), legend.position = "none")) + plot_layout(ncol = 4, widths = c(1,0.2,1,0.2))
        final_plot_3.   <- pAvg_3.1 / pSMs_3.CA + plot_layout(nrow = 2, heights = c(1,0.75))
      }
      
      ### Save plots
      PLOTNAME <- paste0(PLOT.OUT, "_", Samples_3, ADD)
      
      if(PLOT_AVG_ONLY == TRUE){
        ggsave(plot = pAvg_box_3, filename = paste0(PLOTNAME, "_avBox.pdf"), height = 5, width = 6)
        # ggsave(plot = pAvg_3, filename = paste0(PLOTNAME, "_av.pdf"), height = 3, width = 4)
        message("DONE: Plot PLOT_AVG_ONLY - ", Samples_3)
      }
      
      if(PLOT_INDIVIDU == TRUE){
        ggsave(plot = final_plot_3., filename = paste0(PLOTNAME, ".pdf"), height = 8, width = 8)
        # ggsave(plot = final_plot_3., filename = paste0(PLOTNAME, ".png"), height = 8, width = 8)
        message("DONE: Plot PLOT_INDIVIDU - ", Samples_3)
      }
      
      if(PLOT_COMBI == TRUE){
        PLOTNAME <- paste0(PLOT.OUT, "_", Samples_1, Samples_3, ADD)
        ### Combo plot SMF 
        final_plot_3.combi  <- pAvg_3.3 / pSMs_1.CA.l / pSMs_3.CA + plot_layout(nrow = 3, heights = c(1, 0.75, 0.75))
        
        ggsave(plot = final_plot_3., filename = paste0(PLOTNAME, "_combi.pdf"), height = 8, width = 8)
        # ggsave(plot = final_plot_3., filename = paste0(PLOTNAME, ".png"), height = 8, width = 8)
        message("DONE: Plot PLOT_COMBI - ", Samples_1, Samples_3)
      }
    }

}
  
  