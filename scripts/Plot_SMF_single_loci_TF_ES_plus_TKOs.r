#TITLE:   Plotting SMF single locus plots - for TF analysis - ES and TKOs
#AUTHOR:  Elisa Kreibich
#DATE:    15-11-2022
#AIM:     Create single locus SMF plots with average and single molecules plots. 
#         For ES and TKOs.
#         TFBS analysis

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
  SAMPLES   <- c("ES_NO", "TKO_DE", "TETTKO_NO")
  REPOI     <- 1
  TKO_TREAT <- "DE"
  TET_TREAT <- "NO"
  
  #Define file path of QuasR input file with all SMF data sets.
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
  
  canidateSamples <- c(MySamples_1[REPOI], MySamples_2[REPOI], MySamples_3[REPOI])
  print(canidateSamples)

  
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
  ### TKO
  Samples_2 <- Samples[grepl(canidateSamples[2], Samples)]
  print(Samples_2)
  ### TETTKO
  Samples_3 <- Samples[grepl(canidateSamples[3], Samples)]
  print(Samples_3)
  
## TFBS
  MotDBf <- readRDS(paste0(WD, 'data/mapped_jaspar2018_ChIP_score10_inBaits_BANP.rds'))
  # MotDBf <- readRDS('/g/krebs/kreibich/DB/mapped_jaspar2018_ChIP_score10_inBaits_BANP.rds')
  names(MotDBf) <- paste('TFBS_', seq_along(MotDBf), sep = '')
  
## Candidate sites 
  #Load in GRange object including states annotation made in scripts/Make_Final_CA_CMH_data_tibbles_cell_line.Rmd
  INPUT_FILE <- '/g/krebs/kreibich/analysis/SMF/TF/rds/Final_data_tibble_TF_CMH_5mC_cov_ES_NO_TKO_DE_TETTKO_NO_BANP_AllFilter_cO10_2022-03-21_gr.rds'
  candidates <- readRDS(INPUT_FILE)
  candidates <- candidates %>% 
    filter(fraction == "bound")
  names(candidates) <- paste0(candidates$TFBS)
  
# Choose regions of interest (ROIs)  ----------------------------------------
  ROIs <- candidates[grepl("antagonist", candidates$state2)]
  # ROIs <- candidates[grepl("bin_34655$",candidates$bin)]

  #Resize ROI for plotting (401 is the usual plotting width for SMF data)
  ROIs_rsz <- resize(ROIs, 401, fix = "center")

  # Plotting settings ------------------------------------------------------
  #Define which plots you want to create
  WITH_WT        <- TRUE     #Include wildtype (WT) data
  WITH_DNMT_TKO  <- TRUE     #Include DNMT TKO data
  WITH_TET_TKO   <- TRUE     #Include TET TKO data
  WITH_TFBS      <- TRUE    #Sort for TFBS of interest. !!!!! TRUE FOR TF ANALYSIS !!!!!
  PLOT_MOTIFS    <- TRUE    #Plot TF motifs within the genomic window (MotDBf has to be loaded)
  MOTIFS_OI      <- ""       #If PLOT_MOTIFS==TRUE, do you want to plot only a specific TF of interest (e.g. "CTCF")
  
  PLOT_AVG_ONLY  <- TRUE    #Plot only average plots
  PLOT_INDIVIDU  <- TRUE    #Plot individual plots
  PLOT_COMBI     <- TRUE    #Plot WT plus the one TKO
  PLOT_COMBI_ALL <- TRUE    #Plot WT plus both TKOs
  
  
  OUTDIR        <- paste0(WD, 'data_results/plots_single_loci/')
  OUTDIR        <- paste0(WD, 'test/')
  dir.create(OUTDIR)
  
  COV           <- 10     #Coverage cutoff for plotting (10x is normally used, as in the single molecule methylation analysis)
  WDW_SIZE      <- 30     #Genomic window of TFBS analysis. 
  START         <- "SM_plot" #Start of final filename of all the plots
  ADD           <- ""     #Additional info added to the end of the file name (FILENAME, ADD, ".png")
  WHAT          <- "TF"   #Type of analysis "CA" for chromatin accessibility analysis | "TF" for TFBS analysis
  ACC           <- FALSE  #Accessibility
  
# SM plotting -------------------------------------------------------------
# Before you start:
#   - define which features you want to include in the plotname (PLOT.OUT)
#   - define which features you want to include in the plot title (TITLE)
#   - define which features you want to include in the plot caption (CAPTION)

i = 2
for(i in seq_along(ROIs)){
  print(i)
  ICR_name  <- ROIs[i]$names %>% str_replace("/", "-")
  if(WHAT == "CA"){
    PLOT.OUT  <- paste0(OUTDIR, paste(START, ROIs[i]$state2, ROIs[i]$bin, sep = "_"))
  } else if (WHAT == "TF"){
    PLOT.OUT  <- paste0(OUTDIR, paste(START, ROIs[i]$state2, ROIs[i]$representative.motif, ROIs[i]$TFBS, sep = "_"))
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
    BinsCoordinates <- IRanges(start = c(TFBS_center+bins[[1]][1], TFBS_center+bins[[2]][1], TFBS_center+bins[[3]][1]),
                                    end = c(TFBS_center+bins[[1]][2], TFBS_center+bins[[2]][2], TFBS_center+bins[[3]][2]))
    
    TFBSs_range           <- subsetByOverlaps(MotDBf, ROIs_rsz[i])
    TFBSs_range           <- TFBSs[grepl(paste(MOTIFS_OI, collapse = "|"), TFBSs$name)]
    
  } else {
    if(PLOT_MOTIFS == TRUE){
      TFBSs_range           <- subsetByOverlaps(MotDBf, ROIs_rsz[i])
      TFBSs_range           <- TFBSs[grepl(paste(MOTIFS_OI, collapse = "|"), TFBSs$name)]
    } else {
      TFBSs_range            <- GRanges()
    }
  }
  
## ES ----------
    MethGR_1 <- try(CallContextMethylation(sampleSheet = Qinput, sample = Samples_1, genome = BSgenome.Mmusculus.UCSC.mm10, coverage = COV, RegionOfInterest = ROIs_rsz[i]))
    if("try-error" %in% class(MethGR_1)) {next}
    if(length(subsetByOverlaps(MethGR_1[[1]][[2]], range_rs, ignore.strand = TRUE)) == 0)  {next}
    Sample_name   <- str_remove(Samples_1, "SMF_MM_")
    
    Sort_me       <- SortReads.TF.CGme.EK(MethSM = MethGR_1, BinsCoordinates, range = ROIs_rsz[i], DE = FALSE, WHAT = WHAT) #Double sorting for CA and then for 5mC -- THIS IS THE PREFERED SORTING!
    SM_count      <- length(unlist(Sort_me))
    
  ### Plot ES only 
    pAvg_1        <- PlotAvgSMF.EK(MethGR = MethGR_1[[1]], TFBSs = TFBSs_range, range = ROIs_rsz[i]) + 
                        scale_y_continuous(limits = c(-0.28, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.25), name = "SMF", sec.axis = sec_axis(~ ., name = "CG methylation", breaks = seq(0, 1, 0.25))) +
                        labs(title = TITLE, caption = CAPTION) +
                        theme(panel.grid.major = element_blank())


    ### Plot sorting according to TF binding + 5mC
    pSM_CG_1.TF      <- PlotSM.CpG.EK(MethSM = MethGR_1, range = ROIs_rsz[i], SortedReads = Sort_me, accessibility = ACC, CGme_sort = TRUE, WHAT = WHAT, title = Sample_name)
    
    pSM_GC_1.TF      <- PlotSM.GpC.EK(MethSM = MethGR_1, range = ROIs_rsz[i], SortedReads = Sort_me, accessibility = ACC, CGme_sort = TRUE, WHAT = WHAT, title = Sample_name)
    
    pStates_1.TF     <- StateQuantificationPlot.GpC.EK(SortedReads = Sort_me, accessibility = ACC) + labs(subtitle = "TF")
    
    pStatesMeth_1.TF <- StateQuantificationPlot.CpG.EK(MethSM = MethGR_1, SortedReads = Sort_me, range = ROIs_rsz[i], range_width = 101, accessibility = ACC, WHAT = "TF")
    
    ### Combine and plot - with double sorting (accessibility & 5mC)
    pAvg_1.2        <- pAvg_1 + plot_spacer() + plot_layout(ncol = 2, widths = c(1, 0.75))
    
    pSMs_1.TF       <- (pSM_GC_1.TF + pStates_1.TF + pSM_CG_1.TF + pStatesMeth_1.TF & theme(plot.background = element_blank())) + plot_layout(ncol = 4, widths = c(1,0.2,1,0.2))
    pSMs_1.TF.l     <- (pSM_GC_1.TF + pStates_1.TF + pSM_CG_1.TF + pStatesMeth_1.TF & theme(plot.background = element_blank(), legend.position = "none")) + plot_layout(ncol = 4, widths = c(1,0.2,1,0.2)) #Version without legend when combined with TKO plots

    final_plot_1.TF <- pAvg_1.2 / pSMs_1.TF + plot_layout(nrow = 2, heights = c(1, 0.75))
    # final_plot_1.TF
  
  ### Save plots
    PLOTNAME <- paste0(PLOT.OUT, "_", Samples_1, ADD)
    ggsave(plot = final_plot_1.TF, filename = paste0(PLOTNAME, ".pdf"), height = 8, width = 8)

## TKO -----------------------------
#For DNMT TKOs we will not plot the 5mC (and also not sort for it), because it is 0 for "NO" and shows accessibility in "DE" datasets.
    if(WITH_DNMT_TKO == TRUE){
      MethGR_2 <- CallContextMethylation(sampleSheet = Qinput, sample = Samples_2, genome = BSgenome.Mmusculus.UCSC.mm10, coverage = COV, RegionOfInterest = ROIs_rsz[i])
      
      if (TKO_TREAT == "DE"){DE_INPUT <- TRUE} else if (TKO_TREAT == "NO"){DE_INPUT <- FALSE}
      Sample_name2   <- str_remove(Samples_2, "SMF_MM_")
      
      ### Sorting SM  
      Sort_2        <- SortReads.TF.CGme.EK(MethSM = MethGR_2, BinsCoordinates, range = ROIs_rsz[i], DE = DE_INPUT, WHAT = WHAT, plusCGme = FALSE)
      SM_count_2    <- length(Sort_2[[1]]) + length(Sort_2[[2]])
      
      ### Plot TKO - only TKO
      pAvg_2    <- PlotAvgSMF.EK(MethGR = MethGR_1[[1]], TFBSs = TFBSs_range, range = ROIs_rsz[i],
                                 TKO = TRUE, DE = DE_INPUT, MethGR_TKO = MethGR_2[[1]], TKOonly = TRUE) +
        scale_y_continuous(limits = c(-0.28, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.25), name = "SMF", sec.axis = sec_axis(~ ., name = "CG methylation", breaks = seq(0, 1, 0.25))) +
        labs(title = TITLE, caption = CAPTION) +
        theme(panel.grid.major = element_blank())
      
      ### Plot TKO - ES + TKO
      pAvg_2.2  <- PlotAvgSMF.EK(MethGR = MethGR_1[[1]], TFBSs = TFBSs_range, range = ROIs_rsz[i],
                                 TKO = TRUE, DE = DE_INPUT, MethGR_TKO = MethGR_2[[1]]) +
        scale_y_continuous(limits = c(-0.28, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.25), name = "SMF", sec.axis = sec_axis(~ ., name = "CG methylation", breaks = seq(0, 1, 0.25))) +
        labs(title = TITLE, caption = CAPTION) +
        theme(panel.grid.major = element_blank())
      
      ### Plot TKO SM
      pSM_GC_2.TF        <- PlotSM.GpC.EK(MethSM = MethGR_2, range = ROIs_rsz[i], SortedReads = Sort_2, DE = DE_INPUT, accessibility = ACC, CGme_sort = FALSE, WHAT = WHAT, title = Sample_name2)
      pStates_2.TF       <- StateQuantificationPlot.GpC.EK(Sort_2, accessibility = ACC) + labs(subtitle = "TF")
      
      ### Plot TKO 5mC Box
      pBox_2            <- Plot.5mC.Box.EK(MethGR1 = MethGR_1[[1]], n_samples = 2, range = ROIs_rsz[i], MethGR2 = MethGR_2[[1]], DE2 = DE_INPUT) +
        labs(title = TITLE)
      
      ### Combine and plot - AVERAGE PLOT
      pAvg_2.1        <- pAvg_2 + plot_spacer() + plot_layout(ncol = 2, widths = c(1,0.75))
      pAvg_2.3        <- pAvg_2.2 + plot_spacer() + plot_layout(ncol = 2, widths = c(1,0.75))
      pAvg_box_2      <- pBox_2 + (pAvg_2.2 + theme(title = element_blank())) + plot_layout(ncol = 1, height = c(0.1, 1))
      
      
      ### Combine and plot
        ###  SMF PLOT withOUT double sorting & withOUT 5mC
        pSMs_2.TF       <- (pSM_GC_2.TF + pStates_2.TF + plot_spacer() + plot_spacer() & theme(plot.background = element_blank())) + plot_layout(ncol = 4, widths = c(1,0.2,1,0.2))
        pSMs_2.TF.l     <- pSMs_2.TF & theme(legend.position = "none")
        final_plot_2.   <- pAvg_2.1 / pSMs_2.TF + plot_layout(nrow = 2, heights = c(1,0.75))
      
      ### Save plots
      PLOTNAME <- paste0(PLOT.OUT, "_", Samples_2, ADD)
      
      if(PLOT_AVG_ONLY == TRUE){
        ggsave(plot = pAvg_box_2, filename = paste0(PLOTNAME, "_avBox.pdf"), height = 5, width = 6)
        # ggsave(plot = pAvg_2, filename = paste0(PLOTNAME, "_av.pdf"), height = 3, width = 4)
        message("DONE: Plot PLOT_AVG_ONLY - ", Samples_2)
      }

      if(PLOT_INDIVIDU == TRUE){
        ggsave(plot = final_plot_2., filename = paste0(PLOTNAME, ".pdf"), height = 8, width = 8)
        # ggsave(plot = final_plot_2., filename = paste0(PLOTNAME, ".png"), height = 8, width = 8)
        message("DONE: Plot PLOT_INDIVIDU - ", Samples_2)
      }
      
      if(PLOT_COMBI == TRUE){
        PLOTNAME <- paste0(PLOT.OUT, "_", Samples_1, Samples_2, ADD)
        ### Combo plot SMF 
        pSMs_2.TF           <- (pSM_GC_2.TF + pStates_2.TF + plot_spacer() + plot_spacer() & theme(plot.background = element_blank())) + plot_layout(ncol = 4, widths = c(1,0.2,1,0.2))
        pSMs_2.TF.l         <- pSMs_2.TF & theme(legend.position = "none")
        final_plot_2.combi  <- pAvg_2.3 / pSMs_1.TF.l / pSMs_2.TF + plot_layout(nrow = 3, heights = c(1, 0.75, 0.75))
        
        ggsave(plot = final_plot_2., filename = paste0(PLOTNAME, "_combi.pdf"), height = 8, width = 8)
        # ggsave(plot = final_plot_2., filename = paste0(PLOTNAME, ".png"), height = 8, width = 8)
        message("DONE: Plot PLOT_COMBI - ", Samples_1, Samples_2)
      }
    }
    
## TET TKO -------
    if(WITH_TET_TKO == TRUE){
      MethGR_3 <- CallContextMethylation(sampleSheet = Qinput, sample = Samples_3, genome = BSgenome.Mmusculus.UCSC.mm10, coverage = COV, RegionOfInterest = ROIs_rsz[i])
      
      if (TET_TREAT == "DE"){DE_INPUT_3 <- TRUE} else if (TET_TREAT == "NO"){DE_INPUT_3 <- FALSE}
      Sample_name3   <- str_remove(Samples_3, "SMF_MM_")
      
      ### Sorting SM  
      Sort_3        <- SortReads.TF.CGme.EK(MethSM = MethGR_3, BinsCoordinates, range = ROIs_rsz[i], DE = DE_INPUT_3, WHAT = WHAT)
      SM_count_3    <- length(Sort_3[[1]]) + length(Sort_3[[2]])
      
      if (TET_TREAT != "DE"){
      Sort_me_3   <- SortReads.TF.CGme.EK(MethSM = MethGR_3, BinsCoordinates, range = ROIs_rsz[i], DE = DE_INPUT_3, WHAT = WHAT, plusCGme = FALSE)
      }
      
      ### Plot TET - only TET
      pAvg_3    <- PlotAvgSMF.EK(MethGR = MethGR_1[[1]], TFBSs = TFBSs_range, range = ROIs_rsz[i],
                                 TETTKO = TRUE, DE = DE_INPUT_3, MethGR_TETTKO = MethGR_3[[1]], TETonly = TRUE) +
        scale_y_continuous(limits = c(-0.28, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.25), name = "SMF", sec.axis = sec_axis(~ ., name = "CG methylation", breaks = seq(0, 1, 0.25))) +
        labs(title = TITLE, caption = CAPTION) +
        theme(panel.grid.major = element_blank())
      
      ### Plot TET - ES + TET
      pAvg_3.2    <- PlotAvgSMF.EK(MethGR = MethGR_1[[1]], TFBSs = TFBSs_range, range = ROIs_rsz[i],
                                   TETTKO = TRUE, DE = DE_INPUT_3, MethGR_TETTKO = MethGR_3[[1]]) +
        scale_y_continuous(limits = c(-0.28, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.25), name = "SMF", sec.axis = sec_axis(~ ., name = "CG methylation", breaks = seq(0, 1, 0.25))) +
        labs(title = TITLE, caption = CAPTION) +
        theme(panel.grid.major = element_blank())
      
      ### Plot TET SM
      if (TET_TREAT != "DE"){
        ### Plot TET SM - with 5mC
        pSM_CG_3.TF        <- PlotSM.CpG.EK(MethSM = MethGR_3, range = ROIs_rsz[i], SortedReads = Sort_me_3, accessibility = ACC, WHAT = WHAT, title = Sample_name3)
        pStatesMeth_3.TF   <- StateQuantificationPlot.CpG.EK(MethSM = MethGR_3, SortedReads = Sort_me_3, range = ROIs_rsz[i], range_width = 201, accessibility = ACC, WHAT = WHAT)
        
        pSM_GC_3.TF       <- PlotSM.GpC.EK(MethSM = MethGR_3, range = ROIs_rsz[i], SortedReads = Sort_me_3, accessibility = ACC, WHAT = WHAT, title = Sample_name3)
        pStates_3.TF      <- StateQuantificationPlot.GpC.EK(Sort_me_3, accessibility = ACC) + labs(subtitle = "TF")
        
      } else if (TET_TREAT == "DE"){
        ### Plot TET SM - withOUT 5mC
        pSM_GC_3.TF        <- PlotSM.GpC.EK(MethSM = MethGR_3, range = ROIs_rsz[i], SortedReads = Sort_3, DE = DE_INPUT_3, accessibility = ACC, title = Sample_name3)
        pStates_3.TF       <- StateQuantificationPlot.GpC.EK(Sort_3, accessibility = ACC) + labs(subtitle = "TF")
      }
      
      ### Plot TET 5mC Box
      pBox_3            <- Plot.5mC.Box.EK(MethGR1 = MethGR_1[[1]], n_samples = 2, range = ROIs_rsz[i], MethGR2 = MethGR_3[[1]], DE2 = DE_INPUT_3) +
        labs(title = TITLE)
      
      ### Combine and plot - AVERAGE PLOT
      pAvg_3.1          <- pAvg_3 + plot_spacer() + plot_layout(ncol = 2, widths = c(1,0.75))
      pAvg_3.3          <- pAvg_3.2 + plot_spacer() + plot_layout(ncol = 2, widths = c(1,0.75))
      pAvg_box_3        <- pBox_3 + (pAvg_3.2 + theme(title = element_blank())) + plot_layout(ncol = 1, height = c(0.1, 1))
      
      ### Combine and plot
      if (TET_TREAT != "DE"){
        ###  SMF PLOT with double sorting
        pSMs_3.TF       <- (pSM_GC_3.TF + pStates_3.TF + pSM_CG_3.TF + pStatesMeth_3.TF & theme(plot.background = element_blank())) + plot_layout(ncol = 4, widths = c(1,0.2,1,0.2))
        pSMs_3.TF.l     <- (pSM_GC_3.TF + pStates_3.TF + pSM_CG_3.TF + pStatesMeth_3.TF & theme(plot.background = element_blank(), legend.position = "none")) + plot_layout(ncol = 4, widths = c(1,0.2,1,0.2))
        final_plot_3.   <- pAvg_3.1 / pSMs_3.TF + plot_layout(nrow = 2, heights = c(1,0.75))
        
      } else if (TET_TREAT == "DE"){
        ###  SMF PLOT withOUT double sorting & withOUT 5mC
        pSMs_3.TF       <- (pSM_GC_3.TF + pStates_3.TF + plot_spacer() + plot_spacer() & theme(plot.background = element_blank())) + plot_layout(ncol = 4, widths = c(1,0.2,1,0.2))
        pSMs_3.TF.l     <- (pSM_GC_3.TF + pStates_3.TF + plot_spacer() + plot_spacer() & theme(plot.background = element_blank(), legend.position = "none")) + plot_layout(ncol = 4, widths = c(1,0.2,1,0.2))
        final_plot_3.   <- pAvg_3.1 / pSMs_3.TF + plot_layout(nrow = 2, heights = c(1,0.75))
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
        final_plot_3.combi  <- pAvg_3.3 / pSMs_1.TF.l / pSMs_3.TF + plot_layout(nrow = 3, heights = c(1, 0.75, 0.75))
        
        ggsave(plot = final_plot_3., filename = paste0(PLOTNAME, "_combi.pdf"), height = 8, width = 8)
        # ggsave(plot = final_plot_3., filename = paste0(PLOTNAME, ".png"), height = 8, width = 8)
        message("DONE: Plot PLOT_COMBI - ", Samples_1, Samples_3)
      }
    }

    
## BOTH TKOs ----------------------------
    if(WITH_DNMT_TKO == TRUE & WITH_TET_TKO == TRUE){
      
      ### Plot ES + TKO + TET
      pAvg_ALL  <- PlotAvgSMF.EK(MethGR = MethGR_1[[1]], TFBSs = TFBSs_range, range = ROIs_rsz[i],
                                 TKO = TRUE, DE = DE_INPUT, MethGR_TKO = MethGR_2[[1]],
                                 TETTKO = TRUE, MethGR_TETTKO = MethGR_3[[1]]) +
        scale_y_continuous(limits = c(-0.28, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.25), name = "SMF", sec.axis = sec_axis(~ ., name = "CG methylation", breaks = seq(0, 1, 0.25))) +
        theme(panel.grid.major = element_blank()) +
        labs(title = TITLE, caption = CAPTION)
      
      pBox_ALL    <- Plot.5mC.Box.EK(MethGR1 = MethGR_1[[1]], n_samples = 3, range = ROIs_rsz[i],
                                     MethGR2 = MethGR_2[[1]], DE2 = DE_INPUT, MethGR3 = MethGR_3[[1]], DE3 = DE_INPUT_3) +
        labs(title = TITLE)
      
      ### Combine and plot ALL
      final_plot_ALL  <- pAvg_ALL + pSMs_1.TF.l + pSMs_2.TF.l + pSMs_3.TF + plot_layout(nrow = 4, heights = c(1.5,1,1,1.3))
      
      pAvg_box_ALL    <- pBox_ALL + (pAvg_ALL + theme(title = element_blank())) + plot_layout(ncol = 1, height = c(0.2, 1))
      
      
      ### Save plots
      PLOTNAME <- paste0(PLOT.OUT, "_", Samples_1, "_", Samples_2, "_", Samples_3, ADD)
      
      if(PLOT_COMBI_ALL == TRUE){
        ggsave(plot = pAvg_box_ALL, filename = paste0(PLOTNAME, "_avBox.pdf"), height = 5, width = 6)
        ggsave(plot = final_plot_ALL, filename = paste0(PLOTNAME, ".pdf"), height = 16, width = 8)
        message("DONE: Plot PLOT_COMBI_ALL - ", Samples_1, Samples_2, Samples_3)
      }
    }
}
  
  