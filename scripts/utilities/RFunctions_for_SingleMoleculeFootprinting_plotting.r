#General Functions -------------------------------------------------
#' Turn sparse single molecule matrix to dense
#'
#' @import magrittr
#' @importFrom testthat equals
#' 
make.SM.mat.dense <- function(binary.matrix){
  
  library(magrittr)
  
  binary.matrix %>%
    as.matrix() %>% 
    replace(equals(.,0), NA) %>% 
    replace(equals(.,1), 0) %>% 
    replace(equals(.,2), 1) -> dense.binary.matrix
  
  return(dense.binary.matrix)
  
}

BinMethylation.EK <- function (MethSM, Bin) {
  binCytosines = colnames(MethSM)[as.numeric(colnames(MethSM)) >= 
                                    start(Bin) & as.numeric(colnames(MethSM)) <= end(Bin)]
  if (length(binCytosines) >= 1) {
    binSummarisedMeth = round((rowMeans(as.matrix(MethSM[, binCytosines]), na.rm = TRUE)))
    binSummarisedMeth = binSummarisedMeth[!(is.na(binSummarisedMeth))]
    return(binSummarisedMeth)
  }
  else if (length(binCytosines) == 0) {
    message(paste0("!!!     [", start(Bin), ";", end(Bin), 
                   "]", " bin overlaps with no covered Cytosines   !!!"))
    return(NA)
  }
}


BinMethylation.OLD <- function (MethSM, Bin) {
  binCytosines = colnames(MethSM)[as.numeric(colnames(MethSM)) >= 
                                    start(Bin) & as.numeric(colnames(MethSM)) <= end(Bin)]
  if (length(binCytosines) > 1) {
    binSummarisedMeth = round(rowMeans(MethSM[, binCytosines],
                                       na.rm = TRUE))
    binSummarisedMeth = binSummarisedMeth[!is.na(binSummarisedMeth)]
  } else if (length(binCytosines) == 1) {
    binSummarisedMeth = MethSM[, binCytosines[1]]
    binSummarisedMeth = binSummarisedMeth[!is.na(binSummarisedMeth)]
  } else if (length(binCytosines) == 0) {
    message(paste0("!!!     [", start(Bin), ";", end(Bin), 
                   "]", " bin overlaps with no covered Cytosines   !!!"))
    return(NA)
  }
}

VectorToMatrix <- function (object, MethGR_object) {
  object <- as.matrix(object)
  colnames(object) <- start(MethGR_object)
  return(object)
  }

VectorizeReads <- function(target,mergedMat){ #vectorize reads from matrices in order to plot in base space
  
  # 			 mergedMat=matList[[i]][naidx,][h$rowInd,]
  GCmet=as.vector(t(mergedMat))
  GCstartP=rep(as.numeric(colnames(mergedMat)),nrow(mergedMat))
  GCreadID=unlist(lapply(seq(nrow(mergedMat)),function(i){rep(i,ncol(mergedMat))}))
  GCgr=GRanges(seqnames(target),IRanges(as.numeric(colnames(mergedMat)),as.numeric(colnames(mergedMat))))
  list(GCstartP,GCreadID,GCmet)
}

OneTFstates   <- function(){
  
  allPos=expand.grid(c(0,1),c(0,1),c(0,1))
  patternStrings=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
  #using only 'pure' states
  states=list(
    bound=patternStrings[6],
    accessible=patternStrings[8],
    closed=patternStrings[1],
    unassigned=patternStrings[!seq_along(patternStrings)%in%c(1,6,8)]
  )
  
  return(states)
  
}

#Plot bulk SMF -------------------------------------------------

PlotAvgSMF.EK <- function (MethGR, TFBSs, range, TKO = NULL, DE = NULL, MethGR_TKO = NULL, TKOonly = NULL, TETTKO = NULL, MethGR_TETTKO = NULL, TETonly = NULL, 
                           CONTEXTS_1 = c("GC WT", "CG WT"),
                           CONTEXTS_2 = c("GC DNMT TKO"),
                           CONTEXTS_3 = c("GC TET TKO", "CG TET TKO"),
                           COLORS_1 = c("black", COLORS_METH[3]),
                           COLORS_2 = c("red3"),
                           COLORS_3 = c("goldenrod", "turquoise3"),
                           SHAPES_1 = c(19, 18),
                           SHAPES_2 = c(19),
                           SHAPES_3 = c(19, 18),
                           SIZE_GC = 2,
                           SIZE_CG = 4) {
  library(tidyverse)
  # CONTEXTS_1 <- c("GC WT", "CG WT")
  # CONTEXTS_2 <- c("GC DNMT TKO")
  # CONTEXTS_3 <- c("GC TET TKO", "CG TET TKO")
  # COLORS_1 <- c("black", COLORS_METH[3])
  # COLORS_2 <- c("red3")
  # COLORS_3 <- c("goldenrod", "turquoise3")
  # SHAPES_1 <- c(19, 18)
  # SHAPES_2 <- c(19)
  # SHAPES_3 <- c(19, 18)
  # SIZE_GC <- 2
  # SIZE_CG <- 4
  
  ## ES
  MethGR_tbl <- lapply(seq(2), function(i){as_tibble(MethGR[[i]])})
  names(MethGR_tbl) <- CONTEXTS_1
  MethGR_tbl <- bind_rows(MethGR_tbl, .id = "context") %>% 
    dplyr::rename(me_mean = ends_with("_MethRate"))
  NAME_1 <- str_remove(names(elementMetadata(MethGR[[1]]))[2], "_MethRate$")
  
  ## TKO
  if(!is.null(TKO)){
    if (DE == TRUE){
    MethGR_TKO      <- MethGR_TKO %>% filter(GenomicContext == "GCH")
    MethGR_TKO_tbl <- as_tibble(MethGR_TKO) %>% 
      mutate(context = case_when(
        GenomicContext == "GCH" ~ CONTEXTS_2[1]
        # GenomicContext == "HCG" ~ CONTEXTS_2[2]
      )) %>% 
      dplyr::rename(me_mean = ends_with("_MethRate"))
    NAME_2 <- str_remove(names(elementMetadata(MethGR_TKO))[2], "_MethRate$")
    
    } else if (DE == FALSE){
      MethGR_TKO_tbl <- lapply(seq(2), function(i){as_tibble(MethGR_TKO[[i]])})
      names(MethGR_TKO_tbl) <- c(CONTEXTS_2, "CG DNMT TKO")
      MethGR_TKO_tbl <- bind_rows(MethGR_TKO_tbl, .id = "context") %>% 
        dplyr::rename(me_mean = ends_with("_MethRate"))
      NAME_2 <- str_remove(names(elementMetadata(MethGR_TKO[[1]]))[2], "_MethRate$")
    }
  }
  ## TET TKO
  if(!is.null(TETTKO)){
    MethGR_TET_tbl <- lapply(seq(2), function(i){as_tibble(MethGR_TETTKO[[i]])})
    names(MethGR_TET_tbl) <- CONTEXTS_3
    MethGR_TET_tbl <- bind_rows(MethGR_TET_tbl, .id = "context") %>% 
      dplyr::rename(me_mean = ends_with("_MethRate"))
    NAME_3 <- str_remove(names(elementMetadata(MethGR_TETTKO[[1]]))[2], "_MethRate$")
  }
  
  plot <- ggplot() +
    ## ES data - GC
    geom_point(data = filter(MethGR_tbl, context == CONTEXTS_1[[1]]), aes(start, 1 - me_mean, color = context, shape = context), size = SIZE_GC) +
    geom_line(data = filter(MethGR_tbl, context == CONTEXTS_1[[1]]), aes(start, 1 - me_mean, color = context)) +
    ## ES data - CG
    geom_point(data = filter(MethGR_tbl, context == CONTEXTS_1[[2]]), aes(start, me_mean, color = context, shape = context), size = SIZE_CG) +
    ## TF motifs
    {if(length(TFBSs) > 0){geom_rect(data = as_tibble(TFBSs), aes(xmin = start, xmax = end, ymin = -0.20 , ymax = -0.26, fill = name), alpha = 0.7, inherit.aes = F)}} +
    ## . ##
    # geom_hline(yintercept = 0, size = 0.2) +
    scale_y_continuous(limits = c(-0.28, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.25), name = "SMF", sec.axis = sec_axis(~ ., name = "CG methylation", breaks = seq(0, 1, 0.25))) +
    scale_x_continuous(limits = c(start(range),end(range)), name = as.character(seqnames(range))) +
    theme_classic() +
    theme(axis.title.y.right = element_text(color = COLORS_METH[3]),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.key.size = unit(0.3, "cm"),
          legend.title = element_blank(),
          legend.margin = margin(-2, 0, 4, 0),
          legend.box.margin = margin(-2, -5, 4, -5),
          axis.text = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA),
          axis.line = element_blank(),
          panel.grid.major = element_line(colour = "grey90"),
          plot.margin = unit(c(10, 5, 0, 5), "pt")) +
    labs(title = paste(range$name2, range$TFBS, range$candidate, sep = " | "))
  
  
  if(is.null(TKO) & is.null(TETTKO) & is.null(TKOonly) & is.null(TETonly)){
    CONTEXTS <- c(CONTEXTS_1)
    COLORS <- c(COLORS_1)
    names(COLORS) <- CONTEXTS
    SHAPES <- c(SHAPES_1)
    names(SHAPES) <- CONTEXTS
    
    plot + 
      ## . ##
      scale_color_manual(values = COLORS) +
      scale_shape_manual(values = SHAPES) +
      labs(subtitle = paste(str_remove(names(elementMetadata(MethGR[[1]]))[2], "_MethRate$")))
    
  } else if (!is.null(TKO) & is.null(TETTKO) & is.null(TKOonly) & is.null(TETonly)) {
    CONTEXTS <- c(CONTEXTS_1, CONTEXTS_2)
    COLORS <- c(COLORS_1, COLORS_2)
    names(COLORS) <- CONTEXTS
    SHAPES <- c(SHAPES_1, SHAPES_2)
    names(SHAPES) <- CONTEXTS
    
    plot + 
      ## TKO data - GC
      geom_point(data = filter(MethGR_TKO_tbl, context == CONTEXTS_2[[1]]), aes(start, 1 - me_mean, color = context, shape = context), size = SIZE_GC) +
      geom_line(data = filter(MethGR_TKO_tbl, context == CONTEXTS_2[[1]]), aes(start, 1 - me_mean, color = context)) +
      ## . ##
      scale_color_manual(values = COLORS) +
      scale_shape_manual(values = SHAPES) +
      labs(subtitle = paste(NAME_1, NAME_2, sep = ' | '))
    
  } else if (is.null(TKO) & !is.null(TETTKO) & is.null(TKOonly) & is.null(TETonly)) {
    CONTEXTS <- c(CONTEXTS_1, CONTEXTS_3)
    COLORS <- c(COLORS_1, COLORS_3)
    names(COLORS) <- CONTEXTS
    SHAPES <- c(SHAPES_1, SHAPES_3)
    names(SHAPES) <- CONTEXTS
    
    plot + 
      # TETTKO data - GC
      geom_point(data = filter(MethGR_TET_tbl, context == CONTEXTS_3[[1]]), aes(start, 1 - me_mean, color = context, shape = context), size = SIZE_GC) +
      geom_line(data = filter(MethGR_TET_tbl, context == CONTEXTS_3[[1]]), aes(start, 1 - me_mean, color = context)) +
      ## TETTKO data - CG
      geom_point(data = filter(MethGR_TET_tbl, context == CONTEXTS_3[[2]]), aes(start, me_mean, color = context, shape = context), size = SIZE_CG) +
      ## . ##
      scale_color_manual(values = COLORS) +
      scale_shape_manual(values = SHAPES) +
      labs(subtitle = paste(NAME_1, NAME_3, sep = ' | '))
    
    
    
  } else if (!is.null(TKO) & !is.null(TETTKO) & is.null(TKOonly) & is.null(TETonly)) {
    CONTEXTS <- c(CONTEXTS_1, CONTEXTS_2, CONTEXTS_3)
    COLORS <- c(COLORS_1, COLORS_2, COLORS_3)
    names(COLORS) <- CONTEXTS
    SHAPES <- c(SHAPES_1, SHAPES_2, SHAPES_3)
    names(SHAPES) <- CONTEXTS
    
    plot + 
      ## TKO data - GC
      geom_point(data = filter(MethGR_TKO_tbl, context == CONTEXTS_2[[1]]), aes(start, 1 - me_mean, color = context, shape = context), size = SIZE_GC) +
      geom_line(data = filter(MethGR_TKO_tbl, context == CONTEXTS_2[[1]]), aes(start, 1 - me_mean, color = context)) +
      # TETTKO data - GC
      geom_point(data = filter(MethGR_TET_tbl, context == CONTEXTS_3[[1]]), aes(start, 1 - me_mean, color = context, shape = context), size = SIZE_GC) +
      geom_line(data = filter(MethGR_TET_tbl, context == CONTEXTS_3[[1]]), aes(start, 1 - me_mean, color = context)) +
      ## TETTKO data - CG
      geom_point(data = filter(MethGR_TET_tbl, context == CONTEXTS_3[[2]]), aes(start, me_mean, color = context, shape = context), size = SIZE_CG) +
      ## . ##
      scale_color_manual(values = COLORS) +
      scale_shape_manual(values = SHAPES) +
      labs(subtitle = paste(NAME_1, NAME_2, NAME_3, sep = ' | '))
    
  } else if (!is.null(TKO) & is.null(TETTKO)  & !is.null(TKOonly) & is.null(TETonly)) {
    CONTEXTS <- c(CONTEXTS_2)
    COLORS <- c(COLORS_2)
    names(COLORS) <- CONTEXTS
    SHAPES <- c(SHAPES_2)
    names(SHAPES) <- CONTEXTS
    
    plot <- ggplot() +
      ## TKO data - GC
      geom_point(data = filter(MethGR_TKO_tbl, context == CONTEXTS_2[[1]]), aes(start, 1 - me_mean, color = context, shape = context), size = SIZE_GC) +
      geom_line(data = filter(MethGR_TKO_tbl, context == CONTEXTS_2[[1]]), aes(start, 1 - me_mean, color = context)) +
      ## TF motifs
      {if(length(TFBSs) > 0){geom_rect(data = as_tibble(TFBSs), aes(xmin = start, xmax = end, ymin = -0.20 , ymax = -0.26, fill = name), alpha = 0.7, inherit.aes = F)}} +
      ## . ##
      # geom_hline(yintercept = 0, size = 0.2) +
      scale_y_continuous(limits = c(-0.28, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.25), name = "SMF", sec.axis = sec_axis(~ ., name = "CG methylation", breaks = seq(0, 1, 0.25))) +
      scale_x_continuous(limits = c(start(range),end(range)), name = as.character(seqnames(range))) +
      scale_color_manual(values = COLORS) +
      scale_shape_manual(values = SHAPES) +
      labs(subtitle = paste(NAME_2, sep = ' | ')) +
      theme_classic() +
      theme(axis.title.y.right = element_text(color = COLORS_METH[3]),
            legend.position = "bottom",
            legend.box = "vertical",
            legend.key.size = unit(0.3, "cm"),
            legend.title = element_blank(),
            legend.margin = margin(-2, 0, 4, 0),
            legend.box.margin = margin(-2, -5, 4, -5),
            axis.text = element_text(color = "black"),
            panel.border = element_rect(color = "black", fill = NA),
            axis.line = element_blank(),
            panel.grid.major = element_line(colour = "grey90"),
            plot.margin = unit(c(10, 5, 0, 5), "pt")) +
      labs(title = paste(range$name2, paste0("TFBS_", range$TFBS), range$candidate, sep = " | "))
    
    plot
    
  } else if (is.null(TKO) & !is.null(TETTKO)  & is.null(TKOonly) & !is.null(TETonly)) {
    CONTEXTS <- c(CONTEXTS_3)
    COLORS <- c(COLORS_3)
    names(COLORS) <- CONTEXTS
    SHAPES <- c(SHAPES_3)
    names(SHAPES) <- CONTEXTS
    
    plot <- ggplot() +
      # TETTKO data - GC
      geom_point(data = filter(MethGR_TET_tbl, context == CONTEXTS_3[[1]]), aes(start, 1 - me_mean, color = context, shape = context), size = SIZE_GC) +
      geom_line(data = filter(MethGR_TET_tbl, context == CONTEXTS_3[[1]]), aes(start, 1 - me_mean, color = context)) +
      ## TETTKO data - CG
      geom_point(data = filter(MethGR_TET_tbl, context == CONTEXTS_3[[2]]), aes(start, me_mean, color = context, shape = context), size = SIZE_CG) +
      ## TF motifs
      {if(length(TFBSs) > 0){geom_rect(data = as_tibble(TFBSs), aes(xmin = start, xmax = end, ymin = -0.20 , ymax = -0.26, fill = name), alpha = 0.7, inherit.aes = F)}} +
      ## . ##
      # geom_hline(yintercept = 0, size = 0.2) +
      scale_y_continuous(limits = c(-0.28, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.25), name = "SMF", sec.axis = sec_axis(~ ., name = "CG methylation", breaks = seq(0, 1, 0.25))) +
      scale_x_continuous(limits = c(start(range),end(range)), name = as.character(seqnames(range))) +
      scale_color_manual(values = COLORS) +
      scale_shape_manual(values = SHAPES) +
      labs(subtitle = paste(NAME_3, sep = ' | ')) +
      theme_classic() +
      theme(axis.title.y.right = element_text(color = COLORS_METH[3]),
            legend.position = "bottom",
            legend.box = "vertical",
            legend.key.size = unit(0.3, "cm"),
            legend.title = element_blank(),
            legend.margin = margin(-2, 0, 4, 0),
            legend.box.margin = margin(-2, -5, 4, -5),
            axis.text = element_text(color = "black"),
            panel.border = element_rect(color = "black", fill = NA),
            axis.line = element_blank(),
            panel.grid.major = element_line(colour = "grey90"),
            plot.margin = unit(c(10, 5, 0, 5), "pt")) 
      # labs(title = paste(range$name2, paste0("TFBS_", range$TFBS), range$candidate, sep = " | "))
    
    plot
  }
}

PlotAvgSMF.REP.EK <- function (MethGR, TFBSs, range, MethGR_REP2 = NULL, REP3 = NULL, MethGR_REP3 = NULL) {
  library(tidyverse)
  CONTEXTS_1 <- c("GC WT REP1", "CG WT REP1")
  CONTEXTS_2 <- c("GC WT REP2", "CG WT REP2")
  CONTEXTS_3 <- c("GC WT REP3", "CG WT REP3")
  
  COLORS_1 <- c("black", COLORS_METH[1])
  COLORS_2 <- c("grey40",COLORS_METH[2])
  COLORS_2 <- c("grey60", COLORS_METH[3])
  
  SHAPES_1 <- c(19, 18)

  if(is.null(REP3)){ 
    CONTEXTS <- c(CONTEXTS_1, CONTEXTS_2)
    COLORS <- c(COLORS_1, COLORS_2)
    SHAPES <- c(SHAPES_1, SHAPES_1)
    names(COLORS) <- CONTEXTS
    names(SHAPES) <- CONTEXTS
  } else {
    CONTEXTS <- c(CONTEXTS_1, CONTEXTS_2, CONTEXTS_3)
    COLORS <- c(COLORS_1, COLORS_2, COLORS_3)
    SHAPES <- c(SHAPES_1, SHAPES_1, SHAPES_1)
    names(COLORS) <- CONTEXTS
    names(SHAPES) <- CONTEXTS
  }
  
  ## REP 1
    MethGR_tbl <- lapply(seq(2), function(i){as_tibble(MethGR[[i]])})
    names(MethGR_tbl) <- CONTEXTS_1
    MethGR_tbl <- bind_rows(MethGR_tbl, .id = "context") %>% 
      dplyr::rename(me_mean = ends_with("_MethRate"))
    
  ## REP 2
    MethGR_REP2_tbl <- lapply(seq(2), function(i){as_tibble(MethGR_REP2[[i]])})
    names(MethGR_REP2_tbl) <- CONTEXTS_2
    MethGR_REP2_tbl <- bind_rows(MethGR_REP2_tbl, .id = "context") %>% 
      dplyr::rename(me_mean = ends_with("_MethRate"))
    
  ## REP 3
    if(!is.null(REP3)){
      MethGR_REP3_tbl <- lapply(seq(2), function(i){as_tibble(MethGR_REP3[[i]])})
      names(MethGR_REP3_tbl) <- CONTEXTS_3
      MethGR_REP3_tbl <- bind_rows(MethGR_REP3_tbl, .id = "context") %>% 
        dplyr::rename(me_mean = ends_with("_MethRate"))
    }
  
  ggplot() +
    ## GC data
    ### REP1
    geom_point(data = filter(MethGR_tbl, context == CONTEXTS_1[[1]]), aes(start, 1 - me_mean, color = context, shape = context), size = 1.3) +
    geom_line(data = filter(MethGR_tbl, context == CONTEXTS_1[[1]]), aes(start, 1 - me_mean, color = context)) +
    ### REP2
    geom_point(data = filter(MethGR_REP2_tbl, context == CONTEXTS_2[[1]]), aes(start, 1 - me_mean, color = context, shape = context), size = 1.3) +
    geom_line(data = filter(MethGR_REP2_tbl, context == CONTEXTS_2[[1]]), aes(start, 1 - me_mean, color = context)) +
    ### REP3
    {if(!is.null(REP3)){geom_point(data = filter(MethGR_REP3_tbl, context == CONTEXTS_3[[1]]), aes(start, 1 - me_mean, color = context, shape = context), size = 1.3)}} +
    {if(!is.null(REP3)){geom_line(data = filter(MethGR_REP3_tbl, context == CONTEXTS_3[[1]]), aes(start, 1 - me_mean, color = context))}} +
    ## CG data
    ### REP1
    geom_point(data = filter(MethGR_tbl, context == CONTEXTS_1[[2]]), aes(start, me_mean, color = context, shape = context), size = 3) +
    ### REP2
    geom_point(data = filter(MethGR_REP2_tbl, context == CONTEXTS_2[[2]]), aes(start, me_mean, color = context, shape = context), size = 3) +
    ### REP3
    {if(!is.null(REP3)){geom_point(data = filter(MethGR_REP3_tbl, context == CONTEXTS_3[[2]]), aes(start, me_mean, color = context, shape = context), size = 3)}} +
    ## TF motifs
    {if(length(TFBSs) > 0){geom_rect(data = as_tibble(TFBSs), aes(xmin = start, xmax = end, ymin = -0.20 , ymax = -0.26, fill = name), alpha=0.7, inherit.aes = F)}} +
    ## . ##
    geom_hline(yintercept = 0, size = 0.2) +
    scale_y_continuous(limits = c(-0.28, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.25), name = "SMF", sec.axis = sec_axis(~ ., name = "CG methylation", breaks = seq(0, 1, 0.25))) +
    scale_x_continuous(limits = c(start(range),end(range)), name = as.character(seqnames(range))) +
    scale_color_manual(values = COLORS) +
    scale_shape_manual(values = SHAPES) +
    theme_classic() +
    theme(axis.title.y.right = element_text(color = COLORS_METH[3]),
          legend.position = "bottom",
          legend.box="vertical",
          legend.key.size = unit(0.3, "cm"),
          legend.title = element_blank(),
          legend.margin = margin(-2, 0, 4, 0),
          legend.box.margin = margin(-2, -5, 4, -5),
          panel.border = element_rect(color = "black", fill = NA),
          axis.line = element_blank(),
          panel.grid.major = element_line(colour = "grey90"),
          plot.margin = unit(c(10, 5, 0, 5), "pt"), 
          plot.subtitle = element_text(size = 7)) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(title = paste(range$name2, paste0("TFBS_", range$TFBS), range$candidate, sep = " | "),
         subtitle = paste(str_remove(names(elementMetadata(MethGR[[1]]))[2], "_MethRate$"),
                          str_remove(names(elementMetadata(MethGR_REP2_tbl[[1]]))[2], "_MethRate$"), 
                          {if(!is.null(REP3)){str_remove(names(elementMetadata(MethGR_REP3_tbl[[1]]))[2], "_MethRate$")}}, sep = ' | '))
  }

PlotAvgSMF.EK_F1 <- function (MethGR, TFBSs, range, 
                              TKO = NULL, DE = NULL, MethGR_TKO = NULL, 
                              TETTKO = NULL, MethGR_TETTKO = NULL, 
                              CONTEXTS_1 = c("GC WT", "CG WT"), CONTEXTS_2 = c("GC DNMT TKO"), CONTEXTS_3 = c("GC TET TKO", "CG TET TKO"),
                              COLORS_1 = c("black", "#01497C"), COLORS_2 = c("#5e0017"), COLORS_3 = c("#c16441", "#468FAF")) {
  # CONTEXTS_1 <- c("GC WT", "CG WT")
  # CONTEXTS_2 <- c("GC DNMT TKO")
  # CONTEXTS_3 <- c("GC TET TKO", "CG TET TKO")
  # COLORS_1 <- c("black", "royalblue3")
  # COLORS_2 <- c("red3")
  # COLORS_3 <- c("goldenrod", "turquoise3")
  SHAPES_1 <- c(19, 18)
  SHAPES_2 <- c(19)
  SHAPES_3 <- c(19, 18)
  SIZE_GC <- 2
  SIZE_CG <- 4
  
  ## ES
  MethGR_tbl <- lapply(seq(2), function(i){as_tibble(MethGR[[i]])})
  names(MethGR_tbl) <- CONTEXTS_1
  MethGR_tbl <- bind_rows(MethGR_tbl, .id = "context") %>% 
    dplyr::rename(me_mean = ends_with("_MethRate"))
  NAME_1 <- str_remove(names(elementMetadata(MethGR[[1]]))[2], "_MethRate$")
  
  ## TKO
  if(!is.null(TKO)){
    if (DE == TRUE){
      MethGR_TKO_tbl <- as_tibble(MethGR_TKO) %>% 
        mutate(context = case_when(
          GenomicContext == "GC" ~ CONTEXTS_2[1],
          GenomicContext == "HCG" ~ CONTEXTS_2[2]
        )) %>% 
        dplyr::rename(me_mean = ends_with("_MethRate"))
      NAME_2 <- str_remove(names(elementMetadata(MethGR_TKO))[2], "_MethRate$")
      
    } else if (DE == FALSE){
      MethGR_TKO_tbl <- lapply(seq(2), function(i){as_tibble(MethGR_TKO[[i]])})
      names(MethGR_TKO_tbl) <- c(CONTEXTS_2, "CG DNMT TKO")
      MethGR_TKO_tbl <- bind_rows(MethGR_TKO_tbl, .id = "context") %>% 
        dplyr::rename(me_mean = ends_with("_MethRate"))
      NAME_2 <- str_remove(names(elementMetadata(MethGR_TKO[[1]]))[2], "_MethRate$")
    }
  }
  ## TET TKO
  if(!is.null(TETTKO)){
    MethGR_TET_tbl <- lapply(seq(2), function(i){as_tibble(MethGR_TETTKO[[i]])})
    names(MethGR_TET_tbl) <- CONTEXTS_3
    MethGR_TET_tbl <- bind_rows(MethGR_TET_tbl, .id = "context") %>% 
      dplyr::rename(me_mean = ends_with("_MethRate"))
    NAME_3 <- str_remove(names(elementMetadata(MethGR_TETTKO[[1]]))[2], "_MethRate$")
  }
  
  plot <- ggplot() +
    ## ES data - GC
    geom_point(data = filter(MethGR_tbl, context == CONTEXTS_1[[1]]), aes(start, 1 - me_mean, color = context, shape = context), size = SIZE_GC) +
    geom_line(data = filter(MethGR_tbl, context == CONTEXTS_1[[1]]), aes(start, 1 - me_mean, color = context)) +
    ## ES data - CG
    geom_point(data = filter(MethGR_tbl, context == CONTEXTS_1[[2]]), aes(start, me_mean, color = context, shape = context), size = SIZE_CG) +
    ## TF motifs
    {if(length(TFBSs) > 0){geom_rect(data = as_tibble(TFBSs), aes(xmin = start, xmax = end, ymin = -0.20 , ymax = -0.26, fill = name), alpha = 0.7, inherit.aes = F)}} +
    ## . ##
    # geom_hline(yintercept = 0, size = 0.2) +
    scale_y_continuous(limits = c(-0.28, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.25), name = "SMF", sec.axis = sec_axis(~ ., name = "CG methylation", breaks = seq(0, 1, 0.25))) +
    scale_x_continuous(limits = c(start(range),end(range)), name = as.character(seqnames(range))) +
    theme_classic() +
    theme(axis.title.y.right = element_text(color = COLORS_METH[3]),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.key.size = unit(0.3, "cm"),
          legend.title = element_blank(),
          legend.margin = margin(-2, 0, 4, 0),
          legend.box.margin = margin(-2, -5, 4, -5),
          panel.border = element_rect(color = "black", fill = NA),
          axis.line = element_blank(),
          panel.grid.major = element_line(colour = "grey90"),
          plot.margin = unit(c(10, 5, 0, 5), "pt")) +
    labs(title = paste(range$name2, paste0("TFBS_", range$TFBS), range$candidate, sep = " | "))
  
  
  if(is.null(TKO) & is.null(TETTKO)){
    CONTEXTS <- c(CONTEXTS_1)
    COLORS <- c(COLORS_1)
    names(COLORS) <- CONTEXTS
    SHAPES <- c(SHAPES_1)
    names(SHAPES) <- CONTEXTS
    
    plot + 
      ## . ##
      scale_color_manual(values = COLORS) +
      scale_shape_manual(values = SHAPES) +
      labs(subtitle = paste(str_remove(names(elementMetadata(MethGR[[1]]))[2], "_MethRate$")))
    
  } else if (!is.null(TKO) & is.null(TETTKO)) {
    CONTEXTS <- c(CONTEXTS_1, CONTEXTS_2)
    COLORS <- c(COLORS_1, COLORS_2)
    names(COLORS) <- CONTEXTS
    SHAPES <- c(SHAPES_1, SHAPES_2)
    names(SHAPES) <- CONTEXTS
    
    plot + 
      ## TKO data - GC
      geom_point(data = filter(MethGR_TKO_tbl, context == CONTEXTS_2[[1]]), aes(start, 1 - me_mean, color = context, shape = context), size = SIZE_GC) +
      geom_line(data = filter(MethGR_TKO_tbl, context == CONTEXTS_2[[1]]), aes(start, 1 - me_mean, color = context)) +
      ## . ##
      scale_color_manual(values = COLORS) +
      scale_shape_manual(values = SHAPES) +
      labs(subtitle = paste(NAME_1, NAME_2, sep = ' | '))
    
  } else if (is.null(TKO) & !is.null(TETTKO)) {
    CONTEXTS <- c(CONTEXTS_1, CONTEXTS_3)
    COLORS <- c(COLORS_1, COLORS_3)
    names(COLORS) <- CONTEXTS
    SHAPES <- c(SHAPES_1, SHAPES_3)
    names(SHAPES) <- CONTEXTS
    
    plot + 
      # TETTKO data - GC
      geom_point(data = filter(MethGR_TET_tbl, context == CONTEXTS_3[[1]]), aes(start, 1 - me_mean, color = context, shape = context), size = SIZE_GC) +
      geom_line(data = filter(MethGR_TET_tbl, context == CONTEXTS_3[[1]]), aes(start, 1 - me_mean, color = context)) +
      ## TETTKO data - CG
      geom_point(data = filter(MethGR_TET_tbl, context == CONTEXTS_3[[2]]), aes(start, me_mean, color = context, shape = context), size = SIZE_CG) +
      ## . ##
      scale_color_manual(values = COLORS) +
      scale_shape_manual(values = SHAPES) +
      labs(subtitle = paste(NAME_1, NAME_3, sep = ' | '))
    
    
    
  } else if (!is.null(TKO) & !is.null(TETTKO)) {
    CONTEXTS <- c(CONTEXTS_1, CONTEXTS_2, CONTEXTS_3)
    COLORS <- c(COLORS_1, COLORS_2, COLORS_3)
    names(COLORS) <- CONTEXTS
    SHAPES <- c(SHAPES_1, SHAPES_2, SHAPES_3)
    names(SHAPES) <- CONTEXTS
    
    plot + 
      ## TKO data - GC
      geom_point(data = filter(MethGR_TKO_tbl, context == CONTEXTS_2[[1]]), aes(start, 1 - me_mean, color = context, shape = context), size = SIZE_GC) +
      geom_line(data = filter(MethGR_TKO_tbl, context == CONTEXTS_2[[1]]), aes(start, 1 - me_mean, color = context)) +
      # TETTKO data - GC
      geom_point(data = filter(MethGR_TET_tbl, context == CONTEXTS_3[[1]]), aes(start, 1 - me_mean, color = context, shape = context), size = SIZE_GC) +
      geom_line(data = filter(MethGR_TET_tbl, context == CONTEXTS_3[[1]]), aes(start, 1 - me_mean, color = context)) +
      ## TETTKO data - CG
      geom_point(data = filter(MethGR_TET_tbl, context == CONTEXTS_3[[2]]), aes(start, me_mean, color = context, shape = context), size = SIZE_CG) +
      ## . ##
      scale_color_manual(values = COLORS) +
      scale_shape_manual(values = SHAPES) +
      labs(subtitle = paste(NAME_1, NAME_2, NAME_3, sep = ' | '))
    
  }
}


#Plot SM stacks (DONE) -------------------------------------------------

PlotSingleMoleculeStack.EK <- function (MethSM, range, context = NULL, title = NULL) {
  if (is.null(context)){
    message("ERROR: Cytosine context has to be determined ('GC' | 'CG').")
  } else if (context == "CG"){
    STATES        <- c("unmethylated", "methylated")
    COLORS        <- c(COLORS_METH[8], COLORS_METH[2])
    names(COLORS) <- STATES
    WIDTH_TILE    <- 8 
    FRACTION_NAME <- "CpGs"
  } else if (context == "GC"){
    STATES        <- c("accessible", "inaccessible")
    COLORS        <- c("grey", "black")
    names(COLORS) <- STATES
    WIDTH_TILE    <- 5
    FRACTION_NAME <- "GpCs"
  } 
  vR1           <- VectorizeReads(range, MethSM)
  names(vR1)    <- c("coordinate", "reads", "state")
  vR_mat1       <- bind_cols(vR1) 
  if(context == "CG"){
    vR_mat      <- vR_mat1 %>% 
      mutate(state = factor(state, levels = c(1,0), labels = rev(STATES)))
  } else if (context == "GC"){
    vR_mat      <- vR_mat1 %>% 
      mutate(state = factor(state, levels = c(0,1), labels = rev(STATES)))
  }
  count_SM      <- rev(vR1[[2]])[1]
  
  vR_mat %>%
    ggplot(aes(coordinate, reads)) +
    geom_tile(aes(fill = as.factor(state)), height = 1, width = WIDTH_TILE, show.legend = T) +
    scale_fill_manual(values = COLORS, labels = STATES, na.value = "white") +
    scale_x_continuous(limits = c(start(range), end(range)), breaks = seq(round(start(range),-2), round(end(range),-2), length.out = 3), name = "") +
    scale_y_continuous(expand = expansion(mult = 0, add = 0), name = paste(count_SM, "molecules")) +
    labs(fill = "",
         subtitle = paste(title, "|", FRACTION_NAME)) +
    guides(fill = guide_legend(ncol = 1)) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          legend.title = element_blank(),
          legend.key.size = unit(0.3, "cm"),
          text = element_text(color = "black"),
          axis.title.x = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          # plot.margin = unit(c(0,1,0,0.2), "cm"),
          axis.line = element_blank())
}

PlotSM.CpG.EK <- function (MethSM, range, SortedReads, DE = FALSE, accessibility = FALSE, CGme_sort = FALSE, WHAT = NULL, title = NULL) {
  #MethSM = output of CallContextMethylation
  #range = region of interest to plot (e.g. 401 bp around CpG)
  #SortedReads = list of sorted reads into different fractions (0,1); output of SortReads.CA.EK (only CA sorting) or SortReads.CA.CGme.EK (CA + 5mC double sorting) 
  #DE = SMF was performed with double enzyme (DE)? TRUE/FALSE
  #accessibility = sorting was performed by accessibility (e.g. by SortReads.CA.EK or SortReads.CA.CGme.EK)? TRUE/FALSE 
  #CGme_sort = sorting was performed also based on 5mC (SortReads.CA.CGme.EK)? TRUE/FALSE
  #WHAT = Type of analysis. "CA"/ "TF"
  #title = name of sample on top of the plot
  
  library(SingleMoleculeFootprinting)
  if (accessibility == TRUE){
    MethGR.    <- MethSM[[1]][[2]]
    MethSM_2  <- Extract.MethSM(MethSM, range, DE = FALSE, WHAT)[[2]] 
  } 
  else if (accessibility == FALSE){
    #CHECK THIS
    MethGR.   <- MethSM[[1]][[2]]
    MethSM_2  <- MethSM[[2]][[1]][[2]] %>% make.SM.mat.dense()
  }
  if (!is.matrix(MethSM_2)) {
    MethSM_2 <- as.matrix(MethSM_2)
    colnames(MethSM_2) <- start(MethGR.)
  }
  if (is.null(SortedReads)) {
    message("No sorting passed or specified, will plot unsorted reads")
    message("Plotting SM stack for CpGs")
    PlotSingleMoleculeStack.EK(MethSM_2, range, context = "CG", title)
  } 
  else if (is.list(SortedReads)) {
    message("Sorting reads according to passed values before plotting")
    if (accessibility == TRUE){
      message("Inferring sorting was performed by accessibility")
      OrderedReads = SortedReads
      names(OrderedReads) = NULL
    } 
    else if (accessibility == FALSE){
      if (CGme_sort == TRUE){
        message("Inferring sorting was performed by single TF plus CGme sort")
        OrderedReads = SortedReads
        OrderedReads = rev(OrderedReads)
      } 
      else if (CGme_sort == FALSE){
        if (length(SortedReads) <= 8) {
          message("Inferring sorting was performed by single TF")
          OrderedReads = SortedReads[as.character(unlist(OneTFstates()))]
          OrderedReads = rev(OrderedReads)
        } 
        else if (length(SortedReads) > 8 & length(SortedReads) <= 16) {
          message("Inferring sorting was performed by TF pair")
          OrderedReads = SortedReads[as.character(unlist(TFpairStates()))]
        }
      }
    }
    
    message("Plotting SM stack for CpGs")
    MethSM_3 = MethSM_2[match(unlist(OrderedReads), rownames(MethSM_2)), ]
    if (!is.matrix(MethSM_3)) {
      MethSM_3 <- as.matrix(MethSM_3)
      colnames(MethSM_3) <- start(MethGR.)
    }    
    PlotSingleMoleculeStack.EK(MethSM_3, range, context = "CG", title)
  }
}

PlotSM.GpC.EK <- function (MethSM, range, SortedReads, DE = FALSE, accessibility = FALSE, CGme_sort = FALSE, WHAT = NULL, title = NULL) {
 #MethSM = output of CallContextMethylation
 #range = region of interest to plot (e.g. 401 bp around CpG)
 #SortedReads = list of sorted reads into different fractions (0,1); output of SortReads.CA.EK (only CA sorting) or SortReads.CA.CGme.EK (CA + 5mC double sorting) 
 #DE = SMF was performed with double enzyme (DE)? TRUE/FALSE
 #accessibility = sorting was performed by accessibility (e.g. by SortReads.CA.EK or SortReads.CA.CGme.EK)? TRUE/FALSE 
 #CGme_sort = sorting was performed also based on 5mC (SortReads.CA.CGme.EK)? TRUE/FALSE
 #WHAT = Type of analysis. "CA"/ "TF"
 #title = name of sample on top of the plot
 
   library(SingleMoleculeFootprinting)
  MethSM_2 <- Extract.MethSM(MethSM, range, DE = DE, WHAT)[[1]] 
  
  if (is.null(SortedReads)) {
    message("No sorting passed or specified, will plot unsorted reads")
    
    message("Plotting SM stack for GpCs")
    PlotSingleMoleculeStack.EK(MethSM_2, range, context = "GC", title)
  } else if (is.list(SortedReads)) {
    message("Sorting reads according to passed values before plotting")
    if (accessibility == TRUE){
      message("Inferring sorting was performed by accessibility")
      OrderedReads = SortedReads
      names(OrderedReads) = NULL
      # OrderedReads = rev(OrderedReads)
    } else if (accessibility == FALSE){ 
      if (CGme_sort == TRUE){
        message("Inferring sorting was performed by single TF and CGme")
        OrderedReads = SortedReads
        OrderedReads = rev(OrderedReads)
      } else if (CGme_sort == FALSE){
        if (length(SortedReads) <= 8) {
          message("Inferring sorting was performed by single TF")
          OrderedReads = SortedReads[as.character(unlist(OneTFstates()))]
          OrderedReads = rev(OrderedReads)
        } else if (length(SortedReads) > 8 & length(SortedReads) <= 16) {
          message("Inferring sorting was performed by TF pair")
          OrderedReads = SortedReads[as.character(unlist(TFpairStates()))]
        }
      }
    }
    
    message("Plotting SM stack for GpCs")
    MethSM_3 = MethSM_2[match(unlist(OrderedReads), rownames(MethSM_2)), ]
    PlotSingleMoleculeStack.EK(MethSM = MethSM_3, range, context = "GC", title)
  } else if (SortedReads == "HC") {
  message("Perfoming hierarchical clustering on single molecules before plotting")
  MethSM_3 = HierarchicalClustering(MethSM_2)   
  
  message("Plotting SM stack for GpCs")
  PlotSingleMoleculeStack.EK(MethSM_3, range, context = "GC", title)
  }
}

#Plot quantification -------------------------------------------------
PlotStateQuantification.GpC.EK <- function (states, what = NULL, OrderedReads) {
  if (is.null(what)){
    message("ERROR: Plotting context has to be determined ('TF' (single transcription factor) | 'CA' (chromatin accessibility)).")
  }
  else if (what == "CA"){
    GroupedCounts         <- unlist(lapply(seq_along(states), function(i) {length(unlist(OrderedReads[as.character(states[[i]])]))}))
    names(GroupedCounts)  <- names(states)
    Colors                <- c("grey", "black")
    names(Colors)         <- names(states)
    Colors_text           <- c("black", "white")
    names(Colors_text)    <- names(states)
  }
  else if (what == "TF") {
    GroupedCounts         <- unlist(lapply(seq_along(states), function(i) {length(unlist(OrderedReads[states[[i]]]))}))
    names(GroupedCounts)  <- names(states)
    if(exists("COLOR_SMF")){
      Colors              <- COLOR_SMF
    } else {
      Colors              <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(9)[c(4, 3, 2)]
    }
    names(Colors)         <- names(states)
    Colors_text           <- rep(c("white"), length(states))
    names(Colors_text)    <- names(states)
  }
  
  GroupedCounts           <- rev(GroupedCounts)
  states_reads            <- lapply(seq(names(GroupedCounts)), function(x){rep(names(GroupedCounts)[x], GroupedCounts[x])}) %>% unlist()
  
  GroupedCounts_df        <- tibble("reads" = seq(sum(GroupedCounts)), 
                                    "x_coord" = rep(1, sum(GroupedCounts)),
                                    "states" = states_reads) %>% 
    add_count(states, name = "sum") %>% 
    add_count(name = "total") %>% 
    mutate(ratio = sum/total*100) %>% 
    group_by(states) %>% 
    mutate(y_coord = mean(reads))
  
  #Plotting
  GroupedCounts_df %>% 
    ggplot() +
    geom_tile(aes(x_coord, reads, fill = states), height = 1, width = 1, show.legend = T) +
    geom_text(data = distinct(GroupedCounts_df, x_coord, y_coord, states, ratio), 
              aes(x_coord, y_coord, label = paste(round(ratio), "%"), color = states), size = 3, show.legend = FALSE) +
    scale_fill_manual(values = Colors) +
    scale_color_manual(values = Colors_text) +
    scale_x_continuous(expand = expansion(mult = 0, add = 0)) +     
    scale_y_continuous(expand = expansion(mult = 0, add = 0)) +
    labs(subtitle = "states",
         x = "") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 1),
           color = "none") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.box="vertical",
      legend.key.size = unit(0.3, "cm"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      # plot.margin = unit(c(5.5, 0, 5.5, 0), "pt"),
      axis.line = element_blank()
      )
}
StateQuantificationPlot.GpC.EK <- function (SortedReads, accessibility = FALSE) {
  if (accessibility == FALSE) {
    if (length(SortedReads) <= 8) {
      message("Inferring sorting was performed by single TF")
      states = OneTFstates()
      states$closed = c(states$closed, states$unassigned)
      states = states[1:3]
      OrderedReads = SortedReads[as.character(unlist(states))]
      PlotStateQuantification.GpC.EK(states, what = "TF", OrderedReads)
    } 
    
    else if (length(SortedReads) > 8 & length(SortedReads) <= 16) {
      message("Inferring sorting was performed by TF pair")
      states = TFpairStates()
      OrderedReads = SortedReads[as.character(unlist(states))]
      PlotStateQuantification.GpC.EK(states, what = "TF", OrderedReads)
    }
  } 
  else if (accessibility == TRUE) {
    message("Inferring sorting was performed by accessibility")
    states = list("accessible" = c(1), "inaccessible" = c(0))
    OrderedReads = SortedReads[as.character(unlist(states))]
    PlotStateQuantification.GpC.EK(states, what = "CA", OrderedReads)
  }
}

PlotStateQuantification.CpG.EK <- function (states = NULL, OrderedReads = NULL, MethSMst, MethGR){
  CG_pos        <- start(MethGR)
  CG_pos        <- as_tibble(cbind("CG" = CG_pos, "x_coord" = seq(length(CG_pos))))
  if(is.null(states)){
    MethSMst_tbl  <- as_tibble(MethSMst) %>% 
      pivot_longer(cols = everything(), names_to = "CG", values_to = "CGme_mean") %>% 
      mutate(CG = as.numeric(CG)) %>% 
      group_by(CG) %>% 
      mutate(CGme_mean = mean(CGme_mean, na.rm = TRUE)) %>% 
      distinct() %>%
      left_join(., CG_pos)
    
    GroupedCounts <- MethSMst %>% 
      as_tibble() %>%
      mutate(reads = seq(nrow(MethSMst))) %>% 
      pivot_longer(cols = -c(reads), names_to = "CG", values_to = "CGme") %>%
      mutate(CG = as.numeric(CG)) %>%
      left_join(., CG_pos)
    
    GroupedCounts_df <- GroupedCounts %>% 
      mutate(y_coord = mean(reads)) %>% 
      left_join(., MethSMst_tbl)
    
    # Colors <- brewer.pal(9, 'Blues')
    Colors <- COLORS_METH

    
  } 
  else {
    states_tbl    <- as_tibble(states) %>% pivot_longer(cols = everything(), names_to = "states", values_to = "combi") %>% distinct() %>% mutate(combi = as.character(combi))
    MethSMst_tbl  <- bind_rows(MethSMst, .id = "combi") %>% 
      left_join(., states_tbl) %>% 
      select(-combi) %>% 
      pivot_longer(cols = -states, names_to = "CG", values_to = "CGme_mean") %>% 
      mutate(CG = as.numeric(CG)) %>% 
      group_by(states, CG) %>% 
      mutate(CGme_mean = mean(CGme_mean, na.rm = TRUE)) %>% 
      distinct() %>%
      left_join(., CG_pos)
    
    GroupedCounts <- unlist(lapply(seq_along(states), function(i) {length(unlist(OrderedReads[as.character(states[[i]])]))}))
    names(GroupedCounts) <- names(states)
    Colors <- brewer.pal(9, 'Blues')
    GroupedCounts <- rev(GroupedCounts)
    
    states_reads <- lapply(seq(names(GroupedCounts)), function(x){rep(names(GroupedCounts)[x], GroupedCounts[x])}) %>% unlist()
    
    GroupedCounts_df <- tibble("reads" = seq(sum(GroupedCounts)), 
                               "states" = states_reads) %>% 
      group_by(states) %>% 
      mutate(y_coord = mean(reads)) %>% 
      left_join(., MethSMst_tbl)
  }
    #Plotting
    GroupedCounts_df %>% 
      ggplot() +
      geom_tile(aes(x_coord, reads, fill = CGme_mean), height = 1, width = 1, show.legend = T) +
      (if(is.null(states)){geom_text(data = distinct(GroupedCounts_df, x_coord, y_coord, states, CGme_mean), aes(x_coord, y_coord, label = round(CGme_mean,2)), size = 3, show.legend = F)} 
       else {geom_text(data = distinct(GroupedCounts_df, x_coord, y_coord, CGme_mean), aes(x_coord, y_coord, label = round(CGme_mean,2)), size = 3, show.legend = F)}) +
      scale_fill_gradientn(colours = Colors, na.value = "transparent",
                           breaks = c(0, 0.5, 1),
                           limits = c(0,1)) +
      scale_x_continuous(expand = expansion(mult = 0, add = 0)) +     
      scale_y_continuous(expand = expansion(mult = 0, add = 0)) +
      labs(subtitle = "5mC") +
      guides(fill = guide_colorbar(title.theme = element_blank())) +
      theme_classic() +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        # plot.margin = unit(c(5.5, 0, 5.5, 0), "pt"),
        axis.line = element_blank())
  }

StateQuantificationPlot.CpG.EK <- function(MethSM, SortedReads = NULL, range, range_width, accessibility = FALSE, WHAT = NULL) {
  library(SingleMoleculeFootprinting)
  MethGR.     <- MethSM[[1]][[2]]
  MethSM_2    <- Extract.MethSM(MethSM, range, DE = FALSE, WHAT)[[2]] 
  
  if (!is.matrix(MethSM_2)) {
    MethSM_2 <- as.matrix(MethSM_2)
    colnames(MethSM_2) <- start(MethGR.)
  }
  
  #Focus on CGs within range_width around central CpG
  range_wd  <- resize(range, range_width, fix = "center")
  MethGR2   <- subsetByOverlaps(MethGR., range_wd, ignore.strand = T)
  MethSM_3  <- MethSM_2[, colnames(MethSM_2) %in% start(MethGR2)]
  if (!is.matrix(MethSM_3)) {
    MethSM_3 <- as.matrix(MethSM_3)
    colnames(MethSM_3) <- start(MethGR2)
  }
  
  #Plotting according to read sorting
  if (accessibility == FALSE){
    if (is.null(SortedReads)) {
      message("No sorting passed or specified, will plot unsorted reads")
      MethSMst <- MethSM1
      PlotStateQuantification.CpG.EK(MethSMst = MethSMst, MethGR = MethGR2)
      
    } else if (length(SortedReads) <= 8) {
      message("Inferring sorting was performed by single TF")
        states = OneTFstates()
        states$closed = c(states$closed, states$unassigned)
        states = states[1:3]
        if(is.list(SortedReads[[1]]) | is.list(SortedReads[[2]])){
          SortedReads <- lapply(SortedReads, unlist)
        }
        OrderedReads = SortedReads[as.character(unlist(states))]
        
        MethSMst <- lapply(seq(OrderedReads), function(x){
          MethSMst <- MethSM_3[rownames(MethSM_3) %in% OrderedReads[[x]],]
          if (!is.matrix(MethSMst)) { # Fix for when only 1 CpG
            if(sum(grepl("^\\d", names(MethSMst))) > 0){ # Fix for when only 1 read in category
              MethSMst <- t(as.matrix(MethSMst))
              rownames(MethSMst) <- OrderedReads[[x]]
            } else {
            MethSMst <- as.matrix(MethSMst)
            colnames(MethSMst) <- start(MethGR2)
            }}
          MethSMst <- colMeans(MethSMst, na.rm = TRUE)
        })
        names(MethSMst) <- as.character(unlist(states))
        PlotStateQuantification.CpG.EK(states, OrderedReads, MethSMst, MethGR = MethGR2)
        
    } else if (length(SortedReads) > 8 & length(SortedReads) <= 16) {
      message("Inferring sorting was performed by TF pair")
      message("ERROR: Plotting of CpG state for TF pair not implemented!")
      return(NA)
      # states = TFpairStates()
      # OrderedReads = SortedReads[as.character(unlist(states))]
      # # PlotStateQuantification.CpG.EK(states, what = "TF", OrderedReads)
      # PlotStateQuantification.CpG.EK(states, OrderedReads, MethSMst, MethGR = MethGR2)
    }
  } else if (accessibility == TRUE){
    message("Inferring sorting was performed by accessibility")
    states = list("accessible" = c(1), "inaccessible" = c(0))
    if(is.list(SortedReads[[1]]) | is.list(SortedReads[[2]])){
      SortedReads <- lapply(SortedReads, unlist)
    }
    OrderedReads = SortedReads[as.character(unlist(states))]
    
    MethSMst <- lapply(seq(OrderedReads), function(x){
      MethSMst <- MethSM_3[row.names(MethSM_3) %in% OrderedReads[[x]],]
      if (!is.matrix(MethSMst)) {
        MethSMst <- as.matrix(MethSMst)
        colnames(MethSMst) <- start(MethGR2)
      }
      MethSMst <- colMeans(MethSMst, na.rm = TRUE)
    })
    names(MethSMst) <- as.character(unlist(states))
    message("Plot with PlotStateQuantification.CpG.EK")
    PlotStateQuantification.CpG.EK(states, OrderedReads, MethSMst, MethGR = MethGR2)
  } 
}


#Sort reads -------------------------------------------------
Extract.MethSM <- function(MethSM, range, DE = NULL, WHAT = NULL) {
  if(is.null(DE)){message("Define type of SMF (DE) as 'TRUE' or 'FALSE'.")}
  if(is.null(WHAT)){message("Define type of analysis (WHAT) as 'CA' or 'TF'.")}
  if(WHAT == "CA") {WIDTH = 3} else {WIDTH=32}
  if (DE == FALSE){
    MethGR      <- MethSM[[1]][[2]]
    MethSM_CG   <- MethSM[[2]][[1]][[2]] %>% make.SM.mat.dense()
    MethSM_GC   <- MethSM[[2]][[1]][[1]] %>% make.SM.mat.dense()
    range_rs    <- resize(range, width = WIDTH, fix = "center")
    CG_OI       <- start(subsetByOverlaps(MethGR, range_rs, ignore.strand = TRUE))
    
    message("Filter for reads that have information for CG of interest")
    MethSM_CG_nonNA <- MethSM_CG[,grepl(CG_OI, colnames(MethSM_CG))]
    MethSM_CG_nonNA_reads <- names(MethSM_CG_nonNA[!is.na(MethSM_CG_nonNA)])
    MethSM_GC_2 <- MethSM_GC[rownames(MethSM_GC) %in% MethSM_CG_nonNA_reads, ]
    MethSM_CG_2 <- MethSM_CG[rownames(MethSM_CG) %in% MethSM_CG_nonNA_reads, ]
    if(!is.matrix(MethSM_CG_2)){
      MethSM_CG_2 <- as.matrix(MethSM_CG_2)
      colnames(MethSM_CG_2) <- CG_OI
    }
    
    MethSM_2 <- list(MethSM_GC_2, MethSM_CG_2, CG_OI)
    names(MethSM_2) <- c("GC", "CG", "CG_OI")
  } else if (DE == TRUE){
    #Define location of CG of interest
    range_rs    <- resize(range, width = WIDTH, fix = "center")
    CG_OI       <- start(subsetByOverlaps(plyranges::filter(MethSM[[1]], GenomicContext == "HCG"), range_rs, ignore.strand = TRUE))
    
    if(length(CG_OI) > 1){
      CG_OI_test <- lapply(seq_along(CG_OI), function(i){
        seq   <- Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm10, seqnames(range_rs), CG_OI[i]-2, CG_OI[i]+2)
        count <- countPattern(DNAString("NWCGW"), seq, fixed = "subject")      
        if(count != 0){
          return(TRUE)
        } else {
          return(FALSE)
        }
      })
      CG_OI <-  CG_OI[unlist(CG_OI_test)]
    }

    
    # Keep only DGCHN | remove CGs
    MethSM_GC   <- MethSM[[2]][[1]] %>% make.SM.mat.dense()
    GCs         <- start(filter(MethSM[[1]], GenomicContext == "GCH"))
    MethSM_GC_2 <- MethSM_GC[, colnames(MethSM_GC) %in% GCs]

    message("Filter for reads that have information for CG of interest")
    MethSM_CG_nonNA <- MethSM_GC[,grepl(CG_OI, colnames(MethSM_GC))]
    MethSM_CG_nonNA_reads <- names(MethSM_CG_nonNA[!is.na(MethSM_CG_nonNA)])
    MethSM_GC_2 <- MethSM_GC_2[rownames(MethSM_GC_2) %in% MethSM_CG_nonNA_reads, ]
    
    MethSM_2 <- list(MethSM_GC_2, NULL, CG_OI)
    names(MethSM_2) <- c("GC", "CG", "CG_OI")
  }
  return(MethSM_2)
}

SortReads.CA.EK <- function (MethSM, range, BinsCoordinates, DE = NULL, sort = TRUE, WHAT = NULL) {
  MethSM_2    <- Extract.MethSM(MethSM, range, DE, WHAT)
  if(DE == FALSE){
    MethSM_CG   <- MethSM_2$CG[,grepl(MethSM_2$CG_OI, colnames(MethSM_2$CG))]
  }
  MethPattern <- as.factor(c(0,1))
  
  if (sort == TRUE){
    message("Collecting summarized methylation for bin")
    binMethylationList <- BinMethylation.EK(MethSM = MethSM_2$GC, Bin = BinsCoordinates)
    
    message("Splitting reads by pattern")
    if (length(binMethylationList) > 0) {
      sortedReadslist         <- list(names(binMethylationList[binMethylationList == MethPattern[1]]),
                                      names(binMethylationList[binMethylationList == MethPattern[2]]))
      names(sortedReadslist)  <- MethPattern
    } else {
      sortedReadslist = list()
    }
  } else if (sort == FALSE){
    sortedReadslist         <- list(names(MethSM_CG), c())
    names(sortedReadslist)  <- MethPattern
  }
  return(sortedReadslist)
}

SortReads.CA.CGme.EK <- function (MethSM, range, BinsCoordinates, DE = FALSE, accessibility = TRUE, WHAT = NULL) {
  MethSM_2      <- Extract.MethSM(MethSM, range, DE, WHAT)
  MethSM_CG     <- MethSM_2$CG[,grepl(MethSM_2$CG_OI, colnames(MethSM_2$CG))]
  MethPattern   <- as.factor(c(0,1))
  
  if(accessibility == TRUE){
    message("Collecting summarized methylation for bin")
    binMethylationList <- BinMethylation.EK(MethSM = MethSM_2$GC, Bin = BinsCoordinates)

    message("Splitting reads by pattern")
    if (length(binMethylationList) > 0) {
      sortedReadslist         <- list(names(binMethylationList[binMethylationList == MethPattern[1]]),
                                      names(binMethylationList[binMethylationList == MethPattern[2]]))
      names(sortedReadslist)  <- MethPattern
    } else {
      sortedReadslist         <- list()
    }
    
    message("Collecting summarized methylation for CGs in fractions")
    #CG unmethylated
    MethSM_CG_0               <- MethSM_CG[grepl(paste(sortedReadslist[[MethPattern[1]]], collapse = "|"), names(MethSM_CG))]
    sortedReadslist_0         <- list(names(MethSM_CG_0[MethSM_CG_0 == MethPattern[1]]),
                                      names(MethSM_CG_0[MethSM_CG_0 == MethPattern[2]]))
    names(sortedReadslist_0)  <- MethPattern
    
    #CG methylated
    MethSM_CG_1               <- MethSM_CG[grepl(paste(sortedReadslist[[MethPattern[2]]], collapse = "|"), names(MethSM_CG))]
    sortedReadslist_1         <- list(names(MethSM_CG_1[MethSM_CG_1 == MethPattern[2]]),
                                      names(MethSM_CG_1[MethSM_CG_1 == MethPattern[1]]))
    names(sortedReadslist_1)  <- rev(MethPattern)
    
    sortedReadslist_fs        <- list(sortedReadslist_0, sortedReadslist_1)
    names(sortedReadslist_fs) <- MethPattern
    
  } else if (accessibility == FALSE){
    message("Skipping accessibility sorting")
    sortedReadslist           <- list(names(MethSM_CG), c())
    names(sortedReadslist)    <- MethPattern
    
    message("Collecting summarized methylation for CGs in fractions")
    #CG unmethylated
    MethSM_CG_0               <- MethSM_CG[grepl(paste(sortedReadslist[[MethPattern[1]]], collapse = "|"), names(MethSM_CG))]
    sortedReadslist_0         <- list(names(MethSM_CG_0[MethSM_CG_0 == MethPattern[2]]),
                                      names(MethSM_CG_0[MethSM_CG_0 == MethPattern[1]]))
    names(sortedReadslist_0)  <- rev(MethPattern)
    sortedReadslist_fs        <- list(sortedReadslist_0, c())
    names(sortedReadslist_fs) <- MethPattern
  }
  return(sortedReadslist_fs)
}

SortReads.CA.plusTF.EK <- function (MethSM, BinsCoordinates, BinsCoordinates_TFBS = NULL, range, DE = FALSE, WHAT = NULL, plusCGme = TRUE, plusTF = TRUE) {
  MethSM_2      <- Extract.MethSM(MethSM, range, DE, WHAT)
  MethSM_CG     <- MethSM_2$CG[,grepl(MethSM_2$CG_OI, colnames(MethSM_2$CG))]
  MethPattern   <- as.factor(c(0,1))
  
  #Sort reads by chromatin accessibility
  message("Collecting summarized GC methylation for bin - for CA")
  binMethylationList <- BinMethylation.EK(MethSM = MethSM_2$GC, Bin = BinsCoordinates)
  
  message("Splitting reads by pattern - for CA")
  if (length(binMethylationList) > 0) {
    sortedReadslist         <- list(names(binMethylationList[binMethylationList == MethPattern[1]]),
                                    names(binMethylationList[binMethylationList == MethPattern[2]]))
    names(sortedReadslist)  <- MethPattern
  } else {
    sortedReadslist         <- list()
  }
  
  #Sort reads by TF binding
  message("Collecting summarized GC methylation for bins - for TFBS")
  binMethylationList2 <- lapply(seq_along(sortedReadslist), function(j){
    # j = 1
    lapply(seq_along(BinsCoordinates_TFBS), function(k) {
      # k = 3
      MethSM_2_sorted <- MethSM_2$GC[grepl(paste(sortedReadslist[[j]], collapse = "|"), rownames(MethSM_2$GC)),]
      BinMethylation.OLD(MethSM = MethSM_2_sorted, Bin = BinsCoordinates_TFBS[k])
    })
  })
  
  if(sum(is.na(binMethylationList2[[1]])) > 0){
    BinsCoordinates_TFBS <- IRanges(start = c(start(BinsCoordinates_TFBS[1])-10, start(BinsCoordinates_TFBS[2]), start(BinsCoordinates_TFBS[3])+10),
                                    end = c(end(BinsCoordinates_TFBS[1])-10, end(BinsCoordinates_TFBS[2]), end(BinsCoordinates_TFBS[3])+10))  
    
    message("Collecting summarized GC methylation for bins - for TFBS - round 2")
    binMethylationList2 <- lapply(seq_along(sortedReadslist), function(j){
      # j = 1
      lapply(seq_along(BinsCoordinates_TFBS), function(k) {
        # k = 3
        MethSM_2_sorted <- MethSM_2$GC[grepl(paste(sortedReadslist[[j]], collapse = "|"), rownames(MethSM_2$GC)),]
        BinMethylation.OLD(MethSM = MethSM_2_sorted, Bin = BinsCoordinates_TFBS[k])
      })
    })
  }
  if(sum(is.na(binMethylationList2[[1]])) > 0){
    message(paste0("!!!  At leats one bin overlaps with no covered Cytosines  !!! FINAL TERMINATION !!!"))
    return(NA)
  } else {
  
  message("Subsetting those reads that cover all bins - for TFBS")
  ReadsSubset <- lapply(seq_along(binMethylationList2), function(y) {
    # y = 1
    Reduce(intersect, lapply(seq_along(binMethylationList2[[y]]), function(x) {
      # x = 1
      names(binMethylationList2[[y]][[x]])
    }))})
  
  message("Summarizing reads into patterns - for TFBS")
  binMethylationList_subset <- lapply(seq_along(binMethylationList2), function(y) {
    lapply(seq_along(binMethylationList2[[y]]), function(x) {
      as.character(binMethylationList2[[y]][[x]][ReadsSubset[[y]]])
    })
  })
  
  MethPattern <- lapply(seq_along(binMethylationList_subset), function(y) {
    Reduce(paste0, binMethylationList_subset[[y]])
  })
  
  message("Splitting reads by pattern - for TFBS")
  sortedReadslist_l <- lapply(seq_along(ReadsSubset), function(y){
    if (length(ReadsSubset[[y]]) > 0) {
      sortedReadslist <- split(ReadsSubset[[y]], MethPattern[[y]])
    } else {
      sortedReadslist <- list()
    }
    
    if(plusCGme == TRUE){
      states = OneTFstates()
      states$closed = c(states$closed, states$unassigned)
      states = states[1:3]
      names(states) = c(states[[1]], states[[2]], states[[3]][1])
      
      OrderedReads <- sortedReadslist
      OrderedReads_states <- lapply(seq(length(states)), function(y){
        sR_l        <- lapply(seq(length(states[[y]])), function(x){
          sR        <- OrderedReads[as.character(states[[y]][x])]
          names(sR) <- names(states)[y]
          return(sR)
        })
        unlist(sR_l, use.names = F)
      })
      names(OrderedReads_states) <- names(states)
      
      message("Collecting summarized methylation for CGs in fractions - for TFBS")
      CGMethPattern <- c(1,0)
      
      sortedReadslist_fs <- lapply(seq(length(names(states))), function(x){
        # x=1
        OrderedReads_states_state   <- OrderedReads_states[[grep(as.character(names(states)[x]), names(OrderedReads_states))]]
        if(!is.null(OrderedReads_states_state)){
          # OrderedReads_states_state     <- OrderedReads_states[[grep(as.character(unlist(states)[x]), names(OrderedReads_states))]]
          MethSM_CG_state               <- MethSM_CG[grepl(paste(OrderedReads_states_state, collapse = "|"), names(MethSM_CG))]
          sortedReadslist_state         <- list(names(MethSM_CG_state[MethSM_CG_state == CGMethPattern[1]]),
                                                names(MethSM_CG_state[MethSM_CG_state == CGMethPattern[2]]))
          names(sortedReadslist_state)  <- CGMethPattern
        } else {
          sortedReadslist_state <- NULL
        } 
        return(sortedReadslist_state)
      })
      names(sortedReadslist_fs) <- names(states)
      return(sortedReadslist_fs)
    } else {
      return(sortedReadslist)
    }
  })
  return(sortedReadslist_l)
  }
}

SortReads.no.EK <- function (MethSM, range, BinsCoordinates, DE = FALSE, WHAT) {
  MethSM_2      <- Extract.MethSM(MethSM, range, DE, WHAT)
  MethSM_CG     <- MethSM_2$CG[,grepl(MethSM_2$CG_OI, colnames(MethSM_2$CG))]

  sortedReadslist <- list(names(MethSM_CG))
  return(sortedReadslist)
}

SortReads.TF.CGme.EK <- function (MethSM, BinsCoordinates, range, DE = FALSE, plusCGme = TRUE, WHAT = NULL) {
  MethSM_2      <- Extract.MethSM(MethSM, range, DE, WHAT)
  MethSM_CG     <- MethSM_2$CG[,grepl(MethSM_2$CG_OI, colnames(MethSM_2$CG))]
  
  if (DE == FALSE){
    MethSM_GC.sparse  <- MethSM[[2]][[1]][[1]]
  } else if (DE == TRUE){
    MethSM_GC.sparse  <- MethSM[[2]][[1]]
  }

#   MethPattern   <- as.factor(c(0,1))
#   MethSM_GC.sm  <- MethSM[[2]][[1]][[1]]
  
  message("Collecting summarized methylation for bins")
  binMethylationList <- lapply(seq_along(BinsCoordinates),
                               function(i) {
                                 BinMethylation(MethSM = MethSM_GC.sparse, Bin = BinsCoordinates[i])
                               })
  message("Subsetting those reads that cover all bins")
  ReadsSubset <- Reduce(intersect, lapply(binMethylationList,
                                         function(x) {
                                           names(x)
                                         }))
  message("Summarizing reads into patterns")
  binMethylationList_subset <- lapply(binMethylationList, function(x) {
    as.character(x[ReadsSubset])
  })
  MethPattern <- Reduce(paste0, binMethylationList_subset)
  message("Splitting reads by pattern")
  if (length(ReadsSubset) > 0) {
    sortedReadslist <- split(ReadsSubset, MethPattern)
  } else {
    sortedReadslist <- list()
  } 
  if(plusCGme == TRUE){
    states = OneTFstates()
    states$closed = c(states$closed, states$unassigned)
    states = states[1:3]
    names(states) = c(states[[1]], states[[2]], states[[3]][1])
    
    OrderedReads <- sortedReadslist
    OrderedReads_states <- lapply(seq(length(states)), function(y){
      sR_l        <- lapply(seq(length(states[[y]])), function(x){
        sR        <- OrderedReads[as.character(states[[y]][x])]
        names(sR) <- names(states)[y]
        return(sR)
      })
      unlist(sR_l, use.names = F)
    })
    names(OrderedReads_states) <- names(states)
    
    message("Collecting summarized methylation for CGs in fractions")
    CGMethPattern <- c(1,0)
    
    sortedReadslist_fs <- lapply(seq(length(states)), function(x){
      MethSM_CG_state               <- MethSM_CG[grepl(paste(OrderedReads_states[[names(states)[x]]], collapse = "|"), names(MethSM_CG))]
      sortedReadslist_state         <- list(names(MethSM_CG_state[MethSM_CG_state == CGMethPattern[1]]),
                                            names(MethSM_CG_state[MethSM_CG_state == CGMethPattern[2]]))
      names(sortedReadslist_state)  <- CGMethPattern
      return(sortedReadslist_state)
    })
    names(sortedReadslist_fs) <- names(states)
    return(sortedReadslist_fs)
  } else {
    return(sortedReadslist)
  }
}

#Plot 5mC Box -------------------------------------------------
Plot.5mC.Box.EK <- function (MethGR1, range, n_samples = NULL, DE1 = FALSE, 
                             MethGR2 = NULL, DE2 = FALSE, 
                             MethGR3 = NULL, DE3 = FALSE) {
  library(dplyr)
  library(ggplot2)
  
  if(n_samples > 1){
    
    # Compute plotting tibbles
    ## SAMPLE 1
    NAME_1 <- str_remove(names(elementMetadata(MethGR1[[1]]))[2], "_MethRate$") %>% 
      str_remove("SMF_MM_")
    if (DE1 == FALSE){ 
      S1_tbl <- MethGR1[[2]] %>% 
        as_tibble() %>% 
        dplyr::rename(me_mean = ends_with("_MethRate")) %>% 
        mutate(me_mean = round(me_mean*100)) %>% 
        mutate(y_coord_min = 0.2, y_coord_max = 1)
    } else {
      stop("Sample 1 has have CG methylation data")
    }
    
    if(n_samples > 1){
      ## SAMPLE 2
      if (DE2 == FALSE){ 
        S2_tbl <- MethGR2[[2]] %>% 
          as_tibble() %>% 
          dplyr::rename(me_mean = ends_with("_MethRate")) %>% 
          mutate(me_mean = round(me_mean*100)) %>% 
          mutate(y_coord_min = 1.2, y_coord_max = 2)
        NAME_2 <- str_remove(names(elementMetadata(MethGR2[[1]]))[2], "_MethRate$") %>% 
          str_remove("SMF_MM_")
      } else {
        S2_tbl <- S1_tbl %>% 
          mutate(me_mean = rep(0, nrow(S1_tbl))) %>% 
          mutate(y_coord_min = 1.2, y_coord_max = 2)
        NAME_2 <- str_remove(names(elementMetadata(MethGR2))[2], "_MethRate$") %>% 
          str_remove("SMF_MM_")
      }
      ## SAMPLE 3
      if(n_samples > 2){
        if (DE3 == FALSE){ 
          S3_tbl <- MethGR3[[2]] %>% 
            as_tibble() %>% 
            dplyr::rename(me_mean = ends_with("_MethRate")) %>% 
            mutate(me_mean = round(me_mean*100)) %>% 
            mutate(y_coord_min = 2.2, y_coord_max = 3)
          NAME_3 <- str_remove(names(elementMetadata(MethGR3[[1]]))[2], "_MethRate$") %>% 
            str_remove("SMF_MM_")
        } else {
          S3_tbl <- S1_tbl %>% 
            mutate(me_mean = rep(0, nrow(S1_tbl))) %>% 
            mutate(y_coord_min = 2.2, y_coord_max = 3)
          NAME_3 <- str_remove(names(elementMetadata(MethGR3))[2], "_MethRate$") %>% 
            str_remove("SMF_MM_")
        }
      }
      if(n_samples > 3){
        stop("Not more than 3 samples supported right now.")
      }
    }
  }
  # Plotting
  ggplot(mapping = aes(xmin = start-10, xmax = end+10, ymin = y_coord_min , ymax = y_coord_max, fill = me_mean)) +
    geom_rect(data = S1_tbl, color = "black") +
    (if(n_samples > 1){geom_rect(data = S2_tbl, color = "black")}) +
    (if(n_samples > 2){geom_rect(data = S3_tbl, color = "black")}) +
    geom_text(data = S1_tbl, aes(x = start-25, y = 0.6 , label =  paste0(me_mean, "%"))) +
    (if(n_samples > 1){geom_text(data = S2_tbl, aes(x = start-25, y = 1.6 , label =  paste0(me_mean, "%")))}) +
    (if(n_samples > 2){geom_text(data = S3_tbl, aes(x = start-25, y = 2.6 , label =  paste0(me_mean, "%")))}) +
    (if(n_samples == 1){scale_y_reverse(expand = expansion(add = 0.2), breaks = c(0.6), labels = NAME_1)}) +
    (if(n_samples == 2){scale_y_reverse(expand = expansion(add = 0.2), breaks = c(0.6, 1.6), labels = c(NAME_1, NAME_2))}) +
    (if(n_samples == 3){scale_y_reverse(expand = expansion(add = 0.2), breaks = c(0.6, 1.6, 2.6), labels = c(NAME_1, NAME_2, NAME_3))}) +
    scale_x_continuous(limits = c(start(range),end(range)), name = as.character(seqnames(range))) +
    scale_fill_gradientn(colours = rev(c(COLORS_METH, "white")), limits = c(0,100), breaks = seq(0,100,50)) +
    theme_bw() +
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.key.size = unit(0.25, "cm"),
          legend.title = element_blank(),
          legend.margin = margin(-2, 0, 4, 0),
          legend.box.margin = margin(-2, -5, 4, -5),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank()
    )
}
