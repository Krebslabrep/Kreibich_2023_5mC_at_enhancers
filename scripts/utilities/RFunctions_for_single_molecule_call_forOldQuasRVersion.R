#Useful function -------------------------------------------------------------------------
string.split <- function(string, sep, pos){
  unlist(lapply(string, function(x){lapply(strsplit(x, sep), '[', pos)}))
  }



# Extract methylation distribution for single CGs methylation as a function of GC methylation in a window ----------------------

# regDF = regsL[[i]]
# sampleSheet = "/g/krebs/krebs/HTS/SMF/MM/QuasR_input_deduplicated.txt"
# target_range = col_wind

extract_CG_GC_overs_CG_bins <- function(regDF, sampleSheet, target_range){
  #sort reads around CG collection windows
  #requires information in all three bins (more stringent)
  
  #reload libraries
  library(QuasR)
  library("BSgenome.Mmusculus.UCSC.mm10")
  # source('/g/krebs/krebs/Rscripts/Rfunctions/SMF/Single_molecule_manipulation _function.r')
  # source('/g/krebs/krebs/Rscripts/Rfunctions/useful_functionsV1.r')
  
  projAll <- qAlign(sampleSheet, "BSgenome.Mmusculus.UCSC.mm10", paired = "fr", bisulfite = "dir")
  projAll@aligner <- "Rbowtie"
  proj <- projAll[alignments(projAll)$genome$SampleName == regDF[1,4]]
  
  ## Define regions to extract methylation information
  reg <- GRanges(seqnames = Rle(regDF[1,1]), ranges = IRanges(regDF[1,2], end = regDF[1,3]), seqlengths = seqlengths(Mmusculus))
  
  ## Expand the region at the end to catch all the GCs within the window
  regExp <- reg
  # regExp <- resize(regExp, width = width(regExp) + 2000, fix = "start")
  end(regExp) <- pmin(end(regExp) + 2000, seqlengths(regExp)[as.character(seqnames(regExp))])

  ## Map CG excluding all ambiguous contexts
  dm3_DCGHNs_chr_XSV <- matchPattern(DNAString("NWCGW"), Mmusculus[[as.character(seqnames(reg))]], fixed = "subject")
  dm3_DCGHNs_chr <- GRanges(seqnames = Rle(as.character(seqnames(reg))), ranges = IRanges(start(dm3_DCGHNs_chr_XSV), end = end(dm3_DCGHNs_chr_XSV)), strand = "+")
  
  dm3_DCGHNs_chr <- resize(dm3_DCGHNs_chr, 1, fix = "center")
  dm3_DCGHNs_chr_coord <- start(dm3_DCGHNs_chr)
  
  ## Map GC excluding all ambiguous contexts			  
  dm3_DGCHNs_chr_XSV <- matchPattern(DNAString("DGCHN"), Mmusculus[[as.character(seqnames(reg))]], fixed = "subject")
  dm3_DGCHNs_chr <- GRanges(seqnames = Rle(as.character(seqnames(reg))), ranges = IRanges(start(dm3_DGCHNs_chr_XSV), end = end(dm3_DGCHNs_chr_XSV)), strand = "+")
  dm3_DGCHNs_chr <- resize(dm3_DGCHNs_chr, 1, fix = "center")	
  dm3_DGCHNs_chr_coord <- start(dm3_DGCHNs_chr)

  ## Call methylation vectors
  methAln <- qMeth(proj = proj, query = regExp, mode = "allC", reportLevel = "alignment")[[1]]
  
  ## Collapse strands according to context
  methAln_CG <- methAln
  methAln_GC <- methAln	
  #GCs
  methAln_GC$Cid[methAln_GC$strand == "-"] <- methAln_GC$Cid[methAln_GC$strand == "-"] + 1
  selCs_GC <- methAln_GC$Cid %in% dm3_DGCHNs_chr_coord
  #CGs
  methAln_CG$Cid[methAln_CG$strand == "-"] <- methAln_CG$Cid[methAln_CG$strand == "-"] - 1
  selCs_CG <- methAln_CG$Cid %in% dm3_DCGHNs_chr_coord 
  
  ## Generate state separated methylation vectors
  regT <- GRanges(regDF[1,1], IRanges(regDF[1,2], regDF[1,3]))
  TFs <- subsetByOverlaps(target_range, regT, type = 'within')
  
  if(length(methAln[[1]]) > 0){ #control that at least one read is found in the region
    sRs = single_molecule_CG_GC_met_CG_bins(methAln_GC, methAln_CG, selCs_GC, selCs_CG, st=TFs)
    sRs
  }else{
    list()
  }
}


single_molecule_CG_GC_met_CG_bins <- function(methAln_GC, methAln_CG, selCs_GC, selCs_CG, st){
  ## Create  C-ranges
  GCrange <- IRanges(methAln_GC[[2]][selCs_GC], methAln_GC[[2]][selCs_GC])
  CGrange <- IRanges(methAln_CG[[2]][selCs_CG], methAln_CG[[2]][selCs_CG])
  
  bins <- list(st)
  
  ## Find overlaping Cs
  #GCs
  ovsGC <- lapply(seq_along(bins), function(i){
    as.matrix(findOverlaps(ranges(bins[[i]]), GCrange))
  })
  #CGs
  #limit to the central CG
  ovsCG <- lapply(seq_along(bins), function(i){
    as.matrix(findOverlaps(ranges(resize(bins[[i]], 1, fix = "center")), CGrange))
  })
  
  ## Compute the methylation vectors
  #GCs
  MvsGC <- lapply(seq_along(bins), function(i){
    methAln_GC[[4]][selCs_GC][ovsGC[[i]][,2]]
  })
  #CGs
  MvsCG <- lapply(seq_along(bins), function(i){
    methAln_CG[[4]][selCs_CG][ovsCG[[i]][,2]]
  })
  
  ## Compute the readIDs
  #GCs					
  GC_IDvs=lapply(seq_along(bins),function(i){
    paste(names(bins[[1]])[ovsGC[[i]][,1]],methAln_GC[[1]][selCs_GC][ovsGC[[i]][,2]],sep='_')
  })	
  
  #CGs					
  CG_IDvs <- lapply(seq_along(bins), function(i){
    paste(names(bins[[1]])[ovsCG[[i]][,1]],methAln_CG[[1]][selCs_CG][ovsCG[[i]][,2]],sep='_')
  })	
  
  ## Group GCs per region/ per read
  #GCs
  g.Ms <- lapply(seq_along(bins), function(i){
    round(tapply(MvsGC[[i]], GC_IDvs[[i]], mean))
  })
  
  ## Calculate average methylation per molecule in each bin
  #GCs
  g.MsGCs <- lapply(seq_along(bins), function(i){
    tapply(MvsGC[[i]], GC_IDvs[[i]], mean)
  })
  #CGs
  g.MsCGs <- lapply(seq_along(bins), function(i){
    tapply(MvsCG[[i]],CG_IDvs[[i]],mean)
  })
  
  ## Find unique IDs
  u.IDs <- lapply(seq_along(bins), function(i){
    sort(unique(GC_IDvs[[i]]))
  })
  
  intMat <- do.call(cbind, lapply(seq_along(bins), function(i){
    u.IDs[[1]] %in% u.IDs[[i]]
  }))
  
  ## Intersect IDs
  ids <- seq_along(u.IDs[[1]])[rowSums(intMat) == ncol(intMat)]

  ## Binary methylation vectors
  sID <- u.IDs[[1]][ids]
  
  s.M <- lapply(seq_along(bins), function(i){
    g.Ms[[i]][sID]
  })
  
  s.M.GC <- lapply(seq_along(bins), function(i){
    g.MsGCs[[i]][sID]
  })
  
  ## Make sure to consider only CGs on sID reads (reads that are classified)
  s.M.CG <- lapply(seq_along(bins), function(i){
    #create an empty array 
    vect <- rep(NA, length(sID))
    names(vect) <- sID
    
    vals <- g.MsCGs[[i]][names(g.MsCGs[[i]] ) %in% sID]
    vect[names(vals)] <- vals
    vect
  })
  
  ## Make binary classification accessible/inaccessible
  patternMat <- do.call(cbind, s.M)
  pattern <- apply(patternMat, 1, function(x) {paste(as.character(x), sep = '', collapse = '')})
  st.id <- paste(string.split(sID, '_', 1), string.split(sID, '_', 2), sep = '_')
  read.id <- string.split(sID, '_', 3)
  
  
  ## SORTED READ LISTS
  if(length(ids) > 0){
    #Split reads
    sR <- split(read.id, paste(st.id, pattern, sep = '_'))
    sRl <- split(sR, paste(string.split(names(sR), '_', 1), string.split(names(sR), '_', 2), sep = '_'))
    
    sR.GC <- lapply(seq_along(bins), function(i){split(s.M.GC[[i]], paste(st.id, pattern, sep = '_'))})
    sRl.GC <- lapply(seq_along(bins), function(i){split(sR.GC[[i]], paste(string.split(names(sR), '_', 1), string.split(names(sR), '_', 2), sep = '_'))})
    
    sR.CG <- lapply(seq_along(bins), function(i){split(s.M.CG[[i]], paste(st.id, pattern, sep = '_'))})
    sRl.CG <- lapply(seq_along(bins), function(i){split(sR.CG[[i]], paste(string.split(names(sR), '_', 1), string.split(names(sR), '_', 2), sep = '_'))})
    
    output <- list(metGC = sRl.GC, metCG = sRl.CG)
    output
  } else {
    list(metGC = list(), metCG = list())
  }
}

