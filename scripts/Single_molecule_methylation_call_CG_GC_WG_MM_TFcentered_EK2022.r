#TITLE:   WG MM TF single molecule methylation call for CpG and GpC
#AUTHOR:  Arnaud Krebs & Elisa Kreibich
#DATE:    27-10-2022
#VERSION: TFBS analysis strategy
#AIM:     Making the single molecule methylation call for CGs and GCs in a TF-motif-centered with at least one CpG fashion.
#         The output file is a multi-level list object with the following structure: [[C-context (1:GC|2:CG)]] [[TF bins (1:upstream | 2:middle | 3:downstream)]] [[TFBS]] [[binding state (000-111)]] vector of methylation values of single molecule reads.
#					VERSION USES JASPAR 2018 MOTIFS
	
# Load libraries and dependencies ------------------------------------------------------
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

	INPUT_QUASR = args[1]     #QuasR input file for qAlign (tab-delimted file with FileName (location of BAM file) and SampleName) | best only one sample/ replicate per input file | see example file
  WDW_SIZE    = args[2]     #101 bp for CA analysis used in Kreibich et al, 2023
  CHR1        = args[3]     #chromosome to start with (e.g. "1")
  CHR2        = args[4]     #chromosome to end with (e.g. "19")
  DIR_WD      = args[5]     #working directory (Github repository directory of Kreibich_2023_5mC_at_enhancers)
  INPUT_TFBS 	= args[6]     #input GRanges file for TFBS regions - from JASPAR 2018
	replicate 	= ""					#put "", when no selection of replicate, then the script runs through all samples and reps listed in the input file

## Define arguments when NOT running on the cluster
	# INPUT_QUASR = 'data/example/QuasR_aligned_files_ES_NO_R1_R2_R5a6_rmdup.txt'
	# WDW_SIZE    = 30
  # CHR1        = 1
  # CHR2        = 19
  # DIR_WD      = './Kreibich_2023_5mC_at_enhancers/' #location of this repository
  # INPUT_TFBS 	= 'data/GRange_mapped_jaspar2018_ChIP_score10_inBaits_BANP.rds'


# Load environment ---------------------------------------------------------------------
dir.create(path = paste0(DIR_WD, 'data_results/')
dir.create(path = paste0(DIR_WD, 'data_results/single_molecule_call/')
OUTDIR   <- paste0(DIR_WD, 'data_results/single_molecule_call/')


# Load aligned data ---------------------------------------------------------------------
Qproj = qAlign(sampleFile = INPUT_QUASR, 
                  genome = "BSgenome.Mmusculus.UCSC.mm10", 
                  paired = "fr", 
                  bisulfite = "undir")
Qproj@aligner = "Rbowtie"
Qaln = as.data.frame(alignments(Qproj)[[1]])
sampleNames = unique(alignments(Qproj)[[1]][,2])
if(!replicate == ""){
  sampleNames = sampleNames[grepl(replicate, sampleNames)]
}
message("Aligned data loaded")

# Load mapped motifs ---------------------------------------------------------------------
#ChIP validated updated JASPAR OBJECT 2018 - including BANP - only within Baits
MotDBf <- readRDS(INPUT_TFBS)
MotDBfs <- MotDBf#[MotDBf$name=='REST']
wmMapped_SOu <- MotDBfs

message("TF motifs loaded")

# Step 1: Create genomic chunks (based on TF positions) ---------------------------------------------------------------------

## Sort TFobject
wms <- (wmMapped_SOu[order(start(wmMapped_SOu))])
## Split by chr
wms.c <- split(wms, seqnames(wms))
regs <- unlist(GRangesList(lapply(seq_along(wms.c), function(i){
	x <- wms.c[[i]]
	l <- length(x)
	if(l>0){
	#cut by TF position
	seq.i <- seq(1,l,round(l/25))
	st.seq <- start(x[seq.i])
	lowB <- st.seq[1:(length(st.seq)-1)]-1000
	lowB[lowB<0] <- 1
	highB <- st.seq[2:(length(st.seq))]+1000
	highB[highB>seqlengths(Mmusculus)[as.character(seqnames(wms.c[[i]])[1])]] <- seqlengths(Mmusculus)[as.character(seqnames(wms.c[[i]])[1])]
	gr <- GRanges(seqnames(x[1]), IRanges(lowB, highB))
	gr
	}else(GRanges())
})))

message("TF motifs sorted")

		
# Step 2: Perform vectorial read extraction ---------------------------------------------------------------------

## Define the sample you want to analyse
sp.sbs <- seq_along(sampleNames)[seq_along(sampleNames)]
sampleNames[sp.sbs]

## Read extraction
for(sp in sp.sbs){
  # sp=1
  print(sampleNames[sp])
  for(CHR in seq(CHR1, CHR2)){
  # CHR = 19
  print(CHR)
	regsD <- cbind(as.data.frame(regs)[,1:3], sample = sampleNames[sp])
	regsD <- regsD[regsD[,1] %in% paste('chr', 1:19, sep = ''),]
	regsD <- regsD[regsD[,1] %in% paste('chr', CHR, sep = ''),]
	regsL <- split(regsD, 1:nrow(regsD))

	print(as.character(sampleNames[sp]))
	sRs <- mclapply(seq_along(regsL), function(i){
		# i=1
		# regDF=regsL[[i]]
		# target_range=wmMapped_SOu
		# projAll=Qproj
		if(WDW_SIZE == 14){
			sR <- extract.CG.GC.over.TF.states(regsL[[i]], Qproj, wmMapped_SOu, inMs=c(-7,7), upMs=c(-35,-25), doMs=c(25,35))  	
		} else if(WDW_SIZE == 30) {
			sR <- extract.CG.GC.over.TF.states(regsL[[i]], Qproj, wmMapped_SOu, inMs=c(-15,15), upMs=c(-35,-25), doMs=c(25,35))  				
		}
		sR
		}, mc.cores = 10)

	### Merge the genomic bins --> return a list organised by (1) context, (2) bin type
	sRs_t <- lapply(seq(2), function(context){
	            lapply(seq(3), function(bin){
	              sR_g <- lapply(seq(sRs), function(gbin){
	                sRs[[gbin]][[context]][[bin]]
	                })
								do.call(c, sR_g)
								})
	          })
				
	saveRDS(sRs_t, paste0(OUTDIR, 'tmp_', sampleNames[sp], '_CG_TF_jaspar2018_ChIP_score10_inBaits_BANP_sorted_met_vect_inM', WDW_SIZE, '_chr', CHR,'_NWCGW.rds'))
	message(paste(paste0(OUTDIR, 'tmp_', sampleNames[sp], '_CG_TF_jaspar2018_ChIP_score10_inBaits_BANP_sorted_met_vect_inM', WDW_SIZE, '_chr', CHR,'_NWCGW.rds'), 'saved'))
	}
}
message(Sys.time(), ': Script finished. Input arguments were:\nfile: ', INPUT_QUASR, '\nWDW_SIZE: ', WDW_SIZE, '\nchromosomes: ', CHR1, ' to ', CHR2)


