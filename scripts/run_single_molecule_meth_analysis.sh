SCRIPT_R_NAME="/g/krebs/kreibich/GitHub_Kreibich_et_al/Single_molecule_methylation_analysis_CG_GC_WG_MM_CGcentered_EK2022.r"	#CA analysis
# SCRIPT_R_NAME="/g/krebs/kreibich/GitHub_Kreibich_et_al/Single_molecule_methylation_analysis_CG_GC_WG_MM_TFcentered_EK2022.r"	#TFBS analysis

Rscript="/g/funcgen/bin/Rscript-4.1.2"
LOGDIR="/scratch/kreibich/SLURM_LOGS"
mkdir -p $LOGDIR

NCORES="1" 		#nb cores to be used

# Arguments needed:
	# NAME        = args[1]  #celltype name (e.g. "ES", if sample name is "ES_R1_R2", or "ES_NO" if sample name is "ES_NO_R1_R2")
	# REP         = args[2]  #replicates (e.g. "R1_R2", if sample name is "ES_R1_R2")
	# CHR1        = args[3]  #chromosome to start with (e.g. "1")
	# CHR2        = args[4]  #chromosome to end with (e.g. "19")
	# PERFORM_CMH = args[5]  #logical argument whether to perform the statistical test ("TRUE"). If "FALSE" only coverage table will be made.
	# WDW_SIZE    = args[6]  #window size of the genomic bin (e.g. "101" for CA analysis, hard coded 30 bp for TFBS)
	# DIR_WD      = args[7]  #working directory (e.g. './Keibich_2022/')
	# DIR_QUASR   = args[8]  #input directory of QuasR files for qAlign (e.g. '/g/krebs/kreibich/HTS/SMF/MM/QuasR/QuasR_rmdup/')


#Launching multiple samples:
NAME="ES_NO" 			#for sample_file <- paste0('/g/krebs/kreibich/HTS/SMF/MM/QuasR/QuasR_rmdup/QuasR_aligned_files_', NAME, '_', REP, '_rmdup.txt')
REP="R1_R2_R5a6" 		#for sample_file <- paste0('/g/krebs/kreibich/HTS/SMF/MM/QuasR/QuasR_rmdup/QuasR_aligned_files_', NAME, '_', REP, '_rmdup.txt')
CHR1=$1
CHR2=$2
PERFORM_CMH=TRUE
WDW_SIZE=101
DIR_WD='/g/krebs/kreibich/GitHub_Kreibich_et_al/'
DIR_QUASR='/g/krebs/kreibich/HTS/SMF/MM/QuasR/QuasR_rmdup/'

#If launching all chromosomes individually the same time
# for CHR in $(seq 1 19)
#  	do
# 	CMD="${Rscript} ${SCRIPT_R_NAME} ${NAME} ${REP} ${CHR} ${CHR} ${PERFORM_CMH} ${WDW_SIZE} ${DIR_WD} ${DIR_QUASR}"
# 	echo $CMD
# 	sbatch --wrap "${CMD}" --job-name="SMana" -t 6:00:00  -N 1 -n ${NCORES} --mem=100G --mail-type=END,FAIL --mail-user=elisa.kreibich@embl.de -o "${LOGDIR}/slurm.%j.out" -e "${LOGDIR}/slurm.%j.err" 
# 	done

#If launching batches of chromosomes
	CMD="${Rscript} ${SCRIPT_R_NAME} ${NAME} ${REP} ${CHR1} ${CHR2} ${PERFORM_CMH} ${WDW_SIZE} ${DIR_WD} ${DIR_QUASR}"
	sbatch --wrap "${CMD}" --job-name="SMana" -t 6:00:00  -N 1 -n ${NCORES} --mem=100G --mail-type=END,FAIL --mail-user=elisa.kreibich@embl.de -o "${LOGDIR}/slurm.%j.out" -e "${LOGDIR}/slurm.%j.err" 

