## THIS FILE IS: scripts/run_single_molecule_meth_call.sh
SCRIPT_R_NAME="scripts/Single_molecule_methylation_call_CG_GC_WG_MM_CGcentered_EK2022.r"	#CA analysis call
# SCRIPT_R_NAME="scripts/Single_molecule_methylation_call_CG_GC_WG_MM_TFcentered_EK2022.r"	#TFBS analysis call

Rscript="/g/funcgen/bin/Rscript-4.1.2"
LOGDIR="/scratch/kreibich/SLURM_LOGS"
mkdir -p $LOGDIR

NCORES="2" #nb cores to be used

# Arguments needed:
	# INPUT_QUASR = args[1]     #QuasR input file for qAlign (tab-delimted file with FileName (location of BAM file) and SampleName) | best only one sample/ replicate per input file
	# WDW_SIZE    = args[2]     #101 bp for CA analysis used in Kreibich et al, 2023
	# CHR1        = args[3]     #chromosome to start with (e.g. "1")
	# CHR2        = args[4]     #chromosome to end with (e.g. "19")
	# DIR_WD      = args[5]     #working directory (e.g. './Keibich_2022/')
	# INPUT_TFBS  = args[7]     #input GRanges file for TFBS regions - from JASPAR 2018
	# replicate   = ""			#put "", when no selection of replicate, then the script runs through all samples and reps listed in the input file

# #Launching multiple samples:
TABLE_DIR="data/example/QuasR_files/" 	#QuasR input files for BAM files of deduplicated samples
WDW_SIZE=101 #101 | 30  - Kreibich et al used: 101 bp window for CG centric analysis and 30 bp window for TFBS centric analysis
CHR1=$1
CHR2=$2
DIR_WD='Kreibich_2023_5mC_at_enhancers/' #This repsitory
DIR_QUASR='data/example/QuasR_aligned_files_ES_NO_R1_R2_R5a6_rmdup.txt' #QuasR input file for qAlign 
INPUT_ADD='data/MouseMethyl_Bait_merged_mm10.bed' #Baits bed file for CA analysis
# INPUT_ADD='data/GRange_mapped_jaspar2018_ChIP_score10_inBaits_BANP.rds' #TFBS motifs GRange for TFBS analysis

#Launching multiple samples:
for INPUT_QUASR in $(find $TABLE_DIR -name "QuasR_input_*.txt")
 	do
	CMD="${Rscript} ${SCRIPT_R_NAME} ${INPUT_QUASR} ${WDW_SIZE} ${CHR1} ${CHR2} ${DIR_WD} ${INPUT_ADD}"
	echo $CMD 
    
	sbatch --wrap "${CMD}" --job-name="SMcall" -t 1440:00  -N 1 -n ${NCORES} --mem=100G --mail-type=END,FAIL --mail-user=elisa.kreibich@embl.de -o "${LOGDIR}/slurm.%j.out" -e "${LOGDIR}/slurm.%j.err" 
	done
	

#Launching individual sample:
# INPUT_QUASR="data/example/QuasR_aligned_files_ES_NO_R1_R2_R5a6_rmdup.txt"
# CMD="${Rscript} ${SCRIPT_R_NAME} ${INPUT_QUASR} ${WDW_SIZE} ${CHR1} ${CHR2} ${DIR_WD} ${INPUT_ADD}"
# sbatch --wrap "${CMD}" --job-name="SMcall" -t 1440:00  -N 1 -n ${NCORES} --mem=100G --mail-type=END,FAIL --mail-user=elisa.kreibich@embl.de -o "${LOGDIR}/slurm.%j.out" -e "${LOGDIR}/slurm.%j.err" 
