#!/bin/bash -l
#SBATCH -p devcore  -n 1
#SBATCH -t 01:00:00
##SBATCH --qos=short

module load bioinfo-tools

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog "${SAMPLEID}" "Start running cutadapt";

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/seqdata" ]; then
	mkdir $ROOT_PATH/seqdata;
fi

if [ $PLATFORM = "Illumina" ]; then

	# get sequencing tags
	. $SERA_PATH/config/sequencingTags.sh;

	# if file ending not fasta/fastq
	zcat $RAWDATA_PE1 > $SNIC_TMP/pe1.fastq;
	zcat $RAWDATA_PE2 > $SNIC_TMP/pe2.fastq;
	
	RAWDATA_PE1="$SNIC_TMP/pe1.fastq";
	RAWDATA_PE2="$SNIC_TMP/pe2.fastq";
	
	TEMP1_PE1="$SNIC_TMP/pe1.temp1.fastq";
	TEMP1_PE2="$SNIC_TMP/pe2.temp1.fastq";
	
	TEMP2_PE1="$SNIC_TMP/pe1.temp2.fastq";
	TEMP2_PE2="$SNIC_TMP/pe2.temp2.fastq";
	
	TEMP3_PE1="$SNIC_TMP/pe1.temp3.fastq";
	TEMP3_PE2="$SNIC_TMP/pe2.temp3.fastq";

	# If MATE_PAIR is set to true in the input file 
	if [ "$MATE_PAIR" == "true" ]; then
		# Check that output file doesn't exist then run cutAdapt, if it does print error message
		if [[ ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz && ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz || ! -z $FORCE ]]; then
			cutadapt -g file:${ROOT_PATH}/EGFR_path_5ptrim.fa -o $TEMP1_PE1 -p $TEMP1_PE2 $RAWDATA_PE1 $RAWDATA_PE2 --minimum-length 40 -e 0.12 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
			SuccessLog "${SAMPLEID}" "cutadapt -g file:${ROOT_PATH}/EGFR_path_5ptrim.fa -o $TEMP1_PE1 -p $TEMP1_PE2 $RAWDATA_PE1 $RAWDATA_PE2 --minimum-length 40 -e 0.12 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;";
			cutadapt -g file:${ROOT_PATH}/EGFR_path_5ptrim.fa -o $TEMP2_PE1 -p $TEMP2_PE2 $TEMP1_PE2 $TEMP1_PE1 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
			SuccessLog "${SAMPLEID}" "cutadapt -g file:${ROOT_PATH}/EGFR_path_5ptrim.fa -o $TEMP2_PE1 -p $TEMP2_PE2 $TEMP1_PE2 $TEMP1_PE1 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;";
			cutadapt -a file:${ROOT_PATH}/EGFR_path_3ptrim.fa -o $TEMP3_PE1 -p $TEMP3_PE2 $TEMP2_PE2 $TEMP2_PE1 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
			SuccessLog "${SAMPLEID}" "cutadapt -a file:${ROOT_PATH}/EGFR_path_3ptrim.fa -o $TEMP3_PE1 -p $TEMP3_PE2 $TEMP2_PE2 $TEMP2_PE1 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;";
			cutadapt -a file:${ROOT_PATH}/EGFR_path_3ptrim.fa -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz -p ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz $TEMP3_PE2 $TEMP3_PE1 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
			SuccessLog "${SAMPLEID}" "cutadapt -a file:${ROOT_PATH}/EGFR_path_3ptrim.fa -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz -p ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz $TEMP3_PE2 $TEMP3_PE1 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;";
		else 
			ErrorLog "${SAMPLEID}" "${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz and ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz already exists and force was NOT used!";
		fi
	else
		ErrorLog "${SAMPLEID}" "Only implemented for paired-end sequencing!";
	fi
fi


if [ "$?" != "0" ]; then
	ErrorLog "${SAMPLEID}" "Failed in cutadapt...";
else
	SuccessLog "${SAMPLEID}" "Passed cutadapt";
fi		
