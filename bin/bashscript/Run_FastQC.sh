#!/bin/bash
#
# Script for running MosaikBuild on either SOLiD or Illumina samples.
#
#SBATCH -p devcore  -n 1
#SBATCH -t 30:00
##SBATCH --qos=short

# Include functions
. $SERA_PATH/includes/logging.sh

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/FastQC" ]; then
	mkdir $ROOT_PATH/FastQC;
fi

if [ "$MATE_PAIR" == "true" ]; then
	extraLoginfo="mate-pair/pair-ends";
else
	extraLoginfo="single";
fi

SuccessLog "${SAMPLEID}" "Building $extraLoginfo reads Mosaik input file...";

# Check for cutadapt sequences
if [ -e ${ROOT_PATH}/filtereddata/${SAMPLEID}.read1.fastq.gz ]; then
	PE1=${ROOT_PATH}/filtereddata/${SAMPLEID}.read1.fastq.gz;
else
	PE1=$RAWDATA_PE1;
fi
if [ -e ${ROOT_PATH}/filtereddata/${SAMPLEID}.read2.fastq.gz ]; then
	PE2=${ROOT_PATH}/filtereddata/${SAMPLEID}.read2.fastq.gz;
else
	PE2=$RAWDATA_PE2;
fi

SuccessLog "${SAMPLEID}" "Using $PE1 as read1 and $PE2 as read2 as input to FastQC...";

# Run FastQC
# If platform is Illumina
if [ $PLATFORM = "Illumina" ]; then
	# Check that the output file doesn't exist or if force is given
	if [[ ! -e $ROOT_PATH/FastQC/${SAMPLEID}.read1_fastqc.html || ! -e $ROOT_PATH/FastQC/${SAMPLEID}.read2_fastqc.html || ! -z $FORCE ]]; then

		# Do paired end library or not
		if [ "$MATE_PAIR" == "true" ]; then
			# Check that input files exist, if not print error message
			if [[ -e ${PE1} && ${PE2} ]]; then
				fastqc -o FastQC $PE1 $PE2 > $ROOT_PATH/FastQC/fastqc_${SAMPLEID}_output.txt;
			else
				ErrorLog "${SAMPLEID}" "${PE1} and/or ${PE2} do NOT exist!";
			fi
		
		# If single reads are used only care about the first read
		else
			# Check that input file exists, if not print error message
			if [ -e ${PE1} ]; then
				fastqc -o FastQC $PE1 > $ROOT_PATH/FastQC/fastqc_${SAMPLEID}_output.txt;
			else 
				ErrorLog "${SAMPLEID}" "${PE1} does NOT exist!";
			fi
		fi
	else
		ErrorLog "${SAMPLEID}" "FastQC file already existed and force was NOT used!"; 

	fi

else
	ErrorLog "${SAMPLEID}" "Platform $PLATFORM not recognized (should be Illumina).";
	exit 1;
fi
	
# Check if MosaikBuild worked
if [ "$?" != "0" ]; then
	ErrorLog "${SAMPLEID}" "Failed in FastQC...";
else
	SuccessLog "${SAMPLEID}" "Passed FastQC ($extraLoginfo reads)";
fi
