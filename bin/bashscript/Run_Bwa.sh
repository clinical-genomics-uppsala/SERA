#!/bin/bash
#
# Script running bwa
#
#SBATCH -p devcore  -n 3
#SBATCH -t 02:00:00
##SBATCH --qos=short

# Include functions
. $SERA_PATH/includes/logging.sh

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/Bwa" ]; then
	mkdir $ROOT_PATH/Bwa;
fi

SuccessLog "${SAMPLEID}" "Aligning paired-end reads with bwa mem";

# Check for cutadapt sequences
if [ -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz ]; then
	PE1=${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz;
else
	PE1=$RAWDATA_PE1;
fi
if [ -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz ]; then
	PE2=${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz;
else
	PE2=$RAWDATA_PE2;
fi

SuccessLog "${SAMPLEID}" "Using $PE1 as read1 and $PE2 as read2 as input to bwa...";

# Run bwa
# If platform is Illumina
if [ $PLATFORM = "Illumina" ]; then
	# Check that the output file doesn't exist or if force is given
	if [[ ! -e $ROOT_PATH/Bwa/${SAMPLEID}.sam || ! -z $FORCE ]]; then
		# Check that input files exist
		if [ "$MATE_PAIR" == "true" ]; then
            if [[ -e ${PE1} && ${PE2} ]]; then
                # Get the date when the analysis is run
                now=$('date' +"%Y%m%d")
                bwa mem -M -R "@RG\tID:"$now"_${SAMPLEID}\tSM:${SAMPLEID}\tPL:illumina" ${GENOME_REF} ${PE1} ${PE2} -t 3 | samtools view -bS /dev/stdin | samtools sort -@ 3 /dev/stdin $ROOT_PATH/Bwa/${SAMPLEID}.sorted;
                samtools index $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam;
                samtools flagstat $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam > $ROOT_PATH/Bwa/${SAMPLEID}.alignmentStats.txt;


                SuccessLog "${SAMPLEID}" "bwa mem -M -R \"@RG\tID:"$now"_${SAMPLEID}\tSM:${SAMPLEID}\tPL:illumina\" ${GENOME_REF} ${PE1} ${PE2} -t 3 | samtools view -bS /dev/stdin | samtools sort -@ 3 /dev/stdin $ROOT_PATH/Bwa/${SAMPLEID}.sorted;"
                SuccessLog "${SAMPLEID}" samtools flagstat $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam;

            else
                ErrorLog "${SAMPLEID}" "${PE1} and/or ${PE2} do NOT exist!";
            fi
        else
            if [[ -e ${PE1} ]]; then
                # Get the date when the analysis is run
                now=$('date' +"%Y%m%d")
                bwa mem -M -R "@RG\tID:"$now"_${SAMPLEID}\tSM:${SAMPLEID}\tPL:illumina" ${GENOME_REF} ${PE1} -t 3 | samtools view -bS /dev/stdin | samtools sort -@ 3 /dev/stdin $ROOT_PATH/Bwa/${SAMPLEID}.sorted;
                samtools index $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam;
                samtools flagstat $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam > $ROOT_PATH/Bwa/${SAMPLEID}.alignmentStats.txt;


                SuccessLog "${SAMPLEID}" "bwa mem -M -R \"@RG\tID:"$now"_${SAMPLEID}\tSM:${SAMPLEID}\tPL:illumina\" ${GENOME_REF} ${PE1} -t 3 | samtools view -bS /dev/stdin | samtools sort -@ 3 /dev/stdin $ROOT_PATH/Bwa/${SAMPLEID}.sorted;"
                SuccessLog "${SAMPLEID}" samtools flagstat $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam;
            else
                ErrorLog "${SAMPLEID}" "${PE1} do NOT exist!";
            fi
        fi
	else
		ErrorLog "${SAMPLEID}" "$ROOT_PATH/Bwa/${SAMPLEID}.sam already exists and force was NOT used!"; 

	fi

else
	ErrorLog "${SAMPLEID}" "Platform $PLATFORM not recognized (only interpreted for Illumina).";
	exit 1;
fi
	
# Check if bwa worked
if [ "$?" != "0" ]; then
	ErrorLog "${SAMPLEID}" "Failed in bwa alignment...";
else
	SuccessLog "${SAMPLEID}" "Passed bwa alignment";
fi
