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

	# If MATE_PAIR is set to true in the input file 
	if [ "$MATE_PAIR" == "true" ]; then
		# Check that output file doesn't exist then run cutAdapt, if it does print error message
		if [[ ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz && ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz || ! -z $FORCE ]]; then

			cutadapt -a $tTag -A `$SERA_PATH/bin/perlscript/reverseComplement.pl $fTag` -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz -p ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz --minimum-length 1 $RAWDATA_PE1 $RAWDATA_PE2 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;

			SuccessLog "${SAMPLEID}" "cutadapt -a $tTag -A `$SERA_PATH/bin/perlscript/reverseComplement.pl $fTag` -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz -p ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz --minimum-length 1 $RAWDATA_PE1 $RAWDATA_PE2 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;"

		else 
			ErrorLog "${SAMPLEID}" "${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz and ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz already exists and force was NOT used!";
		fi
	else
		if [ ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz ]; then

            cutadapt -a $tTag -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz --minimum-length 1 $RAWDATA_PE1 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
            SuccessLog "${SAMPLEID}" "cutadapt -a $tTag -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz --minimum-length 1 $RAWDATA_PE1 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;"

        else
            ErrorLog "${SAMPLEID}" "${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz already exists and force was NOT used!";
        fi
	fi
fi


if [ "$?" != "0" ]; then
	ErrorLog "${SAMPLEID}" "Failed in cutadapt...";
else
	SuccessLog "${SAMPLEID}" "Passed cutadapt";
fi		
