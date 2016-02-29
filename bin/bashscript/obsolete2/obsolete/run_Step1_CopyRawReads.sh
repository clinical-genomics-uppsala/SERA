#!/bin/bash
#
# Script copys and /usr/bin/pigz -p 1 samples into analysis directory.
#
#SBATCH -p node -n 1
#SBATCH -t 01:00:00

# Include functions
. $SERA_PATH/includes/logging.sh

# Copy Illumina read functionecho "Mate pair: $MATE_PAIR";
function copyIlluminaReads {

	if [ "$(file -ib ${1})" == "application/x-gzip; charset=binary" ]; then
		cp ${1}* $ROOT_PATH/Reads/${SAMPLEID}_${2}.fastq.gz;
	elif [ "$(file -ib ${1})" == "text/plain; charset=us-ascii" ]; then
		cat ${1}* | gzip > $ROOT_PATH/Reads/${SAMPLEID}_${2}.fastq.gz;
	else
		ErrorLog "${SAMPLEID}" "Unknown file format, not gzip nor ascii, no copy done.";
		exit 1;
	fi
	
	# Check if copy worked
	if [ "$?" != "0" ]; then
		ErrorLog "${SAMPLEID}" "Mate-${2} failed in copy raw reads (fastq-file).";
	else
		SuccessLog "${SAMPLEID}" "Mate-${2} passed copy raw reads (fastq-file).";
	fi
}

# Copy SOLiD reads function
function copySolidReads {

	if [ ! -e $ROOT_PATH/Reads/${SAMPLEID}.csfasta.gz ]; then
		
		cat ${1}*.csfasta | gzip 1 > $ROOT_PATH/Reads/${SAMPLEID}_${2}.csfasta.gz;
	
		# Check if copy worked
		if [ "$?" != "0" ]; then
			ErrorLog "${SAMPLEID}" "Failed in copy raw reads (csfasta-file).";
		else
			SuccessLog "${SAMPLEID}" "Passed copy raw reads (csfasta-file).";
		fi
	fi

	if [ ! -e  $ROOT_PATH/Reads/${SAMPLEID}.qual.gz ]; then

		cat ${1}*.qual | gzip > $ROOT_PATH/Reads/${SAMPLEID}_${2}.qual.gz;

		# Check if copy worked
		if [ "$?" != "0" ]; then
			ErrorLog "${SAMPLEID}" "Failed in copy raw reads (qual-file).";
		else
			SuccessLog "${SAMPLEID}" "Passed copy raw reads (qual-file).";
		fi
	fi
}

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/Reads" ]; then
	mkdir $ROOT_PATH/Reads;
fi

# ReadsNoDot run only if the platform is SOLiD
if [ $PLATFORM = "SOLiD" ]; then
	
	copySolidReads $RAWDATA_PE1 1;
	
	if [ "$RAWDATA_PE2" != "false" ]; then
		copySolidReads $RAWDATA_PE2 2;
	fi
fi

# ReadsNoDot run only if the platform is SOLiD
if [ $PLATFORM = "Illumina" ]; then

	if [ ! -e $ROOT_PATH/Reads/${SAMPLEID}_1.fastq.gz ]; then
		copyIlluminaReads $RAWDATA_PE1 1;		
	fi

	if [ "$RAWDATA_PE2" != "false" ]; then
		if [ ! -e $ROOT_PATH/Reads/${SAMPLEID}_2.fastq.gz ]; then
			copyIlluminaReads $RAWDATA_PE2 2;
		fi

		# Simulate single end reads?		
		if [ "$MATE_PAIR" == "false" ]; then
			mv $ROOT_PATH/Reads/${SAMPLEID}_1.fastq.gz $ROOT_PATH/Reads/${SAMPLEID}_1.old.fastq.gz;
			mv $ROOT_PATH/Reads/${SAMPLEID}_2.fastq.gz $ROOT_PATH/Reads/${SAMPLEID}_2.old.fastq.gz;
			zcat $ROOT_PATH/Reads/${SAMPLEID}_1.old.fastq.gz $ROOT_PATH/Reads/${SAMPLEID}_2.old.fastq.gz | gzip > $ROOT_PATH/Reads/${SAMPLEID}_1.fastq.gz;
			rm $ROOT_PATH/Reads/${SAMPLEID}_1.old.fastq.gz;
			rm $ROOT_PATH/Reads/${SAMPLEID}_2.old.fastq.gz;
		fi	
	fi

fi
