#!/bin/bash -l
#SBATCH -p core -n 1 
#SBATCH -t 24:00:00
##SBATCH -C fat
##SBATCH --qos=short

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/reads2Molecules" ]; then
	mkdir $ROOT_PATH/reads2Molecules;
fi
if [ ! -d "$ROOT_PATH/reads2Molecules/seqdata" ]; then
        mkdir $ROOT_PATH/reads2Molecules/seqdata;
fi

	SuccessLog "Counting reads per molecule for ${SAMPLEID}";

if [ $PLATFORM = "Illumina" ]; then

	if [ "$MATE_PAIR" == "true" ]; then
		if [[ -e ${ROOT_PATH}/filtereddata/${SAMPLEID}.read1.fastq.gz && -e ${ROOT_PATH}/filtereddata/${SAMPLEID}.read2.fastq.gz && -e ${ROOT_PATH}/filtereddata/${SAMPLEID}.index.fastq.gz ]]; then
			$ROOT_PATH_READS2MOLECULES/Reads2Molecules.sh "-t ${ROOT_PATH}/filtereddata/${SAMPLEID}.index.fastq.gz -p ${ROOT_PATH}/filtereddata/${SAMPLEID}.read1.fastq.gz -p2 ${ROOT_PATH}/filtereddata/${SAMPLEID}.read2.fastq.gz -o $ROOT_PATH/reads2Molecules/seqdata/${SAMPLEID}.reads2Molecules.read1.fastq -o2 $ROOT_PATH/reads2Molecules/seqdata/${SAMPLEID}.reads2Molecules.read2.fastq -bl 9 -tl 10";

			gzip $ROOT_PATH/reads2Molecules/${SAMPLEID}.reads2Molecules.read1.fastq;
			gzip $ROOT_PATH/reads2Molecules/${SAMPLEID}.reads2Molecules.read2.fastq;
		fi
	fi
fi


if [ "$?" != "0" ]; then
	ErrorLog "${SAMPLEID}" "Failed in cutadapt...";
else
#	basesWritten=`awk '{if(substr(\$0,1,7)=="# bases") {print \$4;}}' $ROOT_PATH/MosaikBuild/MosaikBuild_${SAMPLEID}_output.txt`;
#	SuccessLog "${SAMPLEID}" "Passed cutadapt ($extraLoginfo )  $basesWritten .";
	SuccessLog "${SAMPLEID}" "Passed cutadapt";
fi		
