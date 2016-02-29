#!/bin/bash
#
# Script to export aligment archives as BAM-files.
#
#SBATCH -p devcore  -n 1
#SBATCH -t 01:00:00
##SBATCH --qos=short
# Include functions
. $SERA_PATH/includes/logging.sh;


# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/SamBamFiles" ]; then
	mkdir $ROOT_PATH/SamBamFiles;
fi


$ROOT_PATH/Bwa/${SAMPLEID}.sorted

SuccessLog "${SAMPLEID}" "Exporting BAM files...";
# Create SNPmania header file for each REFSEQ
if [[ ! -e $ROOT_PATH/refFiles/${REFSEQ}.SNPmania.header.sam  && -e $ROOT_PATH/Alignments/${SAMPLEID}.h.sapiens.aligned.bam ]]; then
	samtools view -H $ROOT_PATH/Alignments/${SAMPLEID}.h.sapiens.aligned.bam > $ROOT_PATH/refFiles/${REFSEQ}.SNPmania.header.sam;
fi	


# Check if call type is h.sapiens
if [ ${DESIGN_TYPE} == "PCR" ]; then
	if [ ${CALL_TYPE} == "h.sapiens" ]; then
		# Check if the sorted bam file doesn't exist or if force is used then create it, otherwise 
		if [[ ! -e $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.bam || ! -z $FORCE ]]; then
			samtools reheader $ROOT_PATH/SamBamFiles/${SAMPLEID}.header.sam $ROOT_PATH/Alignments/${SAMPLEID}.h.sapiens.aligned.bam | samtools sort /dev/stdin $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted;
			#Index file
			if [ -e $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.bam ]; then
				samtools index $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.bam $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.bam.bai;
			else
				ErrorLog "Could NOT create "$ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.bam"!";
        	fi
		# Remove the file with header for the sample
		rm $ROOT_PATH/SamBamFiles/${SAMPLEID}.header.sam;
		else
			ErrorLog "${SAMPLEID}" "$ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.bam already exists and force was NOT used!";
		fi
	else
		ErrorLog "At the moment it is only supported for CALL_TYPE=h.sapiens for DESIGN_TYPE=PCR!";
	fi
else 
	ErrorLog "DESIGN_TYPE is unknown (only PCR is supported)!";
fi

# Check if SamBamFilesBuild worked
if [ "$?" != "0" ]; then
	ErrorLog "${SAMPLEID}" "Failed in exporting genome alignment to BAM files.";
else

	SuccessLog "${SAMPLEID}" "Passed exporting genome alignment to BAM files.";
fi
