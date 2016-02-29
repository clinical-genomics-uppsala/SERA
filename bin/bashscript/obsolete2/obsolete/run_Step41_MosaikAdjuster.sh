#!/bin/bash
#
# Script to run MosaikAdjuster
#
#SBATCH -p node -n 1
#SBATCH -t 01:00:00

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog "$SAMPLEID" "Starts MosaikAdjuster ...";

# Run MosaikAdjuster for SOLiD
if [ ${PLATFORM} = "SOLiD" ]; then
	if [ -e $ROOT_PATH/MosaikText/${SAMPLEID}.${CALL_TYPE}.bam.gz ]; then
		zcat $ROOT_PATH/SamBamFiles/${SAMPLEID}.${CALL_TYPE}.corrSel.bam.gz | $ROOT_PATH_JSNPMANIA/MosaikAdjuster.sh -i /dev/stdin -o $ROOT_PATH/SamBamFiles/${SAMPLEID}.${CALL_TYPE}.corrSel.MosaikAdjusted.bam -r $ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq -n -l;
		gzip -f $ROOT_PATH/SamBamFiles/${SAMPLEID}.${CALL_TYPE}.corrSel.MosaikAdjusted.bam;
	else
		ErrorLog "${SAMPLEID}" "doesn't have a fixed genome BAM-file.";
	fi
	
# Illumina doesn't need MosaikAdjuster
elif [ ${PLATFORM} = "Illumina" ]; then
	SuccessLog "${SAMPLEID}" "Skipping MosaikAdjuster for Illumina data";
else
	ErrorLog "${SAMPLEID} platform ${PLATFORM} not recognized (should be SOLiD or Illumina)";
	exit 1;
fi

# Check if it worked
if [ "$?" != "0" ]; then
	ErrorLog "${SAMPLEID}" "Failed in MosaikAdjuster";
else
	SuccessLog "${SAMPLEID}" "Passed MosaikAdjuster";
fi
