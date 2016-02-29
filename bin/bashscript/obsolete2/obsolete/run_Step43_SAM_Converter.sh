#!/bin/bash
#
# Script to export aligment archives as SAM-files.
#
#SBATCH -p node -n 1
#SBATCH -t 01:00:00

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts converting SAM-files ...";

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/SamBamFiles" ]; then
	mkdir $ROOT_PATH/SamBamFiles;
fi

ADJUSTED="";
if [ $PLATFORM = "SOLiD" ]; then
	ADJUSTED=".MosaikAdjusted";
fi

if [[ ! -e $ROOT_PATH/SamBamFiles/${SAMPLEID}.${CALL_TYPE}.corrSel${ADJUSTED}.corrected.bam || ! -z $FORCE ]]; then
	$ROOT_PATH_JSNPMANIA/SAMConverter.sh -i $ROOT_PATH/SamBamFiles/${SAMPLEID}.${CALL_TYPE}.corrSel${ADJUSTED}.bam -o $ROOT_PATH/SamBamFiles/${SAMPLEID}.${CALL_TYPE}.corrSel${ADJUSTED}.corrected.bam -s;
	gzip -f $ROOT_PATH/SamBamFiles/${SAMPLEID}.${CALL_TYPE}.corrSel${ADJUSTED}.corrected.bam;

else
	ErrorLog "${SAMPLEID} platform ${PLATFORM} not recognized (should be SOLiD or Illumina)";
	exit 1;
fi

# Check if it worked
if [ "$?" != "0" ]; then
	ErrorLog ${SAMPLEID} "Failed in SAMConverter";
else
	SuccessLog ${SAMPLEID} "Passed SAMConverter";
fi
