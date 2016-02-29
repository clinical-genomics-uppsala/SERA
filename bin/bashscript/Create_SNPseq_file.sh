#!/bin/bash
#
# Script creates SNPseq file.
#
#SBATCH -p devcore  -n 2
#SBATCH -t 01:00:00

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog "${SAMPLEID}" "Starts creating ampregion SNPseq files ...";

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/refFiles" ]; then
	mkdir $ROOT_PATH/refFiles;
fi

# Check if the SNPmania reference file exists or if force is used
if [[ ! -e $ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq || ! -z $FORCE ]]; then
	# Check that the ampregion reference file exists
	if [ -e $ROOT_PATH/refFiles/${REFSEQ}.ampregion ]; then
		perl $DOWNLOAD2FASTA -th 2 -s $ROOT_PATH/refFiles/${REFSEQ}.ampregion -o $ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq -d $BLAST_DB -t full;
	else
		ErrorLog "$SAMPLEID" "$ROOT_PATH/refFiles/${REFSEQ}.ampregion does NOT exist, run step 0 to create this file first!";
	fi
else
	ErrorLog "$SAMPLEID" "$ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq already exists and force was NOT used!";
fi

# Check if creating ampregion file worked
if [ "$?" != "0" ]; then
	ErrorLog "${REFSEQ}" "Failed in creating ampregion SNPseq file";
else
	SuccessLog "${REFSEQ}" "Passed in creating ampregion SNPseq file";
fi
