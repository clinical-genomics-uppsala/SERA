#!/bin/bash
#
# Script to export aligment archives as SAM-files.
#
#SBATCH -p core -n 2
#SBATCH -t 02:00:00

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Gathering selector properties";

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/junctionHits" ]; then
	mkdir $ROOT_PATH/junctionHits;
fi

# get parameters
if [[ ! -e $ROOT_PATH/junctionHits/${SAMPLEID}.junctionHits.length.gc || ! -z $FORCE ]]; then
	zcat $ROOT_PATH/MosaikRef/${REFSEQ}.selection.seq.gz | awk '{print $0"\t"length($5)}' - | sort | perl $SERA_PATH/bin/perlscript/doGCanalysis.pl -i /dev/stdin -o $ROOT_PATH/refFiles/${REFSEQ}.selection.seq.length.gc -c 5 -s 2 -n 3 -p 25,25;
	paste $ROOT_PATH/junctionHits/${SAMPLEID}.junctionHits $ROOT_PATH/refFiles/${REFSEQ}.selection.seq.length.gc | awk '{print $1"\t"$2"\t"$8"\t"$9"\t"$10}' > $ROOT_PATH/junctionHits/${SAMPLEID}.junctionHits.length.gc;
fi


# Check if worked
if [ "$?" != "0" ]; then
	ErrorLog ${SAMPLEID} "Failed in gathering selector properties";
else
	SuccessLog ${SAMPLEID} "Passed gathering selector properties";
fi
