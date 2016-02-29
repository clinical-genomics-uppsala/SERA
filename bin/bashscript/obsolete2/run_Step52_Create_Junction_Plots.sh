#!/bin/bash
#
# Create plot based on junction hits and selector properties
#
#SBATCH -p core -n 1
#SBATCH -t 01:00:00

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Creating junction plots";

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/junctionPdfs" ]; then
	mkdir $ROOT_PATH/junctionPdfs;
fi

# length bin plot
if [[ ! -e $ROOT_PATH/junctionPdfs/${SAMPLEID}.junctionHits.length.pdf || ! -z $FORCE ]]; then
	cat $ROOT_PATH/junctionHits/${SAMPLEID}.junctionHits.length.gc | perl $SERA_PATH/bin/perlscript/lengthBins.pl -i /dev/stdin -o $ROOT_PATH/junctionHits/${SAMPLEID}.junctionHits.lengthBin -c 3 -h 2 -pdf $ROOT_PATH/junctionPdfs/${SAMPLEID}.junctionHits.length.pdf -t ${SAMPLEID} -type junction;
	if [ "$?" != "0" ]; then
		ErrorLog ${SAMPLEID} "Failed in creating length bin plot";
	else
		SuccessLog ${SAMPLEID} "Passed creating length bin plot";
	fi
fi

# gc bin plot
if [[ ! -e $ROOT_PATH/junctionPdfs/${SAMPLEID}.junctionHits.gc.pdf || ! -z $FORCE ]]; then
	cat $ROOT_PATH/junctionHits/${SAMPLEID}.junctionHits.length.gc | perl $SERA_PATH/bin/perlscript/gcBins.pl -i /dev/stdin -o $ROOT_PATH/junctionHits/${SAMPLEID}.junctionHits.gcBin -c 4 -h 2 -pdf $ROOT_PATH/junctionPdfs/${SAMPLEID}.junctionHits.gc.pdf -t ${SAMPLEID} -type junction -xlabel "Mean fragment GC bins";
	if [ "$?" != "0" ]; then
		ErrorLog ${SAMPLEID} "Failed in creating gc bin plot";
	else
		SuccessLog ${SAMPLEID} "Passed creating gc bin plot";
	fi
fi

# reaction bin plot
if [[ ! -e $ROOT_PATH/junctionPdfs/${SAMPLEID}.junctionHits.reactions.pdf || ! -z $FORCE ]]; then
	cat $ROOT_PATH/junctionHits/${SAMPLEID}.junctionHits.length.gc | perl $SERA_PATH/bin/perlscript/reactionBins.pl -i /dev/stdin -o $ROOT_PATH/junctionHits/${SAMPLEID}.junctionHits.reactionBin -h 2 -pdf $ROOT_PATH/junctionPdfs/${SAMPLEID}.junctionHits.reactions.pdf -t ${SAMPLEID} -type junction -r $SELECTIONFILE -p;
	if [ "$?" != "0" ]; then
		ErrorLog ${SAMPLEID} "Failed in creating reaction bin plot";
	else
		SuccessLog ${SAMPLEID} "Passed creating reaction bin plot";
	fi
fi
