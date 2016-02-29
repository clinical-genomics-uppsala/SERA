#!/bin/bash
#
# Creates junction reference files for Mosaik
#
#SBATCH -p core -n 1
#SBATCH -t 01:00:00

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts extracting all junctions from alignment";

# check whether this is a MDA design
if [[ "$SEQUNCING_TAG" == "false" ]]&&[[ "$DESIGN_TYPE" == "MDA" ]]; then
	
	# Check if the directory exists, if not create it
	if [ ! -d "$ROOT_PATH/junctionHits" ]; then
		mkdir $ROOT_PATH/junctionHits;
	fi

	# extracts all junctions from BAM file
	if [[ ! -e $ROOT_PATH/junctionHits/${SAMPLEID}.junctionHits || ! -z $FORCE ]]; then
		$ROOT_PATH/pdfs/${SAMPLEID}.${CALL_TYPE}.hapmap.pdf
		samtools -h $ROOT_PATH/SamBamFiles/$SAMPLEID.ampregion.uniq.bam -o /dev/stdout | awk '{if(substr($1,0,1)!="@"){print $3}}' | sort | uniq -c | awk '{print $2"\t"$1}' | awk '{FS="[#\t]" } { print $1"\t"$6 }' | sort - | perl $SERA_PATH/bin/perlscript/listAllJunctions.pl -j /dev/stdin -o $ROOT_PATH/junctionHits/${SAMPLEID}.junctionHits -s ${SELECTIONFILE};
	fi

elif [ "$DESIGN_TYPE" != "MDA" ]; then
	SuccessLog ${SAMPLEID} "No junctions included for non MDA designs.";
fi

if [ "$?" != "0" ]; then
	ErrorLog "${SAMPLEID}" "Failed in counting hits per junction";
else
	SuccessLog "${SAMPLEID}" "Passed counting hits per junction";
fi
