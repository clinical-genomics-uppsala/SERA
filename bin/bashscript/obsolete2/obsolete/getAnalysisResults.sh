#!/bin/bash
#SBATCH -p core
#SBATCH -t 00:15:00
#SBATCH --qos=short
#
# Script for running MosaikAligner against a genome and save unaligned reads.
#

# Include functions
#. $SERA_PATH/includes/logging.sh;

SuccessLog ${REFSEQ} "Starts collecting data...";

# Check if the directory exists, 
if [ ! -d "$ROOT_PATH/TextResults" ]; then
	ErrorLog ${REFSEQ} "Unable to check diretory $ROOT_PATH/TextResults";
	exit 1;
fi

if [[ ! -e $ROOT_PATH/TextResults/Sample_results.txt || ! -z $FORCE ]]; then

	types=( "SPEC_REGION_UNIQ" "SPEC_ROI_UNIQ" "SPEC_REGION_ALL" "SPEC_ROI_ALL" "DEPTH_SEQREGION_UNIQ" "DEPTH_SEQROI_UNIQ" "DEPTH_UNIQ" "DEPTH_SEQREGION_ALL" "DEPTH_SEQROI_ALL" "DEPTH_ALL" );

	# create file headers
	head="#Sample";
	for type in ${types[@]}; do
		head="$head\t$type";
	done;
	echo -e $head > $ROOT_PATH/TextResults/Sample_results.txt;



	# loop through all files
	samples=( `grep -lir "SPEC*" --exclude="*_*" $ROOT_PATH/TextResults | xargs` );

	for sample in ${samples[@]}; do

		#get sampleid
		conc="`basename $sample | awk '{split($0,a,"."); print a[1]}'`";

		# loop through all combinations
		for type in ${types[@]}; do

			# check each row and print type
			conc="$conc\t`cat $sample | grep $type | tail -n 1 | awk '{split($0,a,"="); print a[2]}'`";

		done

		# Write results.
		echo -e "$conc" >> $ROOT_PATH/TextResults/Sample_results.txt;

	done;

fi
