#!/bin/bash
#
# Script for running ensembl variant effect prediction from SNPMania output
#
##SBATCH -p node -n 1
##SBATCH -t 05:00:00
#SBATCH --qos=short -t 00:15:00 -p core

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/TextResults" ]; then
	mkdir $ROOT_PATH/TextResults;
fi

# function to sort maintaining the comment at first position
function sortFile {
	head -n 1 ${1} > ${2};
	let lines=`wc -l ${1} | awk '{print $1}'`-1;
	tail -n $lines ${1} | sort ${3} ${4} ${5} >> ${2};
}


# function to grep the most severe consequences predicted
# also greps ampregion/roi if second argument set as ampregion or roi
function grepSevere {

	# define output name, strip extension
	output="`echo "${1}" | cut -d'.' -f1`.severe.txt";

	# keep header
	head -n 1 ${1} > $output;
	let lines=`wc -l ${1} | awk '{print $1}'`-1;

	# extract consequences
#	cons=("stop gained" "stop lost" "frameshift coding" "non-synonymous coding" "splice site");
	cons=("STOP" "FRAMESHIFT" "NON-SYNONYMOUS" "SPLICE");
	for c in ${cons[@]}; do
		tail -n $lines ${1} | grep $c | grep "${2}" >> $output.unsorted;
	done

	cat $output.unsorted | sort | uniq >> $output;
	rm $output.unsorted;
}	


if [[ ! -s "$ROOT_PATH/TextResults/variants_positions.txt" || ! -z $FORCE ]]; then

	files=`ls $ROOT_PATH/effectPrediction/ | grep ".predicted" | wc -l`;

	if [[ $files > 0 ]]; then

#		perl $SERA_PATH/bin/perlscript/gatherMutatedRegions.pl -r $ROOT_PATH/refFiles/${REFSEQ}.ampregion -s $ROIFILE -i $ROOT_PATH/effectPrediction -o $ROOT_PATH/TextResults -h /bubo/home/h20/elinfalk/CLL_analysis/names.txt -g /bubo/home/h20/elinfalk/CLL_analysis/samples.txt -e roi,severe;
		perl $SERA_PATH/bin/perlscript/gatherMutatedRegions.pl -r $ROOT_PATH/refFiles/${REFSEQ}.ampregion -s $ROIFILE -i $ROOT_PATH/effectPrediction -o $ROOT_PATH/TextResults -g /bubo/home/h20/elinfalk/CLL_analysis/samples.txt -h /bubo/home/h20/elinfalk/CLL_analysis/names.txt -e roi,severe,variations

		if [ "$?" != "0" ]; then
			ErrorLog $SAMPLEID "Failed creating data tables for VEP";
		else
			SuccessLog $SAMPLEID "Variant calls summarized for $files samples and reported in table";
		fi

		# extract all files with severe variations
#		grepSevere $ROOT_PATH/TextResults/variants_positions.txt ROI;
#		grepSevere $ROOT_PATH/TextResults/variants_samples.txt ROI;

		# sort all files, keeping line header
#		sortFile $ROOT_PATH/TextResults/variants_positions.txt $ROOT_PATH/TextResults/variants_positions.sorted.txt -k2,2n -k3,3n;
#		sortFile $ROOT_PATH/TextResults/variants_samples.txt $ROOT_PATH/TextResults/variants_samples.sorted.txt -k1,1 -k3,3n -k4,4n;
#		sortFile $ROOT_PATH/TextResults/variants_genes.txt $ROOT_PATH/TextResults/variants_genes.sorted.txt -k6,6nr;

#		sortFile $ROOT_PATH/TextResults/variants_positions.severe.txt $ROOT_PATH/TextResults/variants_positions.severe.sorted.txt -k2,2n -k3,3n;
#		sortFile $ROOT_PATH/TextResults/variants_samples.severe.txt $ROOT_PATH/TextResults/variants_samples.severe.sorted.txt -k1,1 -k3,3n -k4,4n;

		# extract deletions/variations/insertions
#		awk '{if($7~/D/) print $0}' $ROOT_PATH/TextResults/variants_samples.severe.sorted.txt > $ROOT_PATH/TextResults/variants_samples.severe.sorted.deletions.txt
#		awk '{if($7~/I/) print $0}' $ROOT_PATH/TextResults/variants_samples.severe.sorted.txt > $ROOT_PATH/TextResults/variants_samples.severe.sorted.insertions.txt
#		awk '{if($7!~/D/ && $7!~/I/) print $0}' $ROOT_PATH/TextResults/variants_samples.severe.sorted.txt > $ROOT_PATH/TextResults/variants_samples.severe.sorted.variations.txt

	else
		ErrorLog "No variant effect prediction files found, skipping this step.";
	fi

else
	ErrorLog $SAMPLEID "Tables for variant effect prediction already exists.";
fi
