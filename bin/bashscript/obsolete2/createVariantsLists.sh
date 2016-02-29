#!/bin/bash
#
# Script for running ensembl variant effect prediction from SNPMania output
#
##SBATCH -p core -n 1
##SBATCH -t 05:00:00
#SBATCH --qos=short

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/TextResults" ]; then
	mkdir $ROOT_PATH/TextResults;
fi

if [[ ! -s "$ROOT_PATH/TextResults/samples.txt" || ! -z $FORCE ]]; then

	files=`ls $ROOT_PATH/filteredVariants/ | grep "variants$" | wc -l`;

	if [[ $files > 0 ]]; then

		perl $SERA_PATH/bin/perlscript/gatherMutationsSimple.pl -i $ROOT_PATH/filteredVariants -o $ROOT_PATH/TextResults -g $SAMPLES_INFORMATION;

		if [ "$?" != "0" ]; then
			ErrorLog $SAMPLEID "Failed creating data tables for variants";
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
		ErrorLog "No variation files found in specified folder";
	fi

else
	ErrorLog $SAMPLEID "Tables for variant effect prediction already exists.";
fi
