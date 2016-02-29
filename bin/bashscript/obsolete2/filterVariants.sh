#!/bin/bash
#
# Script for running ensembl variant effect prediction from SNPMania output
#
##SBATCH -p core -n 1 -t 05:00:00
#SBATCH -p core -t 00:15:00 --qos=short

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/filteredVariants" ]; then
	mkdir $ROOT_PATH/filteredVariants;
fi

if [[ ! -s "$ROOT_PATH/filteredVariants/${SAMPLEID}.${CALL_TYPE}.variants" || ! -z $FORCE ]]; then

	# check if normal sample is provided and include it in the flags
	if [[ -e $NORMAL_SAMPLEID || $NORMAL_SAMPLEID != "false" ]]; then
		echo "SNVFILTER_FLAGS=\"$SNVFILTER_FLAGS -n $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.${CALL_TYPE}\";";	
		SNVFILTER_FLAGS="$SNVFILTER_FLAGS -n $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.${CALL_TYPE}";
	fi

	# include the previous output file if running this sample again
	if [ -e $ROOT_PATH/filteredVariants/${SAMPLEID}.${CALL_TYPE}.variants ]; then
		SNVFILTER_FLAGS="$SNVFILTER_FLAGS -p $ROOT_PATH/filteredVariants/${SAMPLEID}.${CALL_TYPE}.variants";
	fi
	if [[ -e ${GENEFILE} && ${GENEFILE}!="false" ]]; then 
		SNVFILTER_FLAGS="$SNVFILTER_FLAGS -g ${GENEFILE}";
	fi
	# handy when extracting e.g. 10% of median, change -ds and -di to $MEDIAN, however, remember that the NORMAL sample will also be filtered at this depth
#	MEDIAN=`$HOME/testHash/calcMedian.pl -i $ROOT_PATH/SNPmania/${SAMPLEID}.${CALL_TYPE}.variations -r 0.1 -d 20 -f /bubo/home/h20/elinfalk/tmp/LarryRebeqa.withUceChanged.roi`;

	# extract variants
#	echo "zcat $ONEKGENOMEPROJECT_AF | $SERA_PATH/bin/perlscript/filterSNPs.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.${CALL_TYPE}.variations -i $ROOT_PATH/SNPmania/${SAMPLEID}.${CALL_TYPE}.insertions -d $ROOT_PATH/SNPmania/${SAMPLEID}.${CALL_TYPE}.deletions -r $ROIFILE -o $ROOT_PATH/filteredVariants/${SAMPLEID}.${CALL_TYPE}.variants -a /dev/stdin -f $ROOT_PATH/filteredVariants/${SAMPLEID}.${CALL_TYPE}.variants.log $SNVFILTER_FLAGS;";
	
	zcat $ONEKGENOMEPROJECT_AF | $SERA_PATH/bin/perlscript/filterSNPs.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.${CALL_TYPE}.variations -i $ROOT_PATH/SNPmania/${SAMPLEID}.${CALL_TYPE}.insertions -d $ROOT_PATH/SNPmania/${SAMPLEID}.${CALL_TYPE}.deletions -r $ROIFILE -o $ROOT_PATH/filteredVariants/${SAMPLEID}.${CALL_TYPE}.variants -a /dev/stdin -f $ROOT_PATH/filteredVariants/${SAMPLEID}.${CALL_TYPE}.variants.log $SNVFILTER_FLAGS;
		

	if [ "$?" != "0" ]; then
		ErrorLog $SAMPLEID "Failed in filtering variants";
	else
		SuccessLog $SAMPLEID "Passed filtering variants";
	fi

else
	SuccessLog $SAMPLEID "Filtered files already exists, skipping step";
fi

