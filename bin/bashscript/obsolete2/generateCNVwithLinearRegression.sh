#!/bin/bash
#
#
#SBATCH -p core -n 1
#SBATCH -t 02:00:00 
##SBATCH -p core -t 00:15:00 --qos=short

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/cnvAnalysisLinearRegression" ]; then
	mkdir $ROOT_PATH/cnvAnalysisLinearRegression;
fi

if [[ ! -z $NORMAL_SAMPLEID  && $NORMAL_SAMPLEID != "false" ]]; then

	SuccessLog $SAMPLEID "Creating CNV file for $SAMPLEID with sample $NORMAL_SAMPLEID as normal pair";

	if [[ ! -e $ROOT_PATH/cnvAnalysisLinearRegression/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.cnvLinearRegression.pdf || ! -z $FORCE ]]; then
		awk '{if($1>0){print $0}}' $ROOT_PATH/SNPmania/${SAMPLEID}.${CALL_TYPE}.variations > $ROOT_PATH/SNPmania/${SAMPLEID}.${CALL_TYPE}.fixed.variations;
		awk '{if($1>0){print $0}}' $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.${CALL_TYPE}.variations > $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.${CALL_TYPE}.fixed.variations;
		$SERA_PATH/bin/Rscript/cnvPlot_with_linearRegression.R $ROOT_PATH/SNPmania/${SAMPLEID}.${CALL_TYPE}.fixed.variations $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.${CALL_TYPE}.fixed.variations $GENEFILE $CNV_REGRESSION_FLAGS $ROOT_PATH/cnvAnalysisLinearRegression/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.cnvLinearRegression.pdf $ROOT_PATH/cnvAnalysisLinearRegression/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.cnvLinearRegression.txt "${SAMPLEID} ${NORMAL_SAMPLEID}";
		
		rm $ROOT_PATH/SNPmania/${SAMPLEID}.${CALL_TYPE}.fixed.variations;
		rm $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.${CALL_TYPE}.fixed.variations;

		# check if passing
		if [ "$?" != "0" ]; then
			ErrorLog $SAMPLEID "Could not create CNV file for tumor/normal pair ${SAMPLEID}/${NORMAL_SAMPLEID}";
		else
			SuccessLog $SAMPLEID "CNV files created with sample $NORMAL_SAMPLEID as normal pair";
		fi

	else
		ErrorLog $SAMPLEID "Previous output exists, skipping step";
	fi

else
	SuccessLog $SAMPLEID "No tumor/normal pair specified, skipping step";
fi
