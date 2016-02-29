#!/bin/bash
#
# Script runs base mapping on ampregion.jointSplit.
#
#SBATCH -p core -n 1
#SBATCH -t 02:00:00 
##SBATCH -p core -t 00:15:00 --qos=short

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/cnvAnalysis" ]; then
	mkdir $ROOT_PATH/cnvAnalysis;
fi

if [[ ! -z $NORMAL_SAMPLEID  || $NORMAL_SAMPLEID != "false" ]]; then

	SuccessLog $SAMPLEID "Creating CNV file for $SAMPLEID with sample $NORMAL_SAMPLEID as normal pair";

	if [[ ! -e $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.normalMinCov$MINCOV.cnv_perBase.median.pdf || ! -z $FORCE ]]; then

		MINCOV=20;
		MINALLELEDIFF=0.2;

		# create
		paste $ROOT_PATH/SNPmania/${SAMPLEID}.${CALL_TYPE}.variations $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.${CALL_TYPE}.variations | $SERA_PATH/bin/awkscript/filterNormalCov_var.awk -v min=$MINCOV - | awk '{print $4"\t"$5"\t"$1"\t"$14"\t"$15"\t"$11}' - > $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.variations.normalMinCov$MINCOV;
		paste $ROOT_PATH/SNPmania/${SAMPLEID}.${CALL_TYPE}.variations $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.${CALL_TYPE}.variations | $SERA_PATH/bin/awkscript/filterNormalCov_var.awk -v min=0 - | awk '{print $4"\t"$5"\t"$1"\t"$14"\t"$15"\t"$11}' - > $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.variations.normalMinCov0;

		# create tumor and normal median cov 20 
		perl $SERA_PATH/bin/perlscript/CalculateCNVmean.pl -h 3 -i $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.variations.normalMinCov$MINCOV -o $ROOT_PATH/cnvAnalysis/${SAMPLEID}.${CALL_TYPE}.median.normalMinCov$MINCOV;
		perl $SERA_PATH/bin/perlscript/CalculateCNVmean.pl -h 6 -i $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.variations.normalMinCov$MINCOV -o $ROOT_PATH/cnvAnalysis/${NORMAL_SAMPLEID}.${CALL_TYPE}.median.normalMinCov$MINCOV;

		# perBase median
		sort -k1,2h $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.variations.normalMinCov0 | perl $SERA_PATH/bin/perlscript/makeCNVfile_perBase.pl -i /dev/stdin -mn $ROOT_PATH/cnvAnalysis/${NORMAL_SAMPLEID}.${CALL_TYPE}.median.normalMinCov$MINCOV -mt $ROOT_PATH/cnvAnalysis/${SAMPLEID}.${CALL_TYPE}.median.normalMinCov$MINCOV -o $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.normalMinCov$MINCOV.cnv_perBase.median -pdf $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.normalMinCov$MINCOV.cnv_perBase.median.pdf -t $SAMPLEID -n $NORMAL_SAMPLEID -min $MINCOV -a median;
#		sort -k1,2h $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.variations.normalMinCov0 | perl $SERA_PATH/bin/perlscript/makeCNVfile_perBase.pl -i /dev/stdin -mn $ROOT_PATH/cnvAnalysis/${NORMAL_SAMPLEID}.${CALL_TYPE}.median.normalMinCov0 -mt $ROOT_PATH/cnvAnalysis/${SAMPLEID}.${CALL_TYPE}.median.normalMinCov0 -o $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.normalMinCov0.cnv_perBase.median -pdf $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.normalMinCov0.cnv_perBase.median.pdf -t $SAMPLEID -n $NORMAL_SAMPLEID -min 0 -a median;

		# perRegion median
		sort -k1,2h $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.variations.normalMinCov0 | perl $SERA_PATH/bin/perlscript/makeCNVfile_perRegion.pl -i /dev/stdin -mn $ROOT_PATH/cnvAnalysis/${NORMAL_SAMPLEID}.${CALL_TYPE}.median.normalMinCov$MINCOV -mt $ROOT_PATH/cnvAnalysis/${SAMPLEID}.${CALL_TYPE}.median.normalMinCov$MINCOV -o $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.normalMinCov$MINCOV.cnv_perRegion.median -pdf $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.normalMinCov$MINCOV.cnv_perRegion.median.pdf -t $SAMPLEID -n $NORMAL_SAMPLEID -min $MINCOV -a median;
#		sort -k1,2h $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.variations.normalMinCov0 | perl $SERA_PATH/bin/perlscript/makeCNVfile_perRegion.pl -i /dev/stdin -mn $ROOT_PATH/cnvAnalysis/${NORMAL_SAMPLEID}.${CALL_TYPE}.median.normalMinCov0 -mt $ROOT_PATH/cnvAnalysis/${SAMPLEID}.${CALL_TYPE}.median.normalMinCov0 -o $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.normalMinCov0.cnv_perRegion.median -pdf $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.normalMinCov0.cnv_perRegion.median.pdf -t $SAMPLEID -n $NORMAL_SAMPLEID -min 0 -a median;

		# compress files
#		gzip -f $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.variations.normalMinCov0;
#		gzip -f $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.variations.normalMinCov$MINCOV;
#		gzip -f $ROOT_PATH/cnvAnalysis/${SAMPLEID}_${NORMAL_SAMPLEID}.${CALL_TYPE}.cnv;
	
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
