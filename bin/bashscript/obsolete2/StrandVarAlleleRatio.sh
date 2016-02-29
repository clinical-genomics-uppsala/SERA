#!/bin/bash
#
#
#SBATCH -p core -n 1
#SBATCH -t 15:00
##SBATCH --qos=short

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/strandVariantAlleleRatio" ]; then
        mkdir $ROOT_PATH/strandVariantAlleleRatio;
fi

SuccessLog "${SAMPLEID}" "Calculating variant allele ratio per strand based on strand info...";


# Check that Reads are true
if [ ${READS} == "true" ]; then
	# Check that the CALL_TYPE is h.sapiens
    if [ ${CALL_TYPE} == "h.sapiens" ]; then
		# Check that the output file doesn't exist
		if [[ ! -e $ROOT_PATH/strandVariantAlleleRatio/${SAMPLEID}.$CALL_TYPE.strandVariantAlleleRatio || ! -z $FORCE ]]; then
			if [ -e $ROOT_PATH/SNPmania/${SAMPLEID}.h.sapiens.ampliconmapped.variations ]; then
				perl $SERA_PATH/bin/perlscript/StrandVarAlleleRatio_strandInfo.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.${CALL_TYPE}.ampliconmapped.variations -o $ROOT_PATH/strandVariantAlleleRatio/${SAMPLEID}.$CALL_TYPE.strandVariantAlleleRatio -minRD 10 -vMinRD 2;
			elif [ -e $ROOT_PATH/SNPmania/${SAMPLEID}.h.sapiens.variations ]; then
				perl $SERA_PATH/bin/perlscript/StrandVarAlleleRatio_strandInfo.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.${CALL_TYPE}.variations -o $ROOT_PATH/strandVariantAlleleRatio/${SAMPLEID}.$CALL_TYPE.strandVariantAlleleRatio -minRD 10 -vMinRD 2;
			else
				ErrorLog "${SAMPLEID}" "None of the possible input files exist ($ROOT_PATH/SNPmania/${SAMPLEID}.h.sapiens.ampliconmapped.variations or $ROOT_PATH/SNPmania/${SAMPLEID}.h.sapiens.variations)!";
			fi
		else
			ErrorLog "${SAMPLEID}" "$ROOT_PATH/strandVariantAlleleRatio/${SAMPLEID}.$CALL_TYPE.strandVariantAlleleRatio already exists and -f was not used!"
		fi
	else
		ErrorLog "${SAMPLEID}" "The analysis is only supported for CALL_TYPE h.sapiens so far!";
	fi
else
	ErrorLog "${SAMPLEID}" "READS has to be true to run the analysis!";
fi

if [ "$?" != "0" ]; then
        ErrorLog "${SAMPLEID}" "Failed in calculating variant allele ratio per strand based on strand info!";
else
        SuccessLog "${SAMPLEID}" "Passed calculating variant allele ratio per strand based on strand info!";
fi
