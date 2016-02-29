#!/bin/bash
#
# Script runs base mapping on ampregion.
#
#SBATCH -p core -n 1
#SBATCH -t 15:00
##SBATCH --qos=short

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/log2cov" ]; then
	mkdir $ROOT_PATH/log2cov;
fi

SuccessLog "${SAMPLEID}" "Running sequenced region base mapping...";

# BaseMapping against ampregion all

# BaseMapping against seqregion for reads
if [ ${READS} == "true" ]; then
	if [ ${CALL_TYPE} == "h.sapiens" ]; then
		if [[ ! -e $ROOT_PATH/log2cov/${SAMPLEID}.h.sapiens.seqroi.log2cov.pdf || ! -z $FORCE ]]; then
			if [ -e $ROOT_PATH/refFiles/${REFSEQ}.seqroi ]; then
				if [ -e $ROOT_PATH/SNPmania/${SAMPLEID}.h.sapiens.variations ]; then
					perl $SERA_PATH/bin/perlscript/ExtractSNPsFromSNPmania.pl -s $ROOT_PATH/SNPmania/${SAMPLEID}.h.sapiens.variations -r $ROOT_PATH/refFiles/${REFSEQ}.seqroi -o $ROOT_PATH/log2cov/${SAMPLEID}.h.sapiens.seqroi.variations
					$SERA_PATH/bin/Rscript/SNPmania_to_log2cov.R $ROOT_PATH/log2cov/${SAMPLEID}.h.sapiens.seqroi.variations $ROOT_PATH/log2cov/${SAMPLEID}.h.sapiens.seqroi.log2cov.pdf "${SAMPLEID} reads" 20 1
				elif [ -e $ROOT_PATH/SNPmania/${SAMPLEID}.h.sapiens.ampliconmapped.variations ]; then
					perl $SERA_PATH/bin/perlscript/ExtractSNPsFromSNPmania.pl -s $ROOT_PATH/SNPmania/${SAMPLEID}.h.sapiens.ampliconmapped.variations -r $ROOT_PATH/refFiles/${REFSEQ}.seqroi -o $ROOT_PATH/log2cov/${SAMPLEID}.h.sapiens.seqroi.variations
					$SERA_PATH/bin/Rscript/SNPmania_to_log2cov.R $ROOT_PATH/log2cov/${SAMPLEID}.h.sapiens.seqroi.variations $ROOT_PATH/log2cov/${SAMPLEID}.h.sapiens.seqroi.log2cov.pdf "${SAMPLEID} reads" 20 1
				else
					ErrorLog "${SAMPLEID}" "None of the accepted SNPmania variation files existed ($ROOT_PATH/SNPmania/${SAMPLEID}.h.sapiens.variations or $ROOT_PATH/SNPmania/${SAMPLEID}.h.sapiens.ampliconmapped.variations)!";
				fi
			else
				ErrorLog "${SAMPLEID}" "$ROOT_PATH/refFiles/${REFSEQ}.seqroi or does NOT exist!";
			fi
		else
			ErrorLog "${SAMPLEID}" "$ROOT_PATH/log2cov/${SAMPLEID}.h.sapiens.seqroi.log2cov.pdf already exists and -f was not used!";
		fi
	else
		ErrorLog "${SAMPLEID}" "Only supported for CALL_TYPE h.sapiens so far!";
	fi
else
	ErrorLog "${SAMPLEID}" "READS has to be true to run the analysis!";
fi

if [ "$?" != "0" ]; then
	ErrorLog "${SAMPLEID}" "Failed in log2coverage plot.";
else
	SuccessLog "${SAMPLEID}" "Failed in log2coverage plot.";
fi
