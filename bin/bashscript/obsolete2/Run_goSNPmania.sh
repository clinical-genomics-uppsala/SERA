#!/bin/bash
#
# Script to run goSNPmania
#
#SBATCH -p core -n 1
#SBATCH -t 30:00
##SBATCH --qos=short

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts goSNPmania with flags ${SNPMANIAFLAGS}";

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/goSNPmania" ]; then
	mkdir $ROOT_PATH/goSNPmania;
fi

if [ ${HAPMAP_SNP_REF} != "false" ]; then
	if [ ${READS} == "true" ]; then
		if [ ${CALL_TYPE} == "h.sapiens" ]; then
			if [[ ! -e $ROOT_PATH/goSNPmania/${SAMPLEID}.h.sapiens.variations.goSNPmania.gz || ! -z $FORCE ]]; then
				if [ -e $ROOT_PATH/SNPmania/${SAMPLEID}.h.sapiens.variations ]; then
					perl $SERA_PATH/bin/perlscript/goSNPmania2.pl -i $ROOT_PATH/SNPmania/${SAMPLEID}.h.sapiens.variations -snp_ref ${HAPMAP_SNP_REF},3,4,8 -o $ROOT_PATH/goSNPmania/${SAMPLEID}.h.sapiens.variations.goSNPmania -nc2chr;
					gzip -f $ROOT_PATH/goSNPmania/${SAMPLEID}.h.sapiens.variations.goSNPmania;
				fi
		fi
		elif [ ${CALL_TYPE} == "ampregion" ]; then
			if [[ ! -e $ROOT_PATH/goSNPmania/${SAMPLEID}.ampregion.variations.goSNPmania.gz || ! -z $FORCE ]]; then
				if [ -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampregion.variations ]; then
					perl $SERA_PATH/bin/perlscript/goSNPmania2.pl -i $ROOT_PATH/SNPmania/${SAMPLEID}.ampregion.variations -snp_ref ${HAPMAP_SNP_REF},3,4,8 -o $ROOT_PATH/goSNPmania/${SAMPLEID}.ampregion.variations.goSNPmania -nc2chr;
					gzip -f $ROOT_PATH/goSNPmania/${SAMPLEID}.ampregion.variations.goSNPmania;
				fi
			fi
		fi
	fi
else
	ErrorLog ${SAMPLEID} " - a hapmap file was not stated, skipping step!";
fi

# Check if MosaikTextBuild worked
if [ "$?" != "0" ]; then
	ErrorLog ${SAMPLEID} "Failed in goSNPmania";
else
	SuccessLog ${SAMPLEID} "Passed goSNPmania";
fi

