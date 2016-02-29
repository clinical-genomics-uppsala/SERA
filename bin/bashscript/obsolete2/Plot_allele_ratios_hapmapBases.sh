#!/bin/bash
#
# Script to plot allele ratios
#
#SBATCH -p core -n 1
#SBATCH -t 01:00:00
##SBATCH --qos=short

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts plotAlleleRatio";

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/plotAlleleRatio" ]; then
	mkdir $ROOT_PATH/plotAlleleRatio;
fi
if [ ! -d "$ROOT_PATH/plotPdfs" ]; then
	mkdir $ROOT_PATH/plotPdfs;
fi

if [ ${READS} == "true" ]; then
	if [ ${CALL_TYPE} == "h.sapiens" ]; then
		if [[ -e $ROOT_PATH/goSNPmania/${SAMPLEID}.h.sapiens.variations.goSNPmania.gz ]]; then
			if [[ ! -e $ROOT_PATH/plotPdfs/${SAMPLEID}.h.sapiens.hapmap.pdf || ! -z $FORCE ]]; then
				# extract hapmap snps 
				zcat $ROOT_PATH/goSNPmania/${SAMPLEID}.h.sapiens.variations.goSNPmania.gz | $SERA_PATH/bin/awkscript/extractHapmapSNPs.awk - | sort -k15 -rn - > $ROOT_PATH/plotAlleleRatio/${SAMPLEID}.h.sapiens.variations.Hapmap.plotAlleleRatio;

				# plot files
				$SERA_PATH/bin/Rscript/plotAlleleRatio_withDensity.R $ROOT_PATH/plotAlleleRatio/${SAMPLEID}.h.sapiens.variations.Hapmap.plotAlleleRatio $ROOT_PATH/plotPdfs/${SAMPLEID}.h.sapiens.hapmap.pdf 20 "${SAMPLEID} $CALL_TYPE ${SNPMANIAFLAGS}";
			fi
		fi
	elif [${CALL_TYPE} == "ampregion" ]; then
		if [[ -e $ROOT_PATH/goSNPmania/${SAMPLEID}.ampregion.variations.goSNPmania.gz ]]; then
			if [[ ! -e $ROOT_PATH/plotPdfs/${SAMPLEID}.ampregion.hapmap.pdf || ! -z $FORCE ]]; then
				# extract hapmap snps 
				zcat $ROOT_PATH/goSNPmania/${SAMPLEID}.ampregion.variations.goSNPmania.gz | $SERA_PATH/bin/awkscript/extractHapmapSNPs.awk - | sort -k15 -rn - > $ROOT_PATH/plotAlleleRatio/${SAMPLEID}.ampregion.variations.Hapmap.plotAlleleRatio;

				# plot files
				$SERA_PATH/bin/Rscript/plotAlleleRatio_withDensity.R $ROOT_PATH/plotAlleleRatio/${SAMPLEID}.ampregion.variations.Hapmap.plotAlleleRatio $ROOT_PATH/plotPdfs/${SAMPLEID}.ampregion.hapmap.pdf 20 "${SAMPLEID} $CALL_TYPE ${SNPMANIAFLAGS}";
			fi
		fi
	fi
fi
if [ ${MOLECULES} == "true" ]; then
	if [ ${CALL_TYPE} == "h.sapiens" ]; then
		if [[ -e $ROOT_PATH/goSNPmania/${SAMPLEID}.h.sapiens.molecules.variations.goSNPmania.gz ]]; then
			if [[ ! -e $ROOT_PATH/plotPdfs/${SAMPLEID}.h.sapiens.molecules.hapmap.pdf || ! -z $FORCE ]]; then
				# extract hapmap snps 
				zcat $ROOT_PATH/goSNPmania/${SAMPLEID}.h.sapiens.molecules.variations.goSNPmania.gz | $SERA_PATH/bin/awkscript/extractHapmapSNPs.awk - | sort -k15 -rn - > $ROOT_PATH/plotAlleleRatio/${SAMPLEID}.h.sapiens.molecules.variations.Hapmap.plotAlleleRatio;

				# plot files
				$SERA_PATH/bin/Rscript/plotAlleleRatio_withDensity.R $ROOT_PATH/plotAlleleRatio/${SAMPLEID}.h.sapiens.molecules.variations.Hapmap.plotAlleleRatio $ROOT_PATH/plotPdfs/${SAMPLEID}.h.sapiens.molecules.hapmap.pdf 20 "${SAMPLEID} $CALL_TYPE ${SNPMANIAFLAGS}";
			fi
		fi
	fi
fi

if [ "$?" != "0" ]; then
	ErrorLog ${SAMPLEID} "Failed in plot allele ratio for hapmap basse";
else
	SuccessLog ${SAMPLEID} "Passed plot allele ratio for hapmap bases, with min read depth 20";
fi
