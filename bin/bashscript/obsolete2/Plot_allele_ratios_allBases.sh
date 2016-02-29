#!/bin/bash
#
# Script to plot allele ratios
#
#SBATCH -p core -n 1
#SBATCH -t 02:00:00
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
                        if [[ ! -e $ROOT_PATH/plotPdfs/${SAMPLEID}.h.sapiens.allBases.pdf || ! -z $FORCE ]]; then
                                # extract allBases snps 
                                zcat $ROOT_PATH/goSNPmania/${SAMPLEID}.h.sapiens.variations.goSNPmania.gz | $SERA_PATH/bin/awkscript/extractAllSNPs.awk - | sort -k15 -n - > $ROOT_PATH/plotAlleleRatio/${SAMPLEID}.h.sapiens.variations.allBases.plotAlleleRatio;

                                # plot files
                                $SERA_PATH/bin/Rscript/plotAlleleRatio_withDensity.R $ROOT_PATH/plotAlleleRatio/${SAMPLEID}.h.sapiens.variations.allBases.plotAlleleRatio $ROOT_PATH/plotPdfs/${SAMPLEID}.h.sapiens.allBases.pdf 20 "${SAMPLEID} $CALL_TYPE ${SNPMANIAFLAGS}";
                        fi
                fi
        elif [${CALL_TYPE} == "ampregion" ]; then
                if [[ -e $ROOT_PATH/goSNPmania/${SAMPLEID}.ampregion.variations.goSNPmania.gz ]]; then
                        if [[ ! -e $ROOT_PATH/plotPdfs/${SAMPLEID}.ampregion.allBases.pdf || ! -z $FORCE ]]; then
                                # extract allBases snps 
                                zcat $ROOT_PATH/goSNPmania/${SAMPLEID}.ampregion.variations.goSNPmania.gz | $SERA_PATH/bin/awkscript/extractAllSNPs.awk - | sort -k15 -n - > $ROOT_PATH/plotAlleleRatio/${SAMPLEID}.ampregion.variations.allBases.plotAlleleRatio;

                                # plot files
                                $SERA_PATH/bin/Rscript/plotAlleleRatio_withDensity.R $ROOT_PATH/plotAlleleRatio/${SAMPLEID}.ampregion.variations.allBases.plotAlleleRatio $ROOT_PATH/plotPdfs/${SAMPLEID}.ampregion.allBases.pdf 20 "${SAMPLEID} $CALL_TYPE ${SNPMANIAFLAGS}";
                        fi
                fi
        fi
fi
if [ ${MOLECULES} == "true" ]; then
        if [ ${CALL_TYPE} == "h.sapiens" ]; then
                if [[ -e $ROOT_PATH/goSNPmania/${SAMPLEID}.h.sapiens.molecules.variations.goSNPmania.gz ]]; then
                        if [[ ! -e $ROOT_PATH/plotPdfs/${SAMPLEID}.h.sapiens.molecules.allBases.pdf || ! -z $FORCE ]]; then
                                # extract allBases snps 
                                zcat $ROOT_PATH/goSNPmania/${SAMPLEID}.h.sapiens.molecules.variations.goSNPmania.gz | $SERA_PATH/bin/awkscript/extractAllSNPs.awk - | sort -k15 -n - > $ROOT_PATH/plotAlleleRatio/${SAMPLEID}.h.sapiens.molecules.variations.allBases.plotAlleleRatio;

                                # plot files
                                $SERA_PATH/bin/Rscript/plotAlleleRatio_withDensity.R $ROOT_PATH/plotAlleleRatio/${SAMPLEID}.h.sapiens.molecules.variations.allBases.plotAlleleRatio $ROOT_PATH/plotPdfs/${SAMPLEID}.h.sapiens.molecules.allBases.pdf 20 "${SAMPLEID} $CALL_TYPE ${SNPMANIAFLAGS}";
                        fi
                fi
        fi
fi
if [ "$?" != "0" ]; then
        ErrorLog ${SAMPLEID} " Failed in plot allele ratio all bases";
else
        SuccessLog ${SAMPLEID} " Passed plot allele ratio all bases, with min read depth 20";
fi

