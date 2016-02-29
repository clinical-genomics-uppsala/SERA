#!/bin/bash
#
# Script for running ensembl variant effect prediction from SNPMania output
#
#SBATCH -p node -n 1
#SBATCH -t 10:00:00
##SBATCH -p core -t 00:15:00 --qos=short

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/effectPrediction" ]; then
	mkdir $ROOT_PATH/effectPrediction;
fi

function prepareFiles {

	# transform files to variant format
#	$SERA_PATH/bin/perlscript/snpManiaToVariationFile.pl -i $ROOT_PATH/SNPmania/${1}.${CALL_TYPE}.insertions -d $ROOT_PATH/SNPmania/${1}.${CALL_TYPE}.deletions -v $ROOT_PATH/SNPmania/${1}.${CALL_TYPE}.variations -t 0.1 -m 0.1 -n 0.2 -h 20 -r $HAPMAP_SNP_REF > $ROOT_PATH/effectPrediction/${1}.${CALL_TYPE}.variations.vep;
	
	# very temporary for just variations with MAFs considered
	$HOME/testHash/snpManiaToVariationFile.pl -v $ROOT_PATH/SNPmania/${1}.${CALL_TYPE}.variations -t 0.1 -m 0.1 -n 0.2 -h 20 -r $HAPMAP_SNP_REF -a 0.05 -b $HOME/testHash/1kg.frq > $ROOT_PATH/effectPrediction/${1}.${CALL_TYPE}.variations.vep;

#	$HOME/program/SERA/bin/perlscript/snpManiaToVariationFile.pl -v $HOME/glob/private/b2011080/SNPmania/S98.h.sapiens.uniq.variations -i  -t 0.1 -m 0.1 -n 0.2 -h 20 -r /bubo/home/h20/elinfalk/data/HapMap/NA12802.genotypes_CEU_r27_nr.b36_fwd.lifted2hg19.txt -a 0.05 -b 1kg.frq > S98.vep

}

# make sure that the normal sampleid exists before carrying on
if [[ ! -s "$ROOT_PATH/effectPrediction/${SAMPLEID}.variations.predicted.novel" || ! -z $FORCE ]]; then

	if [[ ! -e "$ROOT_PATH/effectPrediction/${NORMAL_SAMPLEID}.${CALL_TYPE}.variations.vep" && ${NORMAL_SAMPLEID} != "false" ]] || [[ ! -z $FORCE ]]; then

		SuccessLog $SAMPLEID "Preparing normal file...";

		# predict effects of normal sampleid
		prepareFiles $NORMAL_SAMPLEID;

	fi

	SuccessLog $SAMPLEID "Starting to predict effect of SNP and indels...";

	# prepare sample file input
	prepareFiles $SAMPLEID;

	# remove normal variations from tumor/normal file (if set)
	if [ ${NORMAL_SAMPLEID} != "false" ]; then

		# print tumor file and normal file twice, sort and uniq to sort out all variations in normal file from tumor
		cp $ROOT_PATH/effectPrediction/${SAMPLEID}.${CALL_TYPE}.variations.vep $ROOT_PATH/effectPrediction/${SAMPLEID}.${CALL_TYPE}.variations.vep.raw;
		cat $ROOT_PATH/effectPrediction/${SAMPLEID}.${CALL_TYPE}.variations.vep.raw $ROOT_PATH/effectPrediction/${NORMAL_SAMPLEID}.${CALL_TYPE}.variations.vep $ROOT_PATH/effectPrediction/${NORMAL_SAMPLEID}.${CALL_TYPE}.variations.vep | sort | rev | uniq -f1 -u | rev > $ROOT_PATH/effectPrediction/${SAMPLEID}.${CALL_TYPE}.variations.vep;

		SuccessLog $SAMPLEID "Mutual variations in tumor/normal removed";

	fi

	# run variant prediction, should test with --most_severe and --summary --regulatory 
	$VARIANTPREDICTION -i $ROOT_PATH/effectPrediction/${SAMPLEID}.${CALL_TYPE}.variations.vep -o $ROOT_PATH/effectPrediction/${SAMPLEID}.${CALL_TYPE}.variations.predicted --terms ensembl --gene --check_existing --hgnc --force_overwrite --quiet;

	# remove copies
	cp $ROOT_PATH/effectPrediction/${SAMPLEID}.${CALL_TYPE}.variations.predicted $ROOT_PATH/effectPrediction/${SAMPLEID}.${CALL_TYPE}.variations.predicted.raw;
	cat $ROOT_PATH/effectPrediction/${SAMPLEID}.${CALL_TYPE}.variations.predicted.raw | sort | uniq > $ROOT_PATH/effectPrediction/${SAMPLEID}.${CALL_TYPE}.variations.predicted;

	if [ "$?" != "0" ]; then
		ErrorLog $SAMPLEID "Failed in predicting variants";
	else
		SuccessLog $SAMPLEID "Passed predicting variants";
	fi

else
	SuccessLog $SAMPLEID "Prediction files already exists, skipping step";
fi

