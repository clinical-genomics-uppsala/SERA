#!/bin/bash -l
#SBATCH -p core -n 4
#SBATCH -t 03:00:00

# Include functions
. $SERA_PATH/includes/logging.sh

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/MosaikRef" ]; then
	mkdir $ROOT_PATH/MosaikRef;
fi

# check which output file should be expected
if [[ "$SEQUNCING_TAG" == "false" ]]&&[[ "$DESIGN_TYPE" == "MDA" ]]; then
	REGION=${REFSEQ};
elif [[ "$SEQUNCING_TAG" != "false" ]]&&[[ "$DESIGN_TYPE" == "PCR" ]]; then
	REGION=${REFSEQ}_${SEQUNCING_TAG};
elif [[ "$SEQUNCING_TAG" == "false" ]]; then
	REGION=${REFSEQ};
fi

# check if we have an output file
if [[ ! -s $ROOT_PATH/MosaikRef/${REGION}.dat || ! -z $FORCE ]]; then

	# if the sequence have not been downloaded, download it!
	if [[ ! -s $ROOT_PATH/MosaikRef/${REFSEQ}.selection.seq.gz  || ! -z $FORCE ]]; then
		perl $DOWNLOAD2FASTA -th 4 -t full -s $ROOT_PATH/refFiles/${REFSEQ}.selection -o - -d $BLAST_DB | sort -k2,4 | gzip > $ROOT_PATH/MosaikRef/${REFSEQ}.selection.seq.gz;
	fi
		
	# add junctions
	if [[ "$SEQUNCING_TAG" == "false" ]]&&[[ "$DESIGN_TYPE" == "MDA" ]]&&[[ ! -s $ROOT_PATH/MosaikRef/${REFSEQ}.fasta.gz || ! -z $FORCE ]]; then

		SuccessLog "${REFSEQ}" "Creating artificial sequence as reference (MDA design)...";

		$SERA_PATH/bin/awkscript/ampregion.awk $ROOT_PATH/refFiles/${REFSEQ}.ampregion | perl $DOWNLOAD2FASTA -th 1 -t small -s /dev/stdin -o /dev/stdout -d $BLASTDB | awk '{print $1"\n"$2}' > $ROOT_PATH/MosaikRef/${REFSEQ}.ampregion.fasta;

		zcat $ROOT_PATH/MosaikRef/${REFSEQ}.selection.seq.gz | $SERA_PATH/bin/awkscript/selection.awk -v mm=$NUMBER_OF_MM r=$READ_LENGTH /dev/stdin >> $ROOT_PATH/MosaikRef/${REFSEQ[${i}]}.fasta

#		zcat $ROOT_PATH/MosaikRef/${REFSEQ}.selection.seq.gz | awk -f $SERA_PATH/bin/awkscript/addArtificialJoints2Reference.awk | gzip > $ROOT_PATH/MosaikRef/${REFSEQ}.fasta.gz

	# add sequencing tags
	elif [[ "$SEQUNCING_TAG" != "false" ]]&&[[ "$DESIGN_TYPE" == "PCR" ]]&&[[ ! -s $ROOT_PATH/MosaikRef/${REFSEQ}_${SEQUNCING_TAG}.fasta.gz || ! -z $FORCE ]]; then

		SuccessLog "${REFSEQ}" "Adding sequencing index ${SEQUNCING_TAG} to reference (PCR design)...";
		. $SERA_PATH/config/sequencingTags.sh
		zcat $ROOT_PATH/MosaikRef/${REFSEQ}.selection.seq.gz | awk -v r=${READ_LENGTH} -v fivePtag="${fTag[${SEQUNCING_TAG}]}" -v threePtag="${tTag[${SEQUNCING_TAG}]}" -f $SERA_PATH/bin/awkscript/addSequencingTags2Reference.awk | gzip > $ROOT_PATH/MosaikRef/${REFSEQ}_${SEQUNCING_TAG}.fasta.gz;

	fi

	# Check if creating file worked
	if [ "$?" != "0" ]; then
		ErrorLog "${REFSEQ}" "Failed in create amplified region file.";
	else
		SuccessLog "${REFSEQ}" "Passed create amplified region file.";
	fi


	# build mosaik file
	if [[ $(zcat $ROOT_PATH/MosaikRef/${REGION}.fasta.gz | grep ">" | wc -l) != 0 ]]; then

		# Build Mosaik reference file.
		MosaikBuild -fr $ROOT_PATH/MosaikRef/${REGION}.fasta.gz -oa $ROOT_PATH/MosaikRef/${REGION}.dat;

		if [ $PLATFORM = "SOLiD" ]; then
		   MosaikBuild -fr $ROOT_PATH/MosaikRef/${REGION}.fasta.gz -oa $ROOT_PATH/MosaikRef/${REGION}.CS.dat -cs;
		fi

		SuccessLog $SAMPLEID "Passed MosaikBuild.";

	else
		ErrorLog ${REFSEQ} "Something is wrong with the reference fasta file.";
	fi
else
	SuccessLog $SAMPLEID "Reference files for Mosaik already created, no need to run step for this sample";
fi





