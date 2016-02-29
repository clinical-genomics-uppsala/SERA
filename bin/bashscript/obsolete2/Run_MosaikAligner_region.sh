#!/bin/bash
#
# Script for running MosaikAligner against a region and save unaligned reads.
#
#SBATCH -p core -n 5
#SBATCH -t 10:00:00


# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/Alignments" ]; then
	mkdir $ROOT_PATH/Alignments;
fi

FLAGS=${ALIGNERFLAGS_REGION};

# Should we use paired end mapping (PCR design)?
if [ "$MATE_PAIR" == "true" ]; then
	SuccessLog "${SAMPLEID}" "Aligning against region ${REFSEQ} using mate-pair algorithm...";
else
	SuccessLog "${SAMPLEID}" "Aligning against region ${REFSEQ} using single reads algorithm...";
fi

# Is it a new alignment?
if [[ ! -e $ROOT_PATH/Alignments/${SAMPLEID}.ampregion.aligned.dat || ! -z $FORCE ]]; then

	# Do we have a region to use? 
	if [[ -e $ROOT_PATH/MosaikRef/${REGION}.dat ]]; then

		if [ $PLATFORM = "Illumina" ]; then
			if [ -e $ROOT_PATH/MosaikBuild/${SAMPLEID}.dat ]; then
				$ROOT_PATH_MOSAIK/MosaikAligner -in $ROOT_PATH/MosaikBuild/${SAMPLEID}.dat -out $ROOT_PATH/Alignments/${SAMPLEID}.ampregion.aligned -ia $ROOT_PATH/MosaikRef/${REGION}.dat -annse $ROOT_PATH_MOSAIK/2.1.26.se.100.005.ann -annpe $ROOT_PATH_MOSAIK/2.1.26.pe.100.0065.ann $FLAGS > $ROOT_PATH/Alignments/MosaikAligner_${SAMPLEID}_region.txt ;
			else
				ErrorLog "${SAMPLEID} MosaikBuild file does not exist!";
			fi
		else
			ErrorLog "${SAMPLEID}" "Platform $PLATFORM not recognized \(should be Illumina\).";
			exit 1;
		fi

		# Check if MosaikAligner worked
		if [ "$?" != "0" ]; then
			ErrorLog "${SAMPLEID}" "Failed in MosaikAligner against region.";
		fi

    else
		ErrorLog "${SAMPLEID}" "No region reference found, skipping this step.";
    fi

else
	ErrorLog "${SAMPLEID}" "MosaikAligner against region already exists and force was NOT used!.";
fi
