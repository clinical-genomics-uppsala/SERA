#!/bin/bash
#
# Script for running MosaikAligner against a genome and save unaligned reads.
#
#SBATCH -p devcore  -n 8
#SBATCH -t 02:00:00

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/Alignments" ]; then
	mkdir $ROOT_PATH/Alignments;
fi

FLAGS=${ALIGNERFLAGS_GENOME};

# Should we use paired end mapping (PCR design)?
if [ "$MATE_PAIR" == "true" ]; then
	SuccessLog "${SAMPLEID}" "Aligning against the genome using mate-pairs/pair-end algorithm...";
else
	SuccessLog "${SAMPLEID}" "Aligning against the genome using single reads algorithm...";
fi

# New alignment?
if [[ ! -e $ROOT_PATH/Alignments/${SAMPLEID}.h.sapiens.aligned.dat || ! -z $FORCE ]]; then

	# Run MosaikAligner and save unaligned reads in fastq format
	if [ $PLATFORM = "Illumina" ]; then
		if [ -e $ROOT_PATH/MosaikBuild/${SAMPLEID}.dat ]; then
			$ROOT_PATH_MOSAIK/MosaikAligner -in $ROOT_PATH/MosaikBuild/${SAMPLEID}.dat -out $ROOT_PATH/Alignments/${SAMPLEID}.h.sapiens.aligned -ia $GENOME_REF.dat -j $GENOME_REF"_15" -annse $ROOT_PATH_MOSAIK/2.1.26.se.100.005.ann -annpe $ROOT_PATH_MOSAIK/2.1.26.pe.100.0065.ann $FLAGS > $ROOT_PATH/Alignments/MosaikAligner_${SAMPLEID}_h.sapiens.txt ;
		else
			ErrorLog "${SAMPLEID} MosaikBuild file does not exist!";
		fi
	else
		ErrorLog "${SAMPLEID}" "Platform $PLATFORM not recognized (only Illumina supported at the moment...).";
		exit 1;
	fi

	# Check if MosaikAligner worked
	if [ "$?" != "0" ]; then
		ErrorLog "${SAMPLEID}" "Failed in MosaikAligner against genome.";
	else
		SuccessLog "${SAMPLEID}" "Passed MosaikAligner against genome."

	fi

else
	ErrorLog "${SAMPLEID}" "MosaikAligner against genome already exists and force was NOT used!.";
fi
