#!/bin/bash
#
# Script creats BlastDB on mate-pair reads.
#
#SBATCH -p core -n 1
#SBATCH -t 15:00
##SBATCH --qos=short

# Include functions
. $SERA_PATH/includes/logging.sh;

if [ ! -d $ROOT_PATH/readsVSmolecules ]; then
	mkdir $ROOT_PATH/readsVSmolecules
fi

# Check that platform is set to Illumina.
if [ ${PLATFORM} == "Illumina" ]; then
	# Check that both reads and molecules are set true, otherwise print error message
	if [[ ${READS} == "true" && ${MOLECULES} == "true" ]]; then
		if [ ${CALL_TYPE} == "h.sapiens" ]; then
			# Check if the inputfile exists
			if [ -e $ROOT_PATH/readsVSmolecules/${SAMPLEID}.h.sapiens.aligned.amplicon.bed ]; then
				# Plot reads vs molecules
				${SERA_PATH}/bin/Rscript/ReadsPerMolecule_vs_readDepth_bam.R $ROOT_PATH/readsVSmolecules/${SAMPLEID}.h.sapiens.aligned.amplicon.bed  1 $ROOT_PATH/readsVSmolecules/${SAMPLEID}.readsPerMolecule_fromBam.h.sapiens.pdf ${SAMPLEID};
			else
				ErrorLog ${SAMPLEID} "Input file $ROOT_PATH/readsVSmolecules/${SAMPLEID}.h.sapiens.aligned.amplicon.bed does not exist!";
			fi
		else
			ErrorLog ${SAMPLEID} "So far only supported for h.sapiens!";
		fi
	else
		ErrorLog ${SAMPLEID} "Both READS and MOLECULES have to be set to true in the input file!";
	fi
else
	ErrorLog ${SAMPLEID} "So far only supported for Illumina!";
fi
# Check if readsVSmolecules worked
if [ "$?" != "0" ]; then
	ErrorLog "${SAMPLEID}" "Failed in readsVSmolecules...";
else
	SuccessLog "${SAMPLEID}" "Passed readsVSmolecules...";
fi
