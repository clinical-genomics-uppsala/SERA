#!/bin/bash
#
# Script to export aligment archives as BAM-files.
#
#SBATCH -p core -n 3
#SBATCH -t 03:00:00

# Include functions
. $SERA_PATH/includes/logging.sh;

# this function extracts the reads in ampregion from the genome alignment
#function extractAmpregion {
#	samtools sort $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens${1}.corrSel.bam $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens${1}.corrSel.sorted;
#	samtools index $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens${1}.corrSel.sorted.bam;
#	samtools view -b -h $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens${1}.corrSel.sorted.bam `$SERA_PATH/bin/perlscript/convertAmpregionRefID.pl $ROOT_PATH/refFiles/${REFSEQ}.ampregion $CORRSEL_REF_CONV_FILE` -o $ROOT_PATH/SamBamFiles/${SAMPLEID}.ampregion${1}.corrSel.bam;
#	rm $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens${1}.corrSel.sorted*;
#}

# Running Mate Pairs or Not.
if [ "$MATE_PAIR" == "true" ]; then
	PE_FLAGS="-p";
else
	PE_FLAGS="-s";
fi

# MDA or PCR
if [ "$SEQUNCING_TAG" != "false" ]; then
	PE_FLAGS="${PE_FLAGS}";
else
	PE_FLAGS="-j";
fi

SuccessLog "${SAMPLEID}" "Starting CorrSel on BAM files...";

# run with ampregion file (if available)
if [ -e $ROOT_PATH/SamBamFiles/${SAMPLEID}.ampregion.uniq.bam ]; then

	if [[ $ANALYZE_NONUNIQUE = "true" ]]; then

		### Run CorrSel ALL READS
		if [[ ! -e $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.corrSel.bam || ! -z $FORCE ]]; then

			# extract h.sapiens alignment
			$ROOT_PATH_CORRSEL/convertSelectorAlignment.sh $PE_FLAGS -f -m -rf $CORRSEL_REF_CONV_FILE -g $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.bam -r $ROOT_PATH/SamBamFiles/${SAMPLEID}.ampregion.bam -o $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.corrSel.bam >> $ROOT_PATH/SamBamFiles/CorrSel_${SAMPLEID}.txt;
		
			# extract only ampregion reads
#			extractAmpregion "";
		fi

	fi

	### Run CorrSel UNIQUE READS
	if [[ ! -e $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.uniq.corrSel.bam || -z $FORCE ]]; then

		$ROOT_PATH_CORRSEL/convertSelectorAlignment.sh $PE_FLAGS -f -m -rf $CORRSEL_REF_CONV_FILE -g $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.uniq.bam -r $ROOT_PATH/SamBamFiles/${SAMPLEID}.ampregion.uniq.bam -o $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.uniq.corrSel.bam >> $ROOT_PATH/SamBamFiles/CorrSel_${SAMPLEID}.txt;
#		extractAmpregion ".uniq";
	fi

	SuccessLog "${SAMPLEID}" "CorrSel on BAM files (using '${PE_FLAGS} -f -m') done...";

# run without ampregion file (ampregion missing!)	
else
	if [[ $ANALYZE_NONUNIQUE = "true" ]]; then

		if [[ ! -e $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.corrSel.bam || ! -z $FORCE ]]; then	
			$ROOT_PATH_CORRSEL/convertSelectorAlignment.sh $PE_FLAGS -f -c -rf $CORRSEL_REF_CONV_FILE -g $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.bam -o $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.corrSel.bam >> $ROOT_PATH/SamBamFiles/CorrSel_${SAMPLEID}.txt;
		fi

	fi

	if [[ ! -e $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.uniq.corrSel.bam || ! -z $FORCE ]]; then
		$ROOT_PATH_CORRSEL/convertSelectorAlignment.sh $PE_FLAGS -f -c -rf $CORRSEL_REF_CONV_FILE -g $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.uniq.bam -o $ROOT_PATH/SamBamFiles/${SAMPLEID}.h.sapiens.uniq.corrSel.bam >> $ROOT_PATH/SamBamFiles/CorrSel_${SAMPLEID}.txt;
	fi

	SuccessLog "${SAMPLEID}" "No need for CorrSel to merge two files. Done exporting to genomic position only using '$PE_FLAGS -f -c'.";
fi
