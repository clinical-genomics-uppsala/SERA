#!/bin/bash
#
# Script creats BlastDB on mate-pair reads.
#
#SBATCH -p core -n 1
#SBATCH -t 30:00
##SBATCH --qos=short
##module load bioinfo-tools blast/2.2.25

# Include functions
. $SERA_PATH/includes/logging.sh;

if [ ! -d $ROOT_PATH/cnvFromBLAST ]; then
	mkdir $ROOT_PATH/cnvFromBLAST
fi

# Check if BLAST analysis is possible.
if [ $PLATFORM == "Illumina" ]; then
	if [ ${READS} == "true" ]; then
		# Create file with log2ratio
		if [[ ${NORMAL_SAMPLEID}!="FALSE" && -e $ROOT_PATH/BLAST/${SAMPLEID}.blast.hits.txt && -e $ROOT_PATH/BLAST/${NORMAL_SAMPLEID}.blast.hits.txt && ! -e $ROOT_PATH/cnvFromBLAST/${SAMPLEID}"_"${NORMAL_SAMPLEID}".blast.log2ratio" ]]; then
			perl $SERA_PATH/bin/perlscript/cnvAnalysisFromBlastAnalysis.pl -n $ROOT_PATH/BLAST/${NORMAL_SAMPLEID}.blast.hits.txt -t $ROOT_PATH/BLAST/${SAMPLEID}.blast.hits.txt -nMinRD 20 -s ${SELECTIONFILE} | sort -k2,3V > $ROOT_PATH/cnvFromBLAST/${SAMPLEID}"_"${NORMAL_SAMPLEID}".blast.log2ratio";
		fi
	
		#Check if sorted ampregion file exists
		if [ ! -e $ROOT_PATH/refFiles/${REFSEQ}.sorted.ampregion ]; then
			if [ -e $ROOT_PATH/refFiles/${REFSEQ}.ampregion ]; then
				sort -k2,3V $ROOT_PATH/refFiles/${REFSEQ}.ampregion > $ROOT_PATH/refFiles/${REFSEQ}.sorted.ampregion;
			else
				$ROOT_PATH_SEDD/mergeregions.sh ${SELECTIONFILE} | sort -k2,3V - > $ROOT_PATH/refFiles/${REFSEQ}.sorted.ampregion
			fi
		fi

		if [ ! -e $ROOT_PATH/refFiles/${REFSEQ}.enzymeCombinations.txt ]; then
			awk '{split($1,a,"|"); print a[2]}' $ROOT_PATH/refFiles/${REFSEQ}.selection | sort -u > $ROOT_PATH/refFiles/${REFSEQ}.enzymeCombinations.txt
		fi

		#Create cnv plot
		if [[ -e $ROOT_PATH/cnvFromBLAST/${SAMPLEID}_${NORMAL_SAMPLEID}.blast.log2ratio && ! -e $ROOT_PATH/cnvFromBLAST/${SAMPLEID}_${NORMAL_SAMPLEID}.blast.log2ratio.perRE.pdf ]]; then
			#SERA_scripts/bin/Rscript/cnvPlotPerFragment_perRE_inOnePlot.R $ROOT_PATH/cnvFromBLAST/${SAMPLEID}_${NORMAL_SAMPLEID}.blast.log2ratio $ROOT_PATH/refFiles/${REFSEQ}.sorted.ampregion $ROOT_PATH/refFiles/${REFSEQ}.enzymeCombinations.txt $ROOT_PATH/cnvFromBLAST/${SAMPLEID}_${NORMAL_SAMPLEID}.blast.log2ratio.perRE.pdf;

			$SERA_PATH/bin/Rscript/cnvPlotPerFragment_perRE.R $ROOT_PATH/cnvFromBLAST/${SAMPLEID}_${NORMAL_SAMPLEID}.blast.log2ratio $ROOT_PATH/refFiles/${REFSEQ}.sorted.ampregion $ROOT_PATH/refFiles/${REFSEQ}.enzymeCombinations.txt $ROOT_PATH/cnvFromBLAST/${SAMPLEID}_${NORMAL_SAMPLEID}.blast.log2ratio.perRE.pdf;
		else 
			ErrorLog "Could not create log2ratio files!";
		fi
	fi

fi

SuccessLog "${SAMPLEID}" "${0##*/}" "Created log2(T/N).";
