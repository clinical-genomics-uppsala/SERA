#!/bin/bash
#SBATCH -p core
#SBATCH -t 00:15:00
#SBATCH --qos=short
#
# Script for running MosaikAligner against a genome and save unaligned reads.
#

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog ${REFSEQ} "Starts collecting alignment data...";

# Check if the directory exists, 
if [ ! -d "$ROOT_PATH/Alignments" ]; then
	ErrorLog ${REFSEQ} "Unable to check diretory $ROOT_PATH/Alignments.";
	exit 1;
fi

if [ ! -d "$ROOT_PATH/TextResults" ]; then
	ErrorLog ${REFSEQ} "Unable to check diretory $ROOT_PATH/TextResults";
	exit 1;
fi

# File headers
if [[ ! -e $ROOT_PATH/TextResults/Results.txt || ! -z $FORCE ]]; then
	echo -e "#SampleID\tA.Fraction.R\tUniq.Fraction.R\tFiltered.Out.R\tFail.hash.R\t#Reads.R\tA.Fraction.G\tUniq.Fraction.G\tFiltered.Out.G\tFail.hash.G\t#Reads.G\tSpec.seqRegion\tSpec.seqROI\tSpec.seqRegion.Uniq\tSpec.seqROI.Uniq\tDepth.seqRegion\tDepth.seqROI\tDepth\tDepth.seqRegion.Uniq\tDepth.seqROI.Uniq\tDepth.Uniq" > $ROOT_PATH/TextResults/Results.txt;
fi

# Go through all samples in the inputFile
if [[ ! -e $ROOT_PATH/TextResults/Results.txt || ! -z $FORCE ]]; then

	# loop over all found sample txt files, a underscore in sample id messes this up
	samples=( `ls $ROOT_PATH/Alignments/ | grep "h.sapiens" | grep ".txt" | awk '{split($0,a,"_"); print a[3]}' | awk '{split($0,a,"."); print a[1];}' | xargs echo` );
	for sample in ${samples[@]}; do

		# overwrite old id
		SAMPLEID=$sample;

		# Collect data for region.
		if [[ -n $ALIGNERFLAGS_REGION && -e $ROOT_PATH/Alignments/MosaikAligner_region_${SAMPLEID}.txt ]]; then
			r_fraction_aligned=`grep -m 1 "total aligned:" $ROOT_PATH/Alignments/MosaikAligner_region_${SAMPLEID}.txt | awk '{print (\$6/100);}'`;
			r_fraction_uniq=`grep -m 1 "# unique:" $ROOT_PATH/Alignments/MosaikAligner_region_${SAMPLEID}.txt | awk '{print (\$5/100);}'`;
			r_fraction_nonuniq=`grep -m 1 "# non-unique:" $ROOT_PATH/Alignments/MosaikAligner_region_${SAMPLEID}.txt | awk '{print (\$5/100);}'`;
			r_filtered_out=`grep -m 1 "# filtered out:" $ROOT_PATH/Alignments/MosaikAligner_region_${SAMPLEID}.txt | awk '{print (\$6/100);}'`;
			r_fail_hash=`grep -m 1 "# failed hash:" $ROOT_PATH/Alignments/MosaikAligner_region_${SAMPLEID}.txt | awk '{print (\$6/100);}'`;
			r_totreads=`grep -m 1 "total:" $ROOT_PATH/Alignments/MosaikAligner_region_${SAMPLEID}.txt | awk '{print \$2;}'`;
		else
			r_fraction_aligned="-";
			r_fraction_uniq="-";
			r_fraction_nonuniq="-";
			r_filtered_out="-";
			r_fail_hash="-";
			r_totreads="-";
		fi

		# Collect data for genome.
		g_fraction_aligned=`grep -m 1 "total aligned:" $ROOT_PATH/Alignments/MosaikAligner_h.sapiens_${SAMPLEID}.txt | awk '{print (\$6/100);}'`;
		g_fraction_uniq=`grep -m 1 "# unique:" $ROOT_PATH/Alignments/MosaikAligner_h.sapiens_${SAMPLEID}.txt | awk '{print (\$5/100);}'`;
		g_fraction_nonuniq=`grep -m 1 "# non-unique:" $ROOT_PATH/Alignments/MosaikAligner_h.sapiens_${SAMPLEID}.txt | awk '{print (\$5/100);}'`;
		g_filtered_out=`grep -m 1 "# filtered out:" $ROOT_PATH/Alignments/MosaikAligner_h.sapiens_${SAMPLEID}.txt | awk '{print (\$6/100);}'`;
		g_fail_hash=`grep -m 1 "# failed hash:" $ROOT_PATH/Alignments/MosaikAligner_h.sapiens_${SAMPLEID}.txt | awk '{print (\$6/100);}'`;
		g_totreads=`grep -m 1 "total:" $ROOT_PATH/Alignments/MosaikAligner_h.sapiens_${SAMPLEID}.txt | awk '{print \$2;}'`;

		# Collect generic data.
		if [[ $ANALYZE_NONUNIQUE == "true" ]]; then
			spec_region_all=`grep "SPEC_REGION_ALL" $ROOT_PATH/TextResults/${SAMPLEID}.txt | tail -n 1 | awk 'BEGIN{FS="=";}{print (\$2/100);}'`;
			spec_roi_all=`grep "SPEC_ROI_ALL" $ROOT_PATH/TextResults/${SAMPLEID}.txt | tail -n 1 | awk 'BEGIN{FS="=";}{print (\$2/100);}'`;
		else
			spec_region_all="-";
			spec_roi_all="-";
		fi
	
		spec_region_uniq=`grep "SPEC_REGION_UNIQ" $ROOT_PATH/TextResults/${SAMPLEID}.txt | tail -n 1 | awk 'BEGIN{FS="=";}{print (\$2/100);}'`;
		spec_roi_uniq=`grep "SPEC_ROI_UNIQ" $ROOT_PATH/TextResults/${SAMPLEID}.txt | tail -n 1 | awk 'BEGIN{FS="=";}{print (\$2/100);}'`;

		if [[ $ANALYZE_NONUNIQUE == "true" ]]; then
			depth_seqregion_all=`grep "DEPTH_SEQREGION_ALL" $ROOT_PATH/TextResults/${SAMPLEID}.txt | tail -n 1 | awk 'BEGIN{FS="=";}{print \$2;}'`;
			depth_seqroi_all=`grep "DEPTH_SEQROI_ALL" $ROOT_PATH/TextResults/${SAMPLEID}.txt | tail -n 1 | awk 'BEGIN{FS="=";}{print \$2;}'`;
			depth_all=`grep "DEPTH_ALL" $ROOT_PATH/TextResults/${SAMPLEID}.txt | tail -n 1 | awk 'BEGIN{FS="=";}{print \$2;}'`;
		else
			depth_seqregion_all="-";
			depth_seqroi_all="-";
			depth_all="-";
		fi

		depth_seqregion_uniq=`grep "DEPTH_SEQREGION_UNIQ" $ROOT_PATH/TextResults/${SAMPLEID}.txt | tail -n 1 | awk 'BEGIN{FS="=";}{print \$2;}'`;
		depth_seqroi_uniq=`grep "DEPTH_SEQROI_UNIQ" $ROOT_PATH/TextResults/${SAMPLEID}.txt | tail -n 1 | awk 'BEGIN{FS="=";}{print \$2;}'`;
		depth_all_uniq=`grep "DEPTH_UNIQ" $ROOT_PATH/TextResults/${SAMPLEID}.txt | tail -n 1 | awk 'BEGIN{FS="=";}{print \$2;}'`;

		# Write results.
		echo -e "${SAMPLEID}\t$r_fraction_aligned\t$r_fraction_uniq\t$r_filtered_out\t$r_fail_hash\t$r_totreads\t$g_fraction_aligned\t$g_fraction_uniq\t$g_filtered_out\t$g_fail_hash\t$g_totreads\t$spec_region_all\t$spec_roi_all\t$spec_region_uniq\t$spec_roi_uniq\t$depth_seqregion_all\t$depth_seqroi_all\t$depth_all\t$depth_seqregion_uniq\t$depth_seqroi_uniq\t$depth_all_uniq" >> $ROOT_PATH/TextResults/Results.txt;

	done;

fi
