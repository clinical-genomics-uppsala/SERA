#!/bin/bash
#
# Script creats BlastDB on mate-pair reads.
#
#SBATCH -p core  -n 8 -N 1-1
#SBATCH -t 02:00:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

##SBATCH --qos=short

. $SERA_PATH/includes/load_modules.sh


# Include functions
. $SERA_PATH/includes/logging.sh;

if [[ ${METHOD} == "haloplex" ]]; then
    if [[ ! -d $ROOT_PATH/AmpliconCoverage ]]; then
    	mkdir $ROOT_PATH/AmpliconCoverage
    fi

    if [[ ! -d $ROOT_PATH/AmpliconMapped ]]; then
    	mkdir $ROOT_PATH/AmpliconMapped
    fi

    # Check that platform is set to Illumina.
    if [[ ${PLATFORM} == "Illumina" ]]; then
    	# Check if the selection file exists in bed-format, if not create
    	if [[ ! -e $ROOT_PATH/refFiles/$REFSEQ.selection.bed || ! -z $FORCE ]]; then
    		singularity exec -B /data -B /scratch -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY awk -f ${SERA_PATH}/bin/awkscript/sedd2bed.awk -v name="$REFSEQ" -v desc="Selection file" $NC2chr $ROOT_PATH/refFiles/$REFSEQ.selection > $ROOT_PATH/refFiles/$REFSEQ.selection.bed;
    	else
    		SuccessLog "${SAMPLEID}" "$ROOT_PATH/refFiles/$REFSEQ.selection.bed already exists and force was not used!";
    	fi
    	# Check that reads are true and molecules are false
    	if [[ ${READS} == "true" && ${MOLECULES} == "false" ]]; then
            # Check if the output file exists
            if [[ ! -e $ROOT_PATH/AmpliconMapped/${SAMPLEID}.aligned.withAmplicon.bam || ! -z $FORCE ]]; then
                # Check if the aligned file from bwa exists
                if [[ -e $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam ]]; then
                    # Querysort bam-file
                    singularity exec -B /data -B /scratch -B /opt -B /beegfs-storage -B /projects $SERA_SINGULARITY sh -c "samtools view -h -b -F 0x100 $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam | samtools sort -n -@ 8 /dev/stdin -o $ROOT_PATH/AmpliconMapped/${SAMPLEID}.querysorted.bam"
                    if [[ "$?" != "0" ]]; then
			                  ErrorLog "$SAMPLEID" "Failed samtools sort";
		                else
  		                  SuccessLog "${SAMPLEID}" "samtools view -h -b -F 0x100 $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam | samtools sort -n -@ 8 /dev/stdin $ROOT_PATH/AmpliconMapped/${SAMPLEID}.querysorted";
		                fi
                    # Run ampliconmapping
                    singularity exec -B /data -B /scratch -B /opt -B /beegfs-storage -B /projects $AMPLICONMAPPING_SINGULARITY java -Xmx8g -jar /jar/GenomeAnalysisTKLite_molecules.jar -T MapReadToAmpliconsIlluminaReadPair -R $GENOME_FASTA_REF -I $ROOT_PATH/AmpliconMapped/${SAMPLEID}.querysorted.bam -o $ROOT_PATH/AmpliconCoverage/${SAMPLEID}.amplicon.bed -fragments $ROOT_PATH/refFiles/$REFSEQ.selection.bed -ampAnReads $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.unsorted.bam -U ALL -nonunique -allowPotentiallyMisencodedQuals --downsample_to_coverage 90000 -molBarCode 0
                    if [[ "$?" != "0" ]]; then
                        ErrorLog "$SAMPLEID" "Amplicon mapping failed";
                    else
   		                  SuccessLog "${SAMPLEID}" "java -Xmx8g -jar ${SERA_PATH}/bin/java/GenomeAnalysisTKLite_molecules.jar -T MapReadToAmpliconsIlluminaReadPair -R $GENOME_FASTA_REF -I $ROOT_PATH/AmpliconMapped/${SAMPLEID}.querysorted.bam -o $ROOT_PATH/AmpliconCoverage/${SAMPLEID}.amplicon.bed -fragments $ROOT_PATH/refFiles/$REFSEQ.selection.bed -ampAnReads $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.unsorted.bam -U ALL -nonunique -allowPotentiallyMisencodedQuals --downsample_to_coverage 90000 -molBarCode 0";
		                fi
        
                    singularity exec -B /data -B /scratch -B /opt -B /beegfs-storage -B /projects $SERA_SINGULARITY samtools sort -@ 8 $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.unsorted.bam -o $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.bam;
                    if [[ "$?" != "0" ]]; then
                        ErrorLog "$SAMPLEID" "Failed samtools sort";
                    else
                        SuccessLog "${SAMPLEID}" "samtools sort -@ 8 $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.unsorted.bam -o $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.bam;";
		                fi
                    
                    singularity exec -B /data -B /scratch -B /opt -B /beegfs-storage -B /projects $SERA_SINGULARITY samtools index $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.bam;
                    if [[ "$?" != "0" ]]; then
                        ErrorLog "$SAMPLEID" "Failed samtools index";
                    else
                        SuccessLog "${SAMPLEID}" "samtools index $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.bam;";
		                fi
                    
                    singularity exec -B /data -B /scratch -B /opt -B /beegfs-storage -B /projects $SERA_SINGULARITY samtools flagstat $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.bam > $ROOT_PATH/Bwa/${SAMPLEID}.alignmentStats_noDuplicateReads.txt;
            		    if [[ "$?" != "0" ]]; then
                        ErrorLog "$SAMPLEID" "Failed in running flagstat";
                    else
                        SuccessLog "${SAMPLEID}" "samtools flagstat $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.bam > $ROOT_PATH/Bwa/${SAMPLEID}.alignmentStats_noDuplicateReads.txt";
                    fi

                    # Remove the querysorted file
                    rm $ROOT_PATH/AmpliconMapped/${SAMPLEID}.querysorted.bam;
                    rm $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.unsorted.bam;

                else
                    ErrorLog "${SAMPLEID}" "$ROOT_PATH/AmpliconMapped/${SAMPLEID}.sorted.bam does not exist!";
                fi
            else
                ErrorLog "${SAMPLEID}" "$ROOT_PATH/AmpliconMapped/${SAMPLEID}.aligned.withAmplicon.bam already exists and force was not used!";
            fi
        else
            ErrorLog "${SAMPLEID}" "READS has to be true and MOLECULES false in the input file for the analysis to run!";
        fi
    else
    	ErrorLog ${SAMPLEID} "So far only supported for Illumina!";
    fi
elif [[ ${METHOD} == "swift" ]]; then
    SuccessLog "${SAMPLEID}" "Ampliconmapping is not needed for swift samples - skipping!";

else
    ErrorLog "${SAMPLEID}" "halo and swift are the only methods supported so far!";
fi
# Check if readsVSmolecules worked
if [[ "$?" != "0" ]]; then
	ErrorLog "${SAMPLEID}" "Failed in ampliconMapping...";
else
	SuccessLog "${SAMPLEID}" "Passed ampliconMapping...";
fi
