#!/bin/bash -l

#SBATCH -p core  -n 6
#SBATCH --nodes=1
#SBATCH -t 02:00:00
##SBATCH --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

if [[ $GLOBALS == "MORIARTY" ]]; then
    SNIC_TMP="${ROOT_PATH}/tmp";
    if [[ ! -d ${SNIC_TMP} ]]; then
        mkdir ${SNIC_TMP};
    fi
fi

if [[ $PLATFORM = "Illumina" ]]; then
    # If MATE_PAIR is set to true in the input file
    if [[ "$MATE_PAIR" == "true" ]]; then
        if [[ ${METHOD} == "swift" && ( ${CUTADAPT_PREFIX} == "cp288_masterfile_191114" || ${CUTADAPT_PREFIX} == "Accel-Amplicon-Plus_Lung_Cancer_masterfile" || ${CUTADAPT_PREFIX} == "18-2132_EGFR_MID_Masterfile_mod20191002" ) ]]; then

           #CUTADAPT_PREFIX should point to the masterfile
           singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY sh -s "primerclip ${ROOT_PATH}/refFiles/${CUTADAPT_PREFIX}.txt $ROOT_PATH/Bwa/${SAMPLEID}.untrimmed.qsorted.sam $ROOT_PATH/Bwa/${SAMPLEID}.pclip.sam";
           SuccessLog "${SAMPLEID}" "/projects/wp4/nobackup/workspace/jonas_test/primerclip/primerclip ${ROOT_PATH}/refFiles/${CUTADAPT_PREFIX}.txt $ROOT_PATH/Bwa/${SAMPLEID}.untrimmed.qsorted.sam $ROOT_PATH/Bwa/${SAMPLEID}.pclip.sam;"

           singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY sh -c "samtools view -bS $ROOT_PATH/Bwa/${SAMPLEID}.pclip.sam | samtools sort -@ 3 /dev/stdin -f $ROOT_PATH/Bwa/${SAMPLEID}.sorted.sc.bam";
           singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY sh -c "samtools index $ROOT_PATH/Bwa/${SAMPLEID}.sorted.sc.bam";
           singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY sh -c "java -jar ${SERA_PATH}/bin/java/biostar84452.jar $ROOT_PATH/Bwa/${SAMPLEID}.sorted.sc.bam --samoutputformat BAM > $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam";
           singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY sh -c "samtools index $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam";
           singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY sh -c "samtools flagstat $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam > $ROOT_PATH/Bwa/${SAMPLEID}.alignmentStats.txt";
           SuccessLog "${SAMPLEID}" "samtools view -bS $ROOT_PATH/Bwa/${SAMPLEID}.pclip.sam | samtools sort -@ 3 /dev/stdin -o $ROOT_PATH/Bwa/${SAMPLEID}.sorted.sc.bam;"
           SuccessLog "${SAMPLEID}" "samtools index $ROOT_PATH/Bwa/${SAMPLEID}.sorted.sc.bam;"
           SuccessLog "${SAMPLEID}" "java -jar ${SERA_PATH}/bin/java/biostar84452.jar $ROOT_PATH/Bwa/${SAMPLEID}.sorted.sc.bam --samoutputformat BAM > $ROOT_PATH/Bwa/${SAMPLEID}.sorted.sc.bam;"
           SuccessLog "${SAMPLEID}" "samtools flagstat $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam > $ROOT_PATH/Bwa/${SAMPLEID}.alignmentStats.txt;"
           rm $ROOT_PATH/Bwa/${SAMPLEID}*.sam;
           rm $ROOT_PATH/Bwa/${SAMPLEID}.untrimmed.qsorted.bam;
           rm $ROOT_PATH/Bwa/${SAMPLEID}.sorted.sc.bam*;
        fi
    else
        ErrorLog "${SAMPLEID}" "Only implemented for paired-end sequencing!";
    fi
fi

if [[ "$?" != "0" ]]; then
    ErrorLog "${SAMPLEID}" "Failed in primerclip...";
else
    SuccessLog "${SAMPLEID}" "Passed primerclip for swift+breast otherwise did nothing";
fi
