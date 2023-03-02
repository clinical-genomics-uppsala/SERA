#!/bin/bash -l

#SBATCH -p core  -n 6
#SBATCH -t 01:00:00
##SBATCH --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

export TRIM_LOG="${ROOT_PATH}/seqdata/${SAMPLEID}.trim.log"
PREFIX="${SNIC_TMP}/${SAMPLEID}";

cputhreads=12;



ptrim()
{
    fqt1=$1
    fqt2=$2
    tprefix=${fqt1%%_R1_001.fastq}

    cutadaptFile5prim="${ROOT_PATH}/refFiles/${CUTADAPT_PREFIX}_5ptrim.fa";

    #5’trim
    singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY cutadapt \
        -g file:$cutadaptFile5prim \
        -o ${tprefix}_tmpR1.fq -p ${tprefix}_tmpR2.fq \
    $fqt1 $fqt2 --minimum-length 40 -e 0.12 >> $TRIM_LOG;

    #5’trim
    singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY cutadapt \
        -g file:$cutadaptFile5prim \
        -o ${tprefix}_5ptmpR2.fq -p ${tprefix}_5ptmpR1.fq \
        ${tprefix}_tmpR2.fq ${tprefix}_tmpR1.fq --minimum-length 40 -e 0.12 >> $TRIM_LOG;

}

r1reformat()
{
    tr '~' '\n' < $1 > ${1}_R1_001.fastq
}

r2reformat()
{
    tr '~' '\n' < $1 > ${1}_R2_001.fastq
}

export -f ptrim
export -f r1reformat
export -f r2reformat

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/seqdata" ]]; then
    mkdir $ROOT_PATH/seqdata;
fi


fastq_files_r1=($(echo "$RAWDATA_PE1" | tr " " "\n"));

if [[ ${#fastq_files_r1[@]} > 1 ]];
then
    if [ -n "$(find ${ROOT_PATH}/seqdata -name ${SAMPLEID}_S*_L000_R1_001.fastq.gz | head -1)" ];
    then
        READ1=$(ls ${ROOT_PATH}/seqdata/${SAMPLEID}_S*_L000_R1_001.fastq.gz);
    else
        ErrorLog "${SAMPLEID}" "Multiple lanes for sample ${SAMPLEID}, read1, please pre-process data!...";
    fi
    if [ -n "$(find ${ROOT_PATH}/seqdata -name ${SAMPLEID}_S*_L000_R2_001.fastq.gz | head -1)" ];
    then
        READ2=$(ls ${ROOT_PATH}/seqdata/${SAMPLEID}_S*_L000_R2_001.fastq.gz);
    else
        ErrorLog "${SAMPLEID}" "Multiple lanes for sample ${SAMPLEID}, read2, please pre-process data!...";
    fi
else
    READ1=$RAWDATA_PE1;
    READ2=$RAWDATA_PE2;
fi

if [[ $PLATFORM = "Illumina" ]]; then
    # If MATE_PAIR is set to true in the input file
    if [[ "$MATE_PAIR" == "true" ]]; then


        if [[ ${METHOD} == "haloplex" ]]; then
            # Get sequencing tags
            . $SERA_PATH/config/sequencingTags.sh;
            # Check that output file doesn't exist then run cutAdapt, if it does print error message
            if [[ ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz && ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz || ! -z $FORCE ]]; then

            # Run cutadapt
            singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY cutadapt -a $tTag -A `perl $SERA_PATH/bin/perlscript/reverseComplement.pl $fTag` -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz -p ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz --minimum-length 1 ${READ1} ${READ2} > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;

        else
        ErrorLog "${SAMPLEID}" "${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz and ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz already exists and force was NOT used!";
        fi

        elif [[ ${METHOD} == "swift" ]]; then
            # Check that output file doesn't exist then run cutAdapt, if it does print error message
            if [[ ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.tmp.fastq.gz && ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.tmp.fastq.gz || ! -z $FORCE ]]; then

                ILLUMINA_ADAPTER_TRIMMOMATIC="${ROOT_PATH}/refFiles/TruSeq3-PE-2.fa";
                PE1_G_T="${PREFIX}_R1_trimd.fq.gz";
                PE2_G_T="${PREFIX}_R2_trimd.fq.gz";

                singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY java -Xmx32g -Xms16g -jar ${TRIMMOMATIC_PATH}/trimmomatic.jar PE \
                    -threads 12 -trimlog $TRIM_LOG \
                    ${READ1} ${READ2} ${PE1_G_T} ${PREFIX}_unpaired_R1.fq.gz \
                    ${PE2_G_T} ${PREFIX}_unpaired_R2.fq.gz \
                    ILLUMINACLIP:${ILLUMINA_ADAPTER_TRIMMOMATIC}:2:30:10 \
                    MINLEN:30

                singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY gunzip -f ${PE1_G_T};
                singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY gunzip -f ${PE2_G_T};
                PE1_T="${PREFIX}_R1_trimd.fq";
                PE2_T="${PREFIX}_R2_trimd.fq";

                # convert fastq format to one line per record for splitting
                singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY paste - - - - < $PE1_T | tr '\t' '~' > ${PE1_T}.tmp1;
                singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY paste - - - - < $PE2_T | tr '\t' '~' > ${PE2_T}.tmp1;

                # get number of fastq records in sample before converting back to fastq format
                l=$(wc -l ${PE1_T}.tmp1 | awk '{print $1}')
                chunklinecnt=$(( $l / $cputhreads ))

                # split re-formatted fastq files into chunks
                singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY split -d -l $chunklinecnt ${PE1_T}.tmp1 ${PREFIX}_r1split
                singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY split -d -l $chunklinecnt ${PE2_T}.tmp1 ${PREFIX}_r2split

                # 20160722 NOTE: adding the sample-specific prefix to the temporary fastq
                #                files should fix the sample concatenation bug
                # convert each chunk back to fastq format
                singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY parallel r1reformat ::: ${PREFIX}_r1split*
                singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY parallel r2reformat ::: ${PREFIX}_r2split*

                ls ${PREFIX}_r1*.fastq > ${SNIC_TMP}/r1infiles
                ls ${PREFIX}_r2*.fastq > ${SNIC_TMP}/r2infiles

                # run parallel on paired chunks of fastq files with ptrim() function
                singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY parallel --xapply ptrim {1} {2} ::: $(cat ${SNIC_TMP}/r1infiles) ::: $(cat ${SNIC_TMP}/r2infiles)

                # merge trimmed files to one output file
                singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY sh -c "cat ${PREFIX}*_5ptmpR1.fq | gzip -f > ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.tmp.fastq.gz";
                singularity exec -B /data -B /opt -B /beegfs-storage -B /projects -B $SERA_PATH $SERA_SINGULARITY sh -c "cat ${PREFIX}*_5ptmpR2.fq | gzip -f > ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.tmp.fastq.gz";

                rm ${PREFIX}_r[1,2]split*
                rm ${PE1_T}.tmp1
                rm ${PE2_T}.tmp1
            else
                ErrorLog "${SAMPLEID}" "${ROOT_PATH}/seqdata/${SAMPLEID}.read1.tmp.fastq.gz and ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.tmp.fastq.gz already exists and force was NOT used!";
            fi
        else
            ErrorLog "${SAMPLEID}" "Only implemented for method haloplex and swift so far!!!";
        fi
    else
        ErrorLog "${SAMPLEID}" "Only implemented for paired-end sequencing!";
    fi
fi


if [[ "$?" != "0" ]]; then
    ErrorLog "${SAMPLEID}" "Failed in cutadapt...";
else
    SuccessLog "${SAMPLEID}" "Passed cutadapt";
fi
