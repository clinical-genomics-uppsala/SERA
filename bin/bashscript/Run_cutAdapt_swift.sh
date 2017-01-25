#!/bin/bash -l
#SBATCH -p devcore  -n 6
#SBATCH -t 01:00:00
##SBATCH --qos=short

module load bioinfo-tools trimmomatic/0.32 gnuparallel/20150522

# Include functions
. $SERA_PATH/includes/logging.sh;

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog "${SAMPLEID}" "Start running cutadapt";

export TRIM_LOG="${ROOT_PATH}/seqdata/${SAMPLEID}.trim.log"
PREFIX="${SNIC_TMP}/${SAMPLEID}";

cputhreads=12;

ptrim()
{
    fqt1=$1
    fqt2=$2
    tprefix=${fqt1%%_R1_001.fastq}
    #trim5pfa=/path/to/your/5prime_primertrim.fasta
    #trim3pfa=/path/to/your/3prime_primertrim.fasta

    trim5pfa="/proj/b2013222/private/nobackup/20170103_AccelAmplicon_BRCA_ctDNA_all/refFiles/BRCA_5prime.fa";
    trim3pfa="/proj/b2013222/private/nobackup/20170103_AccelAmplicon_BRCA_ctDNA_all/refFiles/BRCA_3prime.fa";

    #5’trim
    cutadapt \
        -g file:$trim5pfa \
        -o ${tprefix}_tmpR1.fq -p ${tprefix}_tmpR2.fq \
    $fqt1 $fqt2 --minimum-length 40 -e 0.12 >> $TRIM_LOG;

    #5’trim
    cutadapt \
        -g file:$trim5pfa \
        -o ${tprefix}_5ptmpR2.fq -p ${tprefix}_5ptmpR1.fq \
        ${tprefix}_tmpR2.fq ${tprefix}_tmpR1.fq --minimum-length 40 -e 0.12 >> $TRIM_LOG;

    #3' trim 
    cutadapt \
        -a file:$trim3pfa \
        -o ${tprefix}_tmp3R1.fq -p ${tprefix}_tmp3R2.fq \
        ${tprefix}_5ptmpR1.fq ${tprefix}_5ptmpR2.fq --minimum-length 40 -e 0.12 >> $TRIM_LOG;

<<'comm'
    #3’trim
    cutadapt \
        -a file:$trim3pfa \
        -o ${tprefix}_R2_primertrimd.fq -p ${tprefix}_R1_primertrimd.fq \
        ${tprefix}_tmp3R2.fq ${tprefix}_tmp3R1.fq --minimum-length 40 -e 0.12 >> $TRIM_LOG;
comm
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

if [[ $PLATFORM = "Illumina" ]]; then

    # get sequencing tags
    . $SERA_PATH/config/sequencingTags.sh;
    ILLUMINA_ADAPTER_TRIMMOMATIC="/proj/b2013222/private/nobackup/20170103_AccelAmplicon_BRCA_ctDNA_all/refFiles/TruSeq3-PE-2.fa";

    # If MATE_PAIR is set to true in the input file
    if [[ "$MATE_PAIR" == "true" ]]; then
        # Check that output file doesn't exist then run cutAdapt, if it does print error message
        if [[ ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.tmp.fastq.gz && ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.tmp.fastq.gz || ! -z $FORCE ]]; then
            PE1_G_T="${PREFIX}_R1_trimd.fq.gz";
            PE2_G_T="${PREFIX}_R2_trimd.fq.gz";

            ## NOTE: custom adapter file for Accel-amplicon Illumina adapter trimming
            ##       (included as attachment in email containing this script)
#            java -Xmx24g -Xms16g -jar ${ROOT_PATH_TRIMMOMATIC}/trimmomatic-0.35.jar PE \
#                -threads 12 -trimlog $TRIM_LOG \
#                $RAWDATA_PE1 $RAWDATA_PE2 ${PE1_G_T} ${PREFIX}_unpaired_R1.fq.gz \
#                ${PE2_G_T} ${PREFIX}_unpaired_R2.fq.gz \
#                ILLUMINACLIP:${ILLUMINA_ADAPTER_TRIMMOMATIC}:2:30:10 \
#                MINLEN:30

             java -Xmx24g -Xms16g -jar $TRIMMOMATIC_HOME/trimmomatic.jar PE \
                -threads 12 -trimlog $TRIM_LOG \
                $RAWDATA_PE1 $RAWDATA_PE2 ${PE1_G_T} ${PREFIX}_unpaired_R1.fq.gz \
                ${PE2_G_T} ${PREFIX}_unpaired_R2.fq.gz \
                ILLUMINACLIP:${ILLUMINA_ADAPTER_TRIMMOMATIC}:2:30:10 \
                MINLEN:30
                rm ${PREFIX}*unpaired*.fq.gz

                gunzip -f ${PE1_G_T};
                gunzip -f ${PE2_G_T};
                PE1_T="${PREFIX}_R1_trimd.fq";
                PE2_T="${PREFIX}_R2_trimd.fq";

                # convert fastq format to one line per record for splitting
                paste - - - - < $PE1_T | tr '\t' '~' > ${PE1_T}.tmp1;
                paste - - - - < $PE2_T | tr '\t' '~' > ${PE2_T}.tmp1;

                # get number of fastq records in sample before converting back to fastq format
#                wc -l ${PE1_T}.tmp1 | tee /dev/tty | awk '{print $1}' > fqcnt
                l=$(wc -l ${PE1_T}.tmp1 | awk '{print $1}')
                chunklinecnt=$(( $l / $cputhreads ))

                # split re-formatted fastq files into chunks
                split -d -l $chunklinecnt ${PE1_T}.tmp1 ${PREFIX}_r1split
                split -d -l $chunklinecnt ${PE2_T}.tmp1 ${PREFIX}_r2split

                # 20160722 NOTE: adding the sample-specific prefix to the temporary fastq
                #                files should fix the sample concatenation bug
                # convert each chunk back to fastq format
                parallel r1reformat ::: ${PREFIX}_r1split*
                parallel r2reformat ::: ${PREFIX}_r2split*

                ls ${PREFIX}_r1*.fastq > ${SNIC_TMP}/r1infiles
                ls ${PREFIX}_r2*.fastq > ${SNIC_TMP}/r2infiles

                # run parallel on paired chunks of fastq files with ptrim() function
                parallel --xapply ptrim {1} {2} ::: $(cat ${SNIC_TMP}/r1infiles) ::: $(cat ${SNIC_TMP}/r2infiles)

                # concatenate primer-trimmed fastq chunks
#                cat ${PREFIX}*_R1_primertrimd.fq | gzip -f > ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.tmp.fastq.gz;
#                cat ${PREFIX}*_R2_primertrimd.fq | gzip -f > ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.tmp.fastq.gz

                cat ${PREFIX}*_tmp3R1.fq | gzip -f > ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.tmp.fastq.gz;
                cat ${PREFIX}*_tmp3R2.fq | gzip -f > ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.tmp.fastq.gz

                rm ${PREFIX}_r[1,2]split*
                #rm ${PREFIX}*_R[1,2]_primertrimd.fq
        else
            ErrorLog "${SAMPLEID}" "${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz and ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz already exists and force was NOT used!";
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
                                                                                                                                       