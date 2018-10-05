#!/bin/bash
#
# Script for running MosaikBuild on either SOLiD or Illumina samples.
#
#SBATCH -p core  -n 1
#SBATCH -t 1:00:00
##SBATCH --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

# Include functions
. $SERA_PATH/includes/logging.sh

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/FastQC" ]]; then
    mkdir $ROOT_PATH/FastQC;
fi

if [[ "$MATE_PAIR" == "true" ]]; then
    extraLoginfo="mate-pair/pair-ends";
else
    extraLoginfo="single";
fi

SuccessLog "${SAMPLEID}" "Building $extraLoginfo reads Mosaik input file...";

# Run FastQC
# If platform is Illumina
if [[ $PLATFORM = "Illumina" ]]; then

    # Check that the output file doesn't exist or if force is given
    # Do paired end library or not
    if [[ "$MATE_PAIR" == "true" ]]; then
        # Check that input files exist, if not print error message
        if [[ -e ${RAWDATA_PE1} && ${RAWDATA_PE2} ]]; then
            fastqc -o FastQC ${RAWDATA_PE1} ${RAWDATA_PE2} > $ROOT_PATH/FastQC/fastqc_${SAMPLEID}_output.txt;
        else
            ErrorLog "${SAMPLEID}" "${RAWDATA_PE1} and/or ${RAWDATA_PE2} do NOT exist!";
        fi

        if [[ ! -e $ROOT_PATH/FastQC/${SAMPLEID}.read1_fastqc.html || ! -e $ROOT_PATH/FastQC/${SAMPLEID}.read2_fastqc.html || ! -z $FORCE ]]; then
            # Check if cutadapted files exist - if so run cutadapt on them
            if [[ -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz && -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz ]]; then
                fastqc -o FastQC ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz > $ROOT_PATH/FastQC/fastqc_${SAMPLEID}_outputCutadapt.txt;
            else
                ErrorLog "${SAMPLEID}" "No cutadapted files exists!"
            fi
        else
            ErrorLog "${SAMPLEID}" "FastQC files already existed for cutadapted output and force was NOT used!"; 
        fi

    # If single reads are used only care about the first read
    else
        # Check that input file exists, if not print error message
        if [[ -e ${RAWDATA_PE1} ]]; then
            fastqc -o FastQC ${RAWDATA_PE1} > $ROOT_PATH/FastQC/fastqc_${SAMPLEID}_output.txt;
        else 
            ErrorLog "${SAMPLEID}" "${RAWDATA_PE1} does NOT exist!";
        fi
        
        # Check if cutadapted files exist - if so run cutadapt on them
        if [[ -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz ]]; then
            fastqc -o FastQC ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz > $ROOT_PATH/FastQC/fastqc_${SAMPLEID}_outputCutadapt.txt;
        else
            ErrorLog "${SAMPLEID}" "No cutadapted file exists!"
        fi
    fi
else
    ErrorLog "${SAMPLEID}" "Platform $PLATFORM not recognized (should be Illumina).";
    exit 1;
fi
    
# Check if MosaikBuild worked
if [[ "$?" != "0" ]]; then
    ErrorLog "${SAMPLEID}" "Failed in FastQC...";
else
    SuccessLog "${SAMPLEID}" "Passed FastQC ($extraLoginfo reads)";
fi
