#!/bin/bash -l
#SBATCH -p devcore  -n 1
#SBATCH -t 01:00:00
##SBATCH --qos=short

module load bioinfo-tools

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog "${SAMPLEID}" "Start running cutadapt";

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/seqdata" ]; then
	mkdir $ROOT_PATH/seqdata;
fi

if [ $PLATFORM = "Illumina" ]; then

	# get sequencing tags
	. $SERA_PATH/config/sequencingTags.sh;

	# if file ending not fasta/fastq
# 	zcat $RAWDATA_PE1 > $SNIC_TMP/pe1.fastq;
# 	zcat $RAWDATA_PE2 > $SNIC_TMP/pe2.fastq;
# 	RAWDATA_PE1="$SNIC_TMP/pe1.fastq";
# 	RAWDATA_PE2="$SNIC_TMP/pe2.fastq";

	# If MATE_PAIR is set to true in the input file 
	if [ "$MATE_PAIR" == "true" ]; then
		# Check that output file doesn't exist then run cutAdapt, if it does print error message
		if [[ ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz && ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz || ! -z $FORCE ]]; then
		    
		    if [[ ${TYPE} == "ffpe" ]]; then
    			cutadapt -a $tTag -A `$SERA_PATH/bin/perlscript/reverseComplement.pl $fTag` -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz -p ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz --minimum-length 1 $RAWDATA_PE1 $RAWDATA_PE2 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;

	       		SuccessLog "${SAMPLEID}" "cutadapt -a $tTag -A `$SERA_PATH/bin/perlscript/reverseComplement.pl $fTag` -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz -p ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz --minimum-length 1 $RAWDATA_PE1 $RAWDATA_PE2 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;"
            elif [[ ${TYPE} == "plasma" ]]; then
                if [ ${CUTADAPT_PREFIX} != "false" ]; then
                
                    TMP1_PE1="$SNIC_TMP/pe1.tmp1.fastq.gz";
                    TMP1_PE2="$SNIC_TMP/pe2.tmp1.fastq.gz";

                    TEMP1_PE1="$SNIC_TMP/pe1.temp1.fastq.gz";
                    TEMP1_PE2="$SNIC_TMP/pe2.temp1.fastq.gz";

                    TEMP2_PE1="$SNIC_TMP/pe1.temp2.fastq.gz";
                    TEMP2_PE2="$SNIC_TMP/pe2.temp2.fastq.gz";

                    TEMP3_PE1="$SNIC_TMP/pe1.temp3.fastq.gz";
                    TEMP3_PE2="$SNIC_TMP/pe2.temp3.fastq.gz";

                    cutadaptFile5prim="${ROOT_PATH}/refFiles/${CUTADAPT_PREFIX}_5ptrim.fa";
                    cutadaptFile3prim="${ROOT_PATH}/refFiles/${CUTADAPT_PREFIX}_3ptrim.fa";

                    #Sequences given by swift
                    #Adapter       AGATCGGAAGAGC ACACGTCTGAACTCCAGTCA
                    #AdapterRead2  AGATCGGAAGAGC GTCGTGTAGGGAAAGAGTGT
                    cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $TMP1_PE1 -p $TMP1_PE2 --minimum-length 1 $RAWDATA_PE1 $RAWDATA_PE2 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
                    SuccessLog "cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $TMP1_PE1 -p $TMP1_PE2 --minimum-length 1 $RAWDATA_PE1 $RAWDATA_PE2 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log";
                    
                    cutadapt -g file:$cutadaptFile5prim -o $TEMP1_PE1 -p $TEMP1_PE2 $TMP1_PE1 $TMP1_PE2 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
                    SuccessLog "cutadapt -g file:$cutadaptFile5prim -o $TEMP1_PE1 -p $TEMP1_PE2 $TMP_PE1 $TMP_PE2 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log";
#                    cutadapt -g file:$cutadaptFile5prim -o $TEMP1_PE1 -p $TEMP1_PE2 ${RAWDATA_PE1} ${RAWDATA_PE2} --minimum-length 40 -e 0.12 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
#                    SuccessLog "${SAMPLEID}" "cutadapt -g file:$cutadaptFile5prim -o $TEMP1_PE1 -p $TEMP1_PE2 ${RAWDATA_PE1} ${RAWDATA_PE2} --minimum-length 40 -e 0.12 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log";


                    cutadapt -g file:$cutadaptFile5prim -o $TEMP2_PE2 -p $TEMP2_PE1 $TEMP1_PE2 $TEMP1_PE1 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
                    SuccessLog "${SAMPLEID}" "cutadapt -g file:$cutadaptFile5prim -o $TEMP2_PE2 -p $TEMP2_PE2 $TEMP1_PE2 $TEMP1_PE1 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log";
                    
                    cutadapt -g file:$cutadaptFile3prim -o $TEMP3_PE1 -p $TEMP3_PE2 $TEMP2_PE1 $TEMP2_PE2 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
                    SuccessLog "${SAMPLEID}" "cutadapt -g file:$cutadaptFile3prim -o $TEMP3_PE1 -p $TEMP3_PE2 $TEMP2_PE1 $TEMP2_PE2 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log";
                    
                    cutadapt -g file:$cutadaptFile3prim -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz -p ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz $TEMP3_PE2 $TEMP3_PE1 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
                    SuccessLog "${SAMPLEID}" "cutadapt -g file:$cutadaptFile3prim -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz -p ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz $TEMP3_PE2 $TEMP3_PE1 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log";
                else
                    ErrorLog "${SAMPLEID}" "CUTADAPT_PREFIX has to be set for plasma in order to be able to run cutadapt!";
                fi
            else
                ErrorLog "${SAMPLEID}" "Only implemented for TYPE ffpe and plasma so far!";
            fi

		else 
			ErrorLog "${SAMPLEID}" "${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz and ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz already exists and force was NOT used!";
		fi
	else
		ErrorLog "${SAMPLEID}" "Only implemented for paired-end sequencing!";
	fi
fi


if [ "$?" != "0" ]; then
	ErrorLog "${SAMPLEID}" "Failed in cutadapt...";
else
	SuccessLog "${SAMPLEID}" "Passed cutadapt";
fi