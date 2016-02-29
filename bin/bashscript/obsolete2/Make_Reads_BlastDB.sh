#!/bin/bash -l
#
# Script creats BlastDB on mate-pair reads.
#
#SBATCH -p core -n 1
#SBATCH -t 1:00:00

module load bioinfo-tools blast/2.2.25

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if BLAST analysis is possible.
if [ $PLATFORM == "SOLiD" ]; then
	ErrorLog "${SAMPLEID}" "Analysis is not supported for SOLiD data.";
	exit 0;
fi
if [ $MATE_PAIR != "true" ]; then
	ErrorLog "${SAMPLEID}" "Mate-pair data needed for this analysis.";
	exit 0;
fi

# Check if directories exists, if not create.
if [ ! -d "$ROOT_PATH/BLAST" ]; then
	mkdir $ROOT_PATH/BLAST;
fi
if [ ! -d "$ROOT_PATH/BLAST/Databases" ]; then
	mkdir $ROOT_PATH/BLAST/Databases;
fi

SuccessLog "${SAMPLEID}" "Starting converting reads to BlastDB.";


# If there is not a file with the sequences for the selectors create it
if [ ! -e $ROOT_PATH/refFiles/${REFSEQ}.selection.seq ]; then
	perl $DOWNLOAD2FASTA -th 4 -t full -s $SELECTIONFILE -o $ROOT_PATH/refFiles/${REFSEQ}.selection.seq -d $BLAST_DB;
fi

if [ ${READS} == "true" ]; then
	# If reads are cut with cutAdapt
	if [[ $ROOT_PATH/filtereddata/${SAMPLEID}.read1.fastq.gz && $ROOT_PATH/filtereddata/${SAMPLEID}.read2.fastq.gz ]]; then
		perl $SERA_PATH/bin/perlscript/seraBlast.pl -m1 $ROOT_PATH/filtereddata/${SAMPLEID}.read1.fastq.gz  -m2 $ROOT_PATH/filtereddata/${SAMPLEID}.read2.fastq.gz -rl 10 -s $ROOT_PATH/refFiles/${REFSEQ}.selection.seq -p $ROOT_PATH/BLAST/Databases/${SAMPLEID};
	elif [[ $ROOT_PATH/seqdata/${SAMPLEID}.read1.fastq.gz && $ROOT_PATH/seqdata/${SAMPLEID}.read2.fastq.gz ]]; then
		 perl $SERA_PATH/bin/perlscript/seraBlast.pl -m1 $ROOT_PATH/seqdata/${SAMPLEID}.read1.fastq.gz  -m2 $ROOT_PATH/seqdata/${SAMPLEID}.read2.fastq.gz -rl 10 -s $ROOT_PATH/refFiles/${REFSEQ}.selection.seq -p $ROOT_PATH/BLAST/Databases/${SAMPLEID};
	elif [[ ${RAWDATA_PE1} && ${RAWDATA_PE2} ]]; then 
		perl $SERA_PATH/bin/perlscript/seraBlast.pl -m1 ${RAWDATA_PE1} -m2 ${RAWDATA_PE2} -rl 10 -s $ROOT_PATH/refFiles/${REFSEQ}.selection.seq -p $ROOT_PATH/BLAST/Databases/${SAMPLEID};
	else
		ErrorLog "${SAMPLEID}" "No reads available!";
	fi
fi

# Creating reads BlastDB
cd $ROOT_PATH/BLAST/Databases;

formatdb -i $ROOT_PATH/BLAST/Databases/${SAMPLEID}.reads.fasta -l $ROOT_PATH/BLAST/Databases/${SAMPLEID}.createDB.log -n ${SAMPLEID}.ReadsDB -p F;

if [ "$?" != "0" ]; then
	ErrorLog "${SAMPLEID}" "Failed to create BlastDB...";
else
	gzip -f $ROOT_PATH/BLAST/Databases/${SAMPLEID}.fragments.fasta;

	SuccessLog "${SAMPLEID}" "Reads converted to BlastDB.";
fi

