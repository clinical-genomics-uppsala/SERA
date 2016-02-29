#!/bin/bash

###########################
### SERVER-SPECIFIC OPTIONS
###########################
# please make sure that all softwares paths (Mosaik, samtools etc) are available from bash without specifying path

# load script modules if necessary
. $SERA_PATH/config/modules.sh

# paths to softwares
export ROOT_PATH_SEDD="$HOME/program/sedd-0.0.5";
export ROOT_PATH_JSNPMANIA="$HOME/program/jSNPmania-0.0.7-SNAPSHOT-amplicons_v4";
export ROOT_PATH_MOSAIK="$HOME/program/MOSAIK-2.1.33-Linux-x64";
export DOWNLOAD2FASTA="$SERA_PATH/bin/perlscript/download2fasta.pl";
export CUTADAPT="$HOME/program/cutadapt-0.9.5/cutadapt";
export ROOT_PATH_BARD="$HOME/program/Bard_20121008";
export ROOT_PATH_ANNOVAR="$HOME/program/annovar_2013.08.23";
export ANNOVAR_MODULES="$HOME/modules";
export ROOT_PATH_PINDEL="$HOME/program/pindel";

# java flags, use quick accessible temp dir for increased performance
export JAVA_FLAGS="-Djava.io.tmpdir=$TMPDIR -Xmx16384m"

#Since you have allocated a whole node you have the whole scratch area on the node to yourself.
#Therefore, use the scratch area for the tmp directory:
export MOSAIK_TMP=$TMPDIR;	# set to empty if not desirable

# Include SLURM submit information
SLURM_MISC_SETTINGS="-A $UPPNEX_PROJECT_ID"; 	#--qos=short";


###########################
### FILE LOCATIONS
###########################
# Blacklist file
export BLACKLIST_FILE="$SERA_PATH/res/blacklist_above300_20150902.txt";

# Log file.
export LOG_FILE="$ROOT_PATH/seralog.xml";

# Genome build hg18 or hg19
export GENOME_BUILD="hg19";

# Aligner reference genome files
export GENOME_REF="$HOME/glob/Mosaik.ref.seqs/h.sapiens_hg19_mosaik_1.1/h.sapiens.ref";

# Reference genome sequence (fasta)
export GENOME_FASTA_REF="$HOME/data/hg19/genome_fasta/hg19.with.mt.fasta"; # Aligner reference genome files

# Blast DB
export BLAST_DB="$HOME/data/hg19/blastdb/human_genomic";	# local hg 19

#NC to chr conversion
export NC2chr="$HOME/data/hg19/hg19_chr2NC.txt";

# File containing allele frequencies in the 1000 Genome Project
export ONEKGENOMEPROJECT_AF="$HOME/data/1000Genomes/alleleFrequencies_interim_phase1.frq.gz";

# BLAST analysis settings
export BLAST_MM=1;  # Allowed missmatches
export BLAST_ML=20; # Minimum match length in bp.
