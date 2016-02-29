#!/bin/bash

###########################
### SERVER-SPECIFIC OPTIONS
###########################
# please make sure that all softwares paths (Mosaik, samtools etc) are available from bash without specifying path

# paths to softwares
export ROOT_PATH_SEDD="/usr/local/sedd";
export ROOT_PATH_CORRSEL="/usr/local/corrsel";
#export ROOT_PATH_JSNPMANIA="/usr/local/jSNPmania-0.0.7-v9";
export ROOT_PATH_JSNPMANIA="/home/tom/jSNPMania/jSNPmania-0.0.7-SNAPSHOT-soft_clip";
export DOWNLOAD2FASTA="$SERA_PATH/bin/perlscript/download2fastaNoThreads.pl";
export VARIANTPREDICTION="/usr/local/bin/variant_effect_predictor.pl";
export CUTADAPT="/usr/local/cutadapt-0.9.5/cutadapt";
export ROOT_PATH_MOSAIK="/usr/local/MOSAIK-2.1.33-Linux-x64";
export ROOT_PATH_READS2MOLECULES="/usr/local/Reads2Molecules-0.0.1-SNAPSHOT";
export ROOT_PATH_BLAST="/usr/bin"
# java flags, use quick accessible temp dir for increased performance
export JAVA_FLAGS="-Xmx16384m"


###########################
### FILE LOCATIONS
###########################

# Log file.
export LOG_FILE="$ROOT_PATH/seralog.xml";

# Aligner reference genome files
export GENOME_REF="/data/hg18/Mosaik.ref.seqs/h.sapiens_mosaik_1.1/h.sapiens.ref";

# Reference genome sequence (fasta)
#export GENOME_FASTA_REF="/data/hg19/Mosaik.ref.seqs/h.sapiens_hg19_mosaik_1.1/hg19_allChromosome_correctFastaHead_20101018.corrected.fasta"; # format chr1#nc#length#-1
export GENOME_FASTA_REF="/data/hg18/genome_fasta/hg19.with.mt.fasta"; # format chr1

# Blast DB
export BLAST_DB="/data/hg18/blastdb/human_genomic";

# CorrSel reference converter file
export CORRSEL_REF_CONV_FILE="$SERA_PATH/config/reference.hg19.Info";  # format chr1

# File containing allele frequencies in the 1000 Genome Project
export ONEKGENOMEPROJECT_AF="/data/1000Genomes/alleleFrequencies_interim_phase1.frq.gz";

# BLAST analysis settings
export BLAST_MM=1;  # Allowed missmatches
export BLAST_ML=20; # Minimum match length in bp.
