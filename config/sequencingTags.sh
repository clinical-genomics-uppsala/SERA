#!/bin/bash
#
# Include file for sequencing tags Illumina.
#
# General motifs
general_IlluminaFP="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
general_IlluminaTP1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
general_IlluminaTP2="ATCTCGTATGCCGTCTTCTGCTTG";

# Create Illumina tags
fTag=$general_IlluminaFP;
tTag="${general_IlluminaTP1}";
#tTag="${general_IlluminaTP1}${SEQUNCING_TAG}${general_IlluminaTP2}";
