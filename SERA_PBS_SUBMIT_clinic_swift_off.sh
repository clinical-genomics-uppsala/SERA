#!/bin/bash -l

# SERA Analysis Software


##############
### PREPARE ENVIRONMENT
##############

# Get script path
export SERA_PATH=$(dirname `readlink -f $0`)
. $SERA_PATH/includes/prepareEnvironment.sh

##############
### PROGRAM GLOBALS
##############

# program title
title="SlurmSERA"

# array structure:
# id brief_explanation default_status associated_script dependency_on_id run_after_samples(0: do not wait, 1: wait and run all samples, 2: wait and run once)
steps=(
1 "Create ampregion, ampROI, seqregion and seqROI Files" off "Create_region_files_python.sh" false 0 \
2 "Remove adapter sequences from reads" off "Run_cutAdapt_swift.sh" false 0 \
3 "Continue remove adapter sequences from reads" off "Run_cutAdapt_swift_part2.sh" "2" 0 \
4 "Run FastQC" off "Run_FastQC.sh" "2" 0 \
5 "Align with BWA against Genome" off "Run_Bwa.sh" "2" 0 \
6 "Amplicon mapping" off "AmpliconMapping.sh" "4" 0 \
10 "Create SNPseq file" off "Create_SNPseq_file.sh" "1" 0 \
11 "Run jSNPmania" off "Run_jSNPmania.sh" "4:5:10" 0 \
12 "SNPmania output to Annovar input" off "jSNPmania2AnnovarInput.sh" "11" 0 \
13 "Run Annovar" off "Run_Annovar.sh" "12" 0 \
14 "Combine Annovar output to one file" off "Combine_annovarOutput_to_one_file.sh" "13" 2 \
15 "Extract info about clinical positions" off "ExtractInfoClinicalPositions.sh" "11" 0 \
20 "Run Pindel" off "Run_Pindel.sh" "4" 0 \
21 "Annotate Pindel with Annovar" off "AnnotatePindel.sh" "20" 0 \
25 "Filter Annovar and Pindelannovar output" off "FilterAnnovar.sh" "13:21" 0 \
26 "Merge info about clinical positions and indels from Pindel" off "MergeAllClinicalInfo.sh" "15:21:25" 0 \
27 "Convert Annovar output to vcf-format" off "ConvertAnnovaroutput2vcf.sh" "25" 0 \
30 "Extract MSI markers" off "ExtractMSI.sh" "13:21:25" 0 \
31 "Combine extracted MSI markers to one file" off "Combine_extracted_MSI.sh" "30" 2 \
32 "Extract EGFR information" off "Extract_T790M.sh" "11" 0 \
33 "Combine extracted EGFR information to one file" off "Combine_extracted_EGFR.sh" "32" 2 \
34 "Extract SNP information" off "ExtractSNPinfo.sh" "11" 0 \
35 "Combine extracted SNP information to one file" off "Combine_extracted_SNPs.sh" "34" 2 \
40 "Calculate amplification ratio" off "CalculateAmplificationRatio.sh" "11" 0 \
50 "Calculate fraction covered" off "FractionCovered.sh" "11" 0 \
51 "Calculate mean cov per gene" off "Calculate_Gene_Cov.sh" "11" 0 \
52 "Calculate mean cov per region (exon)" off "Calculate_Region_Cov.sh" "11" 0 \
60 "Create hits per base file from SNPmania file SeqRegion" off "SNPmania2hitsPerBase_seqregion.sh" "11" 0 \
61 "Create hits per base file from SNPmania file SeqRoi" off "SNPmania2hitsPerBase_seqroi.sh" "11" 0 \
62 "Produce plot-values SeqRegion Against SeqRegion" off "Create_plotValues_seqregion_against_seqregion.sh" "60" 0 \
63 "Produce plot-values SeqRoi Against SeqRoi" off "Create_plotValues_seqroi_against_seqroi.sh" "61" 0 \
64 "Create stenberg nusbaum input file" off "makeGnuplotInputFilePost.sh" "62:63" 2 \
65 "Make stenberg nusbaum plots" off "MakePlots.sh" "64" 2 \
66 "Calculate Mean Sequencing Depth" off "Calculate_mean_depth.sh" "60:61" 0 \
70 "Run Pileup" off "Run_PileUp.sh" "4" 0 \
71 "Separate on- and offtarget bases from pileup" off "PileUp2onOffTarget.sh" "70" 0 \
72 "Extract specificities" off "ExtractSpecificity.sh" "71" 0 \
73 "Calculate variant allele ratio per strand" off "StrandVarAlleleRatio.sh" "11" 0 );
let NUMBOFSTEPS=${#steps[@]}-1;

# adjust steps based off input file
. $SERA_PATH/includes/adjustSteps.sh



##############
### START PROGRAM
##############

# Start interface
if [ -z "$INTERFACE" ]; then
    . $SERA_PATH/includes/commandLineInterface.sh
else
    . $SERA_PATH/includes/dialogInterface.sh
fi

# Run steps for samples
. $SERA_PATH/includes/submitScripts.sh

#submit successfull exit status
#exit 0;
