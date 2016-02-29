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
title="SlurmSERA v.1.0"

# array structure:
# id brief_explanation default_status associated_script dependency_on_id run_after_samples(0: do not wait, 1: wait and run all samples, 2: wait and run once)
steps=(
0 "Create ampregion, ampROI, seqregion and seqROI Files" off "Create_region_files_python.sh" false 2 \
1 "Create Mosaik Reference Files" off "run_Step0_create_MosaikReferenceFiles.sh" 0 0 \
2 "Remove adapter sequences from reads" off "Run_cutAdapt.sh" false 0 \
3 "Build Mosaik Reads Files" off "Run_MosaikBuild.sh" 2 0 \
4 "MosaikAligner (Against Region)" off "Run_MosaikAligner_region.sh" "1:2:3" 0 \
5 "MosaikAligner (Against Genome)" off "Run_MosaikAligner_genome.sh" 3 0 \
6 "Create correct Alignment BAM-format" off "Export_BAMfiles.sh" "4:5" 0 \
7 "Run BarD" off "Run_BarD.sh" "6" 0 \
8 "Amplicon mapping" off "AmpliconMapping.sh" "6:7" 0 \
9 "Count reads and molecules per amplicon" off "ReadsPermolecules_vs_readDepth_bam.sh" "8" 0 \
10 "Run Pileup" off "Run_PileUp.sh" "6" 0 \
11 "Separate on- and offtarget bases from pileup" off "PileUp2onOffTarget.sh" "10" 0 \
12 "Extract specificities" off "ExtractSpecificity.sh" "11" 0 \
20 "Create SNPseq file" off "Create_SNPseq_file.sh" 1 2 \
21 "Run jSNPmania" off "Run_jSNPmania.sh" "20:6:8" 0 \
22 "SNPmania output to Annovar input" off "jSNPmania2AnnovarInput.sh" "21" 0 \
23 "Run Annovar" off "Run_Annovar.sh" "22" 0 \
24 "Combine Annovar output to one file" off "Combine_annovarOutput_to_one_file.sh" "23" 2 \
25 "Calculate variant allele ratio per strand" off "StrandVarAlleleRatio.sh" "21" 0 \
26 "Extract info about clinical positions" off "ExtractInfoClinicalPositions.sh" "21" 0 \
27 "Calculate fraction covered" off "FractionCovered.sh" "21" 0 \
30 "Create hits per base file from SNPmania file SeqRegion" off "SNPmania2hitsPerBase_seqregion.sh" "21" 0 \
31 "Create hits per base file from SNPmania file SeqRoi" off "SNPmania2hitsPerBase_seqroi.sh" "21" 0 \
32 "Produce plot-values SeqRegion Against SeqRegion" off "Create_plotValues_seqregion_against_seqregion.sh" 30 0 \
33 "Produce plot-values SeqRoi Against SeqRoi" off "Create_plotValues_seqroi_against_seqroi.sh" 31 0 \
34 "Create stenberg nusbaum input file" off "makeGnuplotInputFilePost.sh" "32:33" 2 \
35 "Make stenberg nusbaum plots" off "MakePlots.sh" 34 2 \
36 "Calculate Mean Sequencing Depth" off "Calculate_mean_depth.sh" "30:31" 0 \
37 "Calculate mean cov per gene" off "Calculate_Gene_Cov.sh" "21" 0 \
38 "Calculate mean cov per region (exon)" off "Calculate_Region_Cov.sh" 21 0 \
39 "Create log2 read depth plot" off "log2coveragePlot.sh" 21 0 \
40 "Run Pindel" off "Run_Pindel.sh" "6" 0 \
41 "Annotate Pindel with Annovar" off "AnnotatePindel.sh" "40" 0 \
50 "Create Reads BLAST Database" off "Make_Reads_BlastDB.sh" false 0 \
51 "Run BLAST analysis" off "Run_Blast.sh" 40 0 \
52 "Run CNV from BLAST" off "Create_cnvFromBLAST.sh" 41 0 \
53 "Compare CNV from BLAST with SNParray data" off "Create_cnvFromBLAST_SNParrayComparison.sh" 42 0 \
54 "Run CNV from BLAST per restriction reaction" off "Create_cnvFromBLAST_perRE.sh" 41 0 \
55 "Make CNV plots" off "Generate_CNV_files.sh" false 1 \
56 "Make CNV plots w. linear regression" off "generateCNVwithLinearRegression.sh" false 1 \
60 "Extract junctions" off "run_Step50_Extract_Junctions.sh" 5 0 \
61 "Gather selector properties" off "run_Step51_Gather_Selector_Properties.sh" 50 0 \
62 "Create junction plots" off "run_Step52_Create_Junction_Plots.sh" "50:51" 0 \
72 "Collect summarized data" off "Create_data_table.sh" false 2 \
75 "Run goSNP Mania" off "Run_goSNPmania.sh" 21 0 \
76 "Plot allele ratios only HapMap positions" off "Plot_allele_ratios_hapmapBases.sh" "65" 0 \
77 "Plot allele ratios all bases" off "Plot_allele_ratios_allBases.sh" "65" 0 );
let NUMBOFSTEPS=${#steps[@]}-1;

# adjust steps based on input file
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
