#!/bin/bash

# Read input options
ROOT_PATH_INPUT=""; INPUT_FILE=""; INTERFACE="" CNV_FILE=""; FORCE="";
while getopts "h:p:i:u:f" OPTION; do
	case $OPTION in
		h) usage; exit 1 ;;
		p) ROOT_PATH_INPUT=$OPTARG ;;
		i) INPUT_FILE=$OPTARG ;;
		u) INTERFACE="gui" ;;
#		c) CNV_FILE=$OPTARG ;;
		f) export FORCE="true" ;;
		?) usage; exit ;;
     esac
done

# validate input
. $SERA_PATH/includes/validateInput.sh

# Go to output folder
cd $ROOT_PATH_INPUT;
export ROOT_PATH=$(pwd);

# Include and validate  inputfile data
. $SERA_PATH/includes/validateInputFile.sh; 
COUNT=-1;

# Read input file.
. $INPUT_FILE

# Get number of samples
count=${#SAMPLEID_ARR_[@]};
let NUMBEROFSAMPLES=count-1;

# Include home or proj specific globals
if [ $GLOBALS = "HOME" ]; then 
	. $SERA_PATH/config/globalsHome.sh;
else
	. $SERA_PATH/config/globalsProj.sh;
fi
## Include server-specific globals
#if [ "$SOFTWARE" = "SLURM" ]; then
#	. $SERA_PATH/config/globalsUppnex.sh;
#else
#	. $SERA_PATH/config/globalsSmaug.sh;
#fi

# create files for logging
if [[ ! -d "$ROOT_PATH/slurmOutput" && $SOFTWARE = "SLURM" ]]; then
	mkdir slurmOutput;
fi
if [ ! -e "$LOG_FILE" ]; then
	echo -e "<?xml version='1.0'?>\n<document info=\"$title\">\n</document>" > $LOG_FILE;
fi

# inputfile information 
let ADJ_SAMPLENUMBER=NUMBEROFSAMPLES+1;
MP="single end"; if [[ $MATE_PAIR = "true" ]]; then MP="paired end"; fi;
if [[ $SOFTWARE = "SLURM" ]]; then INFO="Running on $SOFTWARE with project $UPPNEX_PROJECT_ID. ";
else INFO="Running in pure bash environment. "; fi
INFO=${INFO}"Analyzing $ADJ_SAMPLENUMBER samples from $DESIGN_TYPE design with $READ_LENGTH bp $PLATFORM $MP sequencing data. ";
if [ ! -z $FORCE ]; then INFO=${INFO}"All previous data will be overwritten for associated steps. "; fi
if [[ $ALIGN_NONUNIQUE == "true" ]]; then INFO=${INFO}"Analyzing both unique and nonunique. "; else INFO=${INFO}"Analyzing only unique reads. "; fi
if [ ! -z $CALL_TYPE ]; then INFO=${INFO}"SNP calling on $CALL_TYPE."; fi

#prepare variable array for concatenation for running post submit scripts
POST="";
