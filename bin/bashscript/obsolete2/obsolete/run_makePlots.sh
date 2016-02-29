#!/bin/bash
#
# Script creates plot data files seqroi vs seqroi.
#
#SBATCH -p node -n 1
#SBATCH -t 01:00:00

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts to create PDF plots ...";

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/plotPdfs" ]; then
	mkdir $ROOT_PATH/plotPdfs;
fi

# Making and plotting data.
perl $SERA_PATH/bin/perlscript/makeGnuplotFile.pl -i plotFile -o $ROOT_PATH/plotPdfs;

# Check if makeGnuplotFile.pl worked
if [ "$?" != "0" ]; then
	ErrorLog $SAMPLEID "Failed in makeGnuplotFile.pl";
else
	SuccessLog $SAMPLEID "Success on makeGnuplotFile.pl";
fi

# Run Gnuplot on all scripts
cd $ROOT_PATH/plotPdfs;
find . -name "*.gnuplot" -exec gnuplot '{}' \;

#cd `dirname $0`;
