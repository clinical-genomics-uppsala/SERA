#!/bin/bash
#
# Script to convert Annovar output to VCF
#SBATCH -p devcore  -n 1
#SBATCH -t 15:00

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts converting Annovar output to vcf...";

# Check if the directory exists, if not create it
if [[ ! -d $ROOT_PATH/vcfOutput ]]; then 
	mkdir $ROOT_PATH/vcfOutput;
fi 

if [[ ${READS} == "true" ]]; then
	# Check for singleSamples without ampliconmapping
	# Check if a combined file exist, then remove it
	if [[ -e $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput ]]; then
		if [[ (! -e $ROOT_PATH/vcfOutput/${SAMPLEID}.vcf) || ($FORCE == "true") ]]; then
			python ${SERA_PATH}/bin/pythonscript/Annovar2vcf.py -v $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput -s $ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq -chr2nc $NC2chr -o $ROOT_PATH/vcfOutput/${SAMPLEID}.vcf

		else 
			ErrorLog "$SAMPLEID" "$ROOT_PATH/vcfOutput/${SAMPLEID}.vcf already exists and Force was not used!";
			
		fi

	elif [[ -e $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput ]]; then
		if [[ (! -e $ROOT_PATH/vcfOutput/${SAMPLEID}.vcf) || ($FORCE == "true") ]]; then
#				if [[ -e $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput ]]; then
#					python ${SERA_PATH}/bin/pythonscript/Annovar2vcf.py -v $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput -s $ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq -chr2nc $NC2chr -p $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput -o $ROOT_PATH/vcfOutput/${SAMPLEID}.vcf
#				else
				
			python ${SERA_PATH}/bin/pythonscript/Annovar2vcf.py -v $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput -s $ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq -chr2nc $NC2chr -o $ROOT_PATH/vcfOutput/${SAMPLEID}.vcf
#				fi
		else 
			ErrorLog "$SAMPLEID" "$ROOT_PATH/vcfOutput/${SAMPLEID}.vcf already exists and Force was not used!";
		fi

	elif [[ -e $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput ]]; then
		if [[ (! -e $ROOT_PATH/vcfOutput/${SAMPLEID}.vcf) || ($FORCE == "true") ]]; then
			python ${SERA_PATH}/bin/pythonscript/Annovar2vcf.py -v $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput -s $ROOT_PATH/Files/${REFSEQ}.ampregion.SNPseq -chr2nc $NC2chr -o $ROOT_PATH/vcfOutput/${SAMPLEID}.vcf
		else
			ErrorLog "$SAMPLEID" "$ROOT_PATH/vcfOutput/${SAMPLEID}.vcf already exists and Force was not used!";
		fi

	elif [[ -e $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput ]]; then
		if [[ (! -e $ROOT_PATH/vcfOutput/${SAMPLEID}.vcf) || ($FORCE == "true") ]]; then
#				if [[ -e $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput ]]; then
#					python ${SERA_PATH}/bin/pythonscript/Annovar2vcf.py -v $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput -s $ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq -chr2nc $NC2chr -p $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput -o $ROOT_PATH/vcfOutput/${SAMPLEID}.vcf
#				else
				
				python ${SERA_PATH}/bin/pythonscript/Annovar2vcf.py -v $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput -s $ROOT_PATH/Files/${REFSEQ}.ampregion.SNPseq -chr2nc $NC2chr -o $ROOT_PATH/vcfOutput/${SAMPLEID}.vcf
		else
			ErrorLog "$SAMPLEID" "$ROOT_PATH/vcfOutput/${SAMPLEID}.vcf already exists and Force was not used!";
		fi
	else
		ErrorLog "$SAMPLEID" "No output file from Annovar was found!";
	fi
else
	ErrorLog "$SAMPLEID" "READS has to be true for the analysis to run!";
	
fi

if [[ "$?" != "0" ]]; then
	ErrorLog "$SAMPLEID" "Failed in converting Annovar output to vcf";
else
	SuccessLog "$SAMPLEID" "Passed converting Annovar output to vcf";
fi
