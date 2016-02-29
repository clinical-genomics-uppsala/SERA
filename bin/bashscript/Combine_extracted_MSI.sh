#!/bin/bash
#
# Script to run jSNPmania
#SBATCH -p devcore -n 1
#SBATCH -t 15:00

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts combining Annovar output ...";

# Check if the directory exists, if not create it
if [ ! -d $ROOT_PATH/MSIanalysis ]; then 
	mkdir $ROOT_PATH/MSIanalysis;
fi 

if [ ${READS} == "true" ]; then
	if [ ${CALL_TYPE} == "h.sapiens" ]; then

		# MSI markers - ampliconmapped
		# Check if a combined file exist, then remove it
		if [ -e $ROOT_PATH/MSIanalysis/MSImarkers.allSamples.ampliconmapped.txt ]; then
			rm $ROOT_PATH/MSIanalysis/MSImarkers.allSamples.ampliconmapped.txt;
		fi
		# If there still exist input files - create combined file
		ls -1 $ROOT_PATH/Extracted_sampleInfo/*.msiMarkers.ampliconmapped.txt > /dev/null 2>&1
		if [  "$?" = "0" ]; then
			awk 'BEGIN{print "#Sample\tChr\tStart\tEnd\tReference_base\tVariant_base\tGene\tType\tExonic_type\tVariant_allele_ratio\t#reference_alleles\t#_variant_alleles\tRead_depth\tRatio_in_1000Genome\tdbSNP_id\tClinically_flagged_dbSNP\tESP_6500\tCosmic\tClinVar_CLNDBN\tClinVar_CLINSIG Strands_A\tStrands_G\tStrands_C\tStrands_T\tStrands_Ins\tStrands_Del\t#variant_+_amplicons\t#variant_-_amplicons\t#reference_+_amplicons\t#reference_-_amplicons\tVariant_ampliconinfo\tReference_ampliconinfo\tTranscripts"} {if($1!~/#Run/){print $0}}' $ROOT_PATH/Extracted_sampleInfo/*.msiMarkers.ampliconmapped.txt > $ROOT_PATH/MSIanalysis/MSImarkers.allSamples.ampliconmapped.txt;
			
		else
			ErrorLog "$SAMPLEID" "No files with ending msiMarkers.ampliconmapped.txt in Extracted_sampleInfo were found!";
		fi

		# MSI markers - without amplicon mapping
		# Check if a combined file exist, then remove it
		if [ -e $ROOT_PATH/MSIanalysis/MSImarkers.allSamples.txt ]; then
			rm $ROOT_PATH/MSIanalysis/MSImarkers.allSamples.txt;
		fi
		# If there still exist input files - create combined file
		ls -1 $ROOT_PATH/Extracted_sampleInfo/*.msiMarkers.txt > /dev/null 2>&1
		if [  "$?" = "0" ]; then
			awk 'BEGIN{print "#Sample\tChr\tStart\tEnd\tReference_base\tVariant_base\tGene\tType\tExonic_type\tVariant_allele_ratio\t#reference_alleles\t#_variant_alleles\tRead_depth\tRatio_in_1000Genome\tdbSNP_id\tClinically_flagged_dbSNP\tESP_6500\tCosmic\tClinVar_CLNDBN\tClinVar_CLINSIG Strands_A\tStrands_G\tStrands_C\tStrands_T\tStrands_Ins\tStrands_Del\t#variant_+_amplicons\t#variant_-_amplicons\t#reference_+_amplicons\t#reference_-_amplicons\tVariant_ampliconinfo\tReference_ampliconinfo\tTranscripts"} {if($1!~/#Run/){print $0}}' $ROOT_PATH/Extracted_sampleInfo/*.msiMarkers.txt > $ROOT_PATH/MSIanalysis/MSImarkers.allSamples.txt;
			
		else
			ErrorLog "$SAMPLEID" "No files with ending msiMarkers.txt in Extracted_sampleInfo were found!";
		fi


	else
		ErrorLog "$SAMPLEID" "Only supported for h.sapiens so far!";
	fi
else
	ErrorLog "$SAMPLEID" "READS has to be true for the analysis to run!";
fi

if [ "$?" != "0" ]; then
	ErrorLog "$SAMPLEID" "Failed in combining EGFR extractions";
else
	SuccessLog "$SAMPLEID" "Passed combining EGFR extractions";
fi
