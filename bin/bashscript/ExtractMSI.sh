#!/bin/bash
#
##SBATCH --qos=short 
#SBATCH -p devcore  -n 1
#SBATCH -t 01:00:00


# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts extracting MSI markers ...";

# Check if the directory exists, if not create it
if [[ ! -d $ROOT_PATH/Extracted_sampleInfo ]]; then 
	mkdir $ROOT_PATH/Extracted_sampleInfo;
fi
 


if [[ ${NORMAL_SAMPLEID} != "false" ]]; then
	if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
		if [[ ${NORMAL_SAMPLEID} == "annovar" ]]; then
			if [[ ${TISSUE} == "colon" ]]; then
				if [[ -e $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput ]]; then
					if [[ ! -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.ampliconmapped.txt  || ! -z $FORCE ]]; then
						python $SERA_PATH/bin/pythonscript/ExtractMSImarkers.py -i $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput -o $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.ampliconmapped.txt --tgfbr2Ratio_1bp 0.1 --tgfbr2Ratio_2bp 0.05 --acvr2aRatio_1bp 0.1 --acvr2aRatio_2bp 0.05;
						echo "python $SERA_PATH/bin/pythonscript/ExtractMSImarkers.py -i $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput -o $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.ampliconmapped.txt --tgfbr2Ratio_1bp 0.1 --tgfbr2Ratio_2bp 0.05 --acvr2aRatio_1bp 0.1 --acvr2aRatio_2bp 0.05";
					else 
						ErrorLog "$SAMPLEID" "$ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.ampliconmapped.txt already exists and force was not used";					
					fi
				elif [[ -e $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput ]]; then
					if [[ ! -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.txt  || ! -z $FORCE ]]; then
						python $SERA_PATH/bin/pythonscript/ExtractMSImarkers.py -i $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput -o $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.txt --tgfbr2Ratio_1bp 0.1 --tgfbr2Ratio_2bp 0.05 --acvr2aRatio_1bp 0.1 --acvr2aRatio_2bp 0.05;
						echo "python $SERA_PATH/bin/pythonscript/ExtractMSImarkers.py -i $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput -o $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.txt --tgfbr2Ratio_1bp 0.1 --tgfbr2Ratio_2bp 0.05 --acvr2aRatio_1bp 0.1 --acvr2aRatio_2bp 0.05";
					else 
						ErrorLog "$SAMPLEID" "$ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.txt already exists and force was not used";
					fi
				else 
					ErrorLog "$SAMPLEID" "Neither $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput nor t$ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput exists. Using AnnovarOutput and PindelAnnovarOutput instead!";
					if [[ -e $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput && -e $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput ]]; then
						if [[ ! -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.ampliconmapped.txt  || ! -z $FORCE ]]; then
							cat $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput | python $SERA_PATH/bin/pythonscript/ExtractMSImarkers.py -i /dev/stdin -o $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.ampliconmapped.txt --tgfbr2Ratio_1bp 0.1 --tgfbr2Ratio_2bp 0.05 --acvr2aRatio_1bp 0.1 --acvr2aRatio_2bp 0.05;
							echo "cat $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput | python $SERA_PATH/bin/pythonscript/ExtractMSImarkers.py -i /dev/stdin -o $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.ampliconmapped.txt --tgfbr2Ratio_1bp 0.1 --tgfbr2Ratio_2bp 0.05 --acvr2aRatio_1bp 0.1 --acvr2aRatio_2bp 0.05";
						else 
						ErrorLog "$SAMPLEID" "$ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.ampliconmapped.txt already exists and force was not used";					
						fi
					elif [[ -e $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput && -e $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput ]]; then
						if [[ ! -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.txt  || ! -z $FORCE ]]; then
							cat $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput | python $SERA_PATH/bin/pythonscript/ExtractMSImarkers.py -i /dev/stdin -o $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.txt --tgfbr2Ratio_1bp 0.1 --tgfbr2Ratio_2bp 0.05 --acvr2aRatio_1bp 0.1 --acvr2aRatio_2bp 0.05;
							echo "cat $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput | python $SERA_PATH/bin/pythonscript/ExtractMSImarkers.py -i /dev/stdin -o $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.txt --tgfbr2Ratio_1bp 0.1 --tgfbr2Ratio_2bp 0.05 --acvr2aRatio_1bp 0.1 --acvr2aRatio_2bp 0.05";
						else 
							ErrorLog "$SAMPLEID" "$ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.txt already exists and force was not used";
						fi
					else
						ErrorLog "$SAMPLEID" "The combination of AnnovarOutput and PindelAnnovarOutput doesn't exist independent of ampliconmapped or not ($ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput, $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput and $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput)!";
					fi
				fi
			else
				ErrorLog "$SAMPLEID" "MSI markers is only extracted for colon cancer!";
			fi
		else
			ErrorLog "$SAMPLEID" "Only run for NORMAL_SAMPLEID annovar";
		fi
				
	else
		ErrorLog "$SAMPLEID" "Only supported for call_type h.sapiens so far!";
	fi
else
	ErrorLog "$SAMPLEID" "Normal_sampleid is false -> annovar is not run!";
fi

if [[ "$?" != "0" ]]; then
	ErrorLog "$SAMPLEID" "Failed in extracting MSI markers";
else
	SuccessLog "$SAMPLEID" "Passed extracting MSI markers";
fi
