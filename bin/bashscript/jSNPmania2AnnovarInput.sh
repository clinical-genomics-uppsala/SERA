#!/bin/bash
#
# Script to run jSNPmania
##SBATCH --qos=short 
#SBATCH -p devcore  -n 1
#SBATCH -t 15:00

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts creating annovar input from jSNPmania ...";

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/Annovar" ]; then
	mkdir $ROOT_PATH/Annovar;
fi

# Check that READS are true
if [ ${READS} == "true" ]; then
	# Check if NORMAL_SAMPLEID is set to false - in that case we are done
	if [ ${NORMAL_SAMPLEID} == "false" ]; then
		SuccessLog "$SAMPLEID" "NORMAL_SAMPLEID is set to false - no annovar input will be created!";

	# Check if the NORMAL_SAMPLEID is set to annovar - then the sample should be run without tumor - normal comparison
	elif [ ${NORMAL_SAMPLEID} == "annovar" ]; then
		# Check if the output contain amplicon information
		if [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations && -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.insertions && -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.deletions ]]; then
			if [[ ! -e $ROOT_PATH/Annovar/${SAMPLEID}.ampliconmapped.annovarInput || ! -z $FORCE ]]; then
				perl $SERA_PATH/bin/perlscript/createAnnovarInputFromSNPmania.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -i $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.insertions -d $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.deletions -o $ROOT_PATH/Annovar/${SAMPLEID}.ampliconmapped.annovarInput -nc2chr $NC2chr -nc_col 2 -chr_col 1 ${ANNOVAR_FLAGS} >> $ROOT_PATH/Annovar/createAnnovar_output.txt
			else
				ErrorLog "$SAMPLEID" "$ROOT_PATH/Annovar/${SAMPLEID}.ampliconmapped.annovarInput already exists and force was not used!"
			fi
	
			# If there is not a file with amplicon information, check if there is one without
		elif [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.variations && -e $ROOT_PATH/SNPmania/${SAMPLEID}.insertions && -e $ROOT_PATH/SNPmania/${SAMPLEID}.deletions ]]; then
			if [[ ! -e $ROOT_PATH/Annovar/${SAMPLEID}.ampliconmapped.annovarInput || ! -z $FORCE ]]; then
				perl $SERA_PATH/bin/perlscript/createAnnovarInputFromSNPmania.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.variations -i $ROOT_PATH/SNPmania/${SAMPLEID}.insertions -d $ROOT_PATH/SNPmania/${SAMPLEID}.deletions -o $ROOT_PATH/Annovar/${SAMPLEID}.annovarInput -nc2chr $NC2chr -nc_col 2 -chr_col 1 ${ANNOVAR_FLAGS} >> $ROOT_PATH/Annovar/createAnnovar_output.txt
			else
				ErrorLog "$SAMPLEID" "$ROOT_PATH/Annovar/${SAMPLEID}.ampliconmapped.annovarInput already exists and force was not used!"
			fi
		else
			ErrorLog "$SAMPLEID" "All files needed to run the analysis is not available!";
		fi
	# Run tumor-normal analysis
	else
		if [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations && -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.insertions && -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.deletions && -e $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.ampliconmapped.variations && -e $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.ampliconmapped.insertions && -e $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.ampliconmapped.deletions ]]; then
			perl $SERA_PATH/bin/perlscript/createAnnovarInputFromSNPmania_tumorNormal.pl -tv $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -ti $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.insertions -td $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.deletions -nv $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.ampliconmapped.variations -ni $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.ampliconmapped.insertions -nd $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.ampliconmapped.deletions -o $ROOT_PATH/Annovar/${SAMPLEID}_${NORMAL_SAMPLEID}.ampliconmapped.annovarInput -nc2chr $NC2chr -nc_col 2 -chr_col 1 ${TUMOR_NORMAL_FLAGS} >> $ROOT_PATH/Annovar/createAnnovar_output.txt;
			
		elif  [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.variations && -e $ROOT_PATH/SNPmania/${SAMPLEID}.insertions && -e $ROOT_PATH/SNPmania/${SAMPLEID}.deletions && -e $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.variations && -e $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.insertions && -e $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.deletions ]]; then
			perl $SERA_PATH/bin/perlscript/createAnnovarInputFromSNPmania_tumorNormal.pl -tv $ROOT_PATH/SNPmania/${SAMPLEID}.variations -ti $ROOT_PATH/SNPmania/${SAMPLEID}.insertions -td $ROOT_PATH/SNPmania/${SAMPLEID}.deletions -nv $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.variations -ni $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.insertions -nd $ROOT_PATH/SNPmania/${NORMAL_SAMPLEID}.deletions -o $ROOT_PATH/Annovar/${SAMPLEID}_${NORMAL_SAMPLEID}.annovarInput -nc2chr $NC2chr -nc_col 2 -chr_col 1 ${TUMOR_NORMAL_FLAGS} >> $ROOT_PATH/Annovar/createAnnovar_output.txt;
		else
			ErrorLog "$SAMPLEID" "All files needed to run the analysis is not available!";
		fi
	fi
else
	ErrorLog "$SAMPLEID" "READS has to be true to run the analysis!";
fi


if [ "$?" != "0" ]; then
	ErrorLog "$SAMPLEID" "Failed in jSNPmania2annovarInput";
else
	SuccessLog "$SAMPLEID" "Passed jSNPmania2annovarInput";
fi
