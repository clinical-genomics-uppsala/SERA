#!/bin/bash
#
# Script creates amproi and ampregion files.
#
##SBATCH -p core -n 1
##SBATCH -t 15:00
#SBATCH --qos=short -t 00:15:00

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [ ! -d "$ROOT_PATH/refFiles" ]; then
	mkdir $ROOT_PATH/refFiles;
fi

SuccessLog "${REFSEQ}" "Creating Ampregion, Seqregion, SeqROI files...";

# copy selection and ROI files
if [[ ! -s ${ROOT_PATH}/refFiles/${REFSEQ}.selection || ! -z $FORCE ]]; then
	# If selection and roi files are in bed-format
	if [ $FILE_FORMAT == "bed" ]; then
		perl $SERA_PATH/bin/perlscript/chr2nc.pl -i ${ROIFILE} -c 1 -chr2nc $NC2chr -o /dev/stdout | awk 'BEGIN{FS="\t"; c=0}{if ($1~/^NC/){c++; print $4"_ROI_"c"\t"$1"\t"$2+1"\t"$3"\t1"}}' > $ROOT_PATH/refFiles/${REFSEQ}.roi;
		perl $SERA_PATH/bin/perlscript/chr2nc.pl -i ${SELECTIONFILE} -c 1 -chr2nc $NC2chr -o /dev/stdout | awk 'BEGIN{FS="\t"}{if ($1~/^NC/){print $1"\t"$2"\t"$3"\t"$6}}' | sort -u | awk 'BEGIN{FS="\t";c=0;}{if ($1~/^NC/){c++; if($4=="+"){strand=1}else{strand=-1} print "Amplicon_"c"\t"$1"\t"$2+1"\t"$3"\t"strand}}' > $ROOT_PATH/refFiles/${REFSEQ}.selection;
	else
		cp ${SELECTIONFILE} $ROOT_PATH/refFiles/${REFSEQ}.selection;
		cp ${ROIFILE} $ROOT_PATH/refFiles/${REFSEQ}.roi;
	fi
else
	ErrorLog "${REFSEQ}" "${ROOT_PATH}/refFiles/${REFSEQ}.selection already existed and force was NOT used!";
fi

# Create ampregion file
if [ ! -e $ROOT_PATH/refFiles/${REFSEQ}.ampregion ]; then
	$ROOT_PATH_SEDD/mergeregions.sh $ROOT_PATH/refFiles/${REFSEQ}.selection | sort -k2,3V | awk '{if($1!~/^#/){print $0"\t1"} else{print $0}}' > $ROOT_PATH/refFiles/${REFSEQ}.ampregion;
	
	# Check if creating ampregion file worked
	if [ "$?" != "0" ]; then
		ErrorLog "${REFSEQ}" "Failed in create amplified region file.";
	else
		SuccessLog "${REFSEQ}" "Passed create amplified region file.";
	fi
else
	ErrorLog "${REFSEQ}" "${ROOT_PATH}/refFiles/${REFSEQ}.ampregion already existed and force was NOT used!";
fi

# Create seqregion file
if [[ ! -e $ROOT_PATH/refFiles/${REFSEQ}.seqregion || ! -z $FORCE ]]; then
	$ROOT_PATH_SEDD/extractends.sh $ROOT_PATH/refFiles/${REFSEQ}.selection $READ_LENGTH | $ROOT_PATH_SEDD/mergeregions.sh | awk '{if($1!~/^#/){print $0"\t1"} else{print $0}}' > $ROOT_PATH/refFiles/${REFSEQ}.seqregion;
	
	# Check if creating ampregion file worked
	if [ "$?" != "0" ]; then
		ErrorLog "${REFSEQ}" "Failed in create sequencible region file.";
	else
		SuccessLog "${REFSEQ}" "Passed create sequencible region file.";
	fi
else
	ErrorLog "${REFSEQ}" "${ROOT_PATH}/refFiles/${REFSEQ}.seqregion already existed and force was NOT used!";
fi

# Create seqroi file
if [[ ! -e $ROOT_PATH/refFiles/${REFSEQ}.seqroi || ! -z $FORCE ]]; then
	$ROOT_PATH_SEDD/intersection.sh $ROOT_PATH/refFiles/${REFSEQ}.seqregion $ROOT_PATH/refFiles/${REFSEQ}.roi | awk '{if($1!~/^#/){print $0"\t1"} else{print $0}}' > $ROOT_PATH/refFiles/${REFSEQ}.seqroi;
	# Check if creating amproi file worked
	if [ "$?" != "0" ]; then
		ErrorLog "${REFSEQ}" "Failed in create sequencible-roi file.";
	else
		SuccessLog "${REFSEQ}" "Passed create sequencible-roi file.";
	fi
else
	ErrorLog "${REFSEQ}" "${ROOT_PATH}/refFiles/${REFSEQ}.seqroi already existed and force was NOT used!";
fi
