#!/bin/bash
#
# Script to run jSNPmania
##SBATCH --qos=short 
#SBATCH -p devcore  -n 1
#SBATCH -t 01:00:00


# Include functions
#. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts copying files to storage ...";


# Check if the directory exists, if not create it
if [[ -d ${STORAGE_PATH} ]]; then
    
    # If STORAGE_PATH ends with a / remove the /
    if [[ ${STORAGE_PATH} == */ ]]; then
        STORAGE_PATH=${STORAGE_PATH:0:-1}
    fi
    # If ROOT_PATH ends with a / remove the /
    if [[ ${ROOT_PATH} == */ ]]; then
        ROOT_PATH=${ROOT_PATH:0:-1}
    fi

    # Copy AmpliconMapped
    if [[ -d "${ROOT_PATH}/AmpliconMapped" ]]; then
        rsync -carvl --progress "${ROOT_PATH}/AmpliconMapped" ${STORAGE_PATH}
    else
        ErrorLog "$SAMPLEID" "${ROOT_PATH}/AmpliconMapped does not exist!";
    fi

   # Go through all bwa files in the Bwa folder
    for bwaFile in "${ROOT_PATH}/Bwa/*"
    do
        parts=$(echo $bwaFile | tr "/" "\n");
        sampleFile=${parts[-1]}; # Extract filename from path
        fileParts=$(echo $sampleFile | tr "." "\n");
        sample=${fileParts[1]}; # Extract sample name
        if [[ ! -e "${ROOT_PATH}/AmpliconMapped/${sample}.*" ]]; then # Check if the sample is not ampliconmapped
            if [[ ! -d "${OUTBOX_PATH}/Bwa" ]]; then # If not check if Bwa folder exists in OUTBOX_PATH
                mkdir -m 664 "${OUTBOX_PATH}/Bwa"; # If not create
            fi
            rsync -carvl --progress $bwaFile "${OUTBOX_PATH}/Bwa"; # Copy the bwa file for the non-ampliconmapped sample
        fi
    done

    # Copy FilteredMutations
    if [[ -d "${ROOT_PATH}/FilteredMutations" ]]; then
        rsync -carvl --progress "${ROOT_PATH}/FilteredMutations" ${STORAGE_PATH}
    else
        ErrorLog "$SAMPLEID" "${ROOT_PATH}/FilteredMutations does not exist!";
    fi

    # Copy FilteredAnnovarOutput
    if [[ -d "${ROOT_PATH}/FilteredAnnovarOutput" ]]; then
        rsync -carvl --progress "${ROOT_PATH}/FilteredAnnovarOutput" ${STORAGE_PATH}
    else
        ErrorLog "$SAMPLEID" "${ROOT_PATH}/FilteredAnnovarOutput does not exist!";
    fi

    # Copy ClinicalPositions
    if [[ -d "${ROOT_PATH}/ClinicalPositions" ]]; then
        rsync -carvl --progress "${ROOT_PATH}/ClinicalPositions" ${STORAGE_PATH}
    else
        ErrorLog "$SAMPLEID" "${ROOT_PATH}/ClinicalPositions does not exist!";
    fi

    # Copy AnnovarOutput
    if [[ -d "${ROOT_PATH}/AnnovarOutput" ]]; then
        rsync -carvl --progress "${ROOT_PATH}/AnnovarOutput" ${STORAGE_PATH}
    else
        ErrorLog "$SAMPLEID" "${ROOT_PATH}/AnnovarOutput does not exist!";
    fi

    # Copy PindelAnnovarOutput
    if [[ -d "${ROOT_PATH}/PindelAnnovarOutput" ]]; then
        rsync -carvl --progress "${ROOT_PATH}/PindelAnnovarOutput" ${STORAGE_PATH}
    else
        ErrorLog "$SAMPLEID" "${ROOT_PATH}/PindelAnnovarOutput does not exist!";
    fi

    # Copy vcfOutput
    if [[ -d "${ROOT_PATH}/vcfOutput" ]]; then
        rsync -carvl --progress "${ROOT_PATH}/vcfOutput" ${STORAGE_PATH}
    else
        ErrorLog "$SAMPLEID" "${ROOT_PATH}/vcfOutput does not exist!";
    fi

    # Copy MSIanalysis
    if [[ -d "${ROOT_PATH}/MSIanalysis" ]]; then
        rsync -carvl --progress "${ROOT_PATH}/MSIanalysis" ${STORAGE_PATH}
    else
        ErrorLog "$SAMPLEID" "${ROOT_PATH}/MSIanalysis does not exist!";
    fi

    # Copy ExtractedSNPs
    if [[ -d "${ROOT_PATH}/ExtractedSNPs" ]]; then
        rsync -carvl --progress "${ROOT_PATH}/ExtractedSNPs" ${STORAGE_PATH}
    else
        ErrorLog "$SAMPLEID" "${ROOT_PATH}/ExtractedSNPs does not exist!";
    fi

    # Copy ExtractedEGFR
    if [[ -d "${ROOT_PATH}/ExtractedEGFR" ]]; then
        rsync -carvl --progress "${ROOT_PATH}/ExtractedEGFR" ${STORAGE_PATH}
    else
        ErrorLog "$SAMPLEID" "${ROOT_PATH}/ExtractedEGFR does not exist!";
    fi

    # Copy Amplification
    if [[ -d "${ROOT_PATH}/Amplification" ]]; then
        rsync -carvl --progress "${ROOT_PATH}/Amplification" ${STORAGE_PATH}
    else
        ErrorLog "$SAMPLEID" "${ROOT_PATH}/Amplification does not exist!";
    fi

    # Copy FastQC
    if [[ -d "${ROOT_PATH}/FastQC" ]]; then
        rsync -carvl --progress "${ROOT_PATH}/FastQC" ${STORAGE_PATH}
    else
        ErrorLog "$SAMPLEID" "${ROOT_PATH}/FastQC does not exist!";
    fi

    # Copy SNPmania
    if [[ -d "${ROOT_PATH}/SNPmania" ]]; then
        rsync -carvl --progress "${ROOT_PATH}/SNPmania" ${STORAGE_PATH}
    else
        ErrorLog "$SAMPLEID" "${ROOT_PATH}/SNPmania does not exist!";
    fi

if [[ "$?" != "0" ]]; then
    ErrorLog "$SAMPLEID" "Failed in copying files to storage";
else
    SuccessLog "$SAMPLEID" "Passed copying files to storage";
fi
