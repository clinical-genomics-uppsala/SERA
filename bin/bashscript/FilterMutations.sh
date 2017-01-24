#!/bin/bash
#
# Script to run jSNPmania
##SBATCH --qos=short 
#SBATCH -p devcore  -n 1
#SBATCH -t 01:00:00


# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts Annovar ...";

# Check if the directory exists, if not create it
if [[ ! -d $ROOT_PATH/AnnovarOutput ]]; then 
    mkdir $ROOT_PATH/AnnovarOutput;
fi
if [[ ! -d $ROOT_PATH/FilteredMutations ]]; then 
    mkdir $ROOT_PATH/FilteredMutations;
fi 

ANNOVARFILE=""
PINDELANNOVARFILE=""
SNPMANIAFILE=""
ampMapped=""

if [[ ${NORMAL_SAMPLEID} != "false" ]]; then
	if [[ ${READS} == "true" ]]; then
        if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
            if [[ ${NORMAL_SAMPLEID} == "annovar" ]]; then
                if [[ -e $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput && $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput && $ROOT_PATH/SNPmania/${SAMPLEID}.variations]]; then
                    ANNOVARFILE="$ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput";
                    PINDELANNOVARFILE="$ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput";
                    SNPMANIAFILE="$ROOT_PATH/SNPmania/${SAMPLEID}.variations";
                
                elif [[ -e $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput && $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.singleSample.annovarOutput && $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations]]; then
                    ANNOVARFILE="$ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput";
                    PINDELANNOVARFILE="$ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput";
                    SNPMANIAFILE="$ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations";
                    ampMapped=" -a ";
                    
                else
                    ErrorLog "$SAMPLEID" "All inputfiles needed can not be found!";
                fi
            
                if [[ -e $ANNOVARFILE && $PINDELANNOVARFILE && $SNPMANIAFILE]]; then
                    if [[ ${TYPE} == "ovarial" ]]; then
                        cat $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/ExtractDepth_new.py -i ${HOTSPOTFILE} -o $ROOT_PATH/FilteredMutations/${SAMPLEID}.singleSample.filteredMutations -s ${SAMPLEID} -chr2nc $NC2chr -g 0.02 ${ampMapped} -t ${MAIN_TRANSCRIPT} -v ${SNPMANIAFILE} ${MUTATION_FLAGS}
                        echo "cat $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/ExtractDepth_new.py -i ${HOTSPOTFILE} -o $ROOT_PATH/FilteredMutations/${SAMPLEID}.singleSample.filteredMutations -s ${SAMPLEID} -chr2nc $NC2chr -g 0.02 ${ampMapped} -t ${MAIN_TRANSCRIPT} -v ${SNPMANIAFILE} ${MUTATION_FLAGS}"
                        
                    else
                        cat $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/ExtractDepth_new.py -i ${HOTSPOTFILE} -o $ROOT_PATH/FilteredMutations/${SAMPLEID}.singleSample.filteredMutations -s ${SAMPLEID} -chr2nc $NC2chr -g 0.01 ${ampMapped} -t ${MAIN_TRANSCRIPT} -v ${SNPMANIAFILE} ${MUTATION_FLAGS}
                        echo "cat $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/ExtractDepth_new.py -i ${HOTSPOTFILE} -o $ROOT_PATH/FilteredMutations/${SAMPLEID}.singleSample.filteredMutations -s ${SAMPLEID} -chr2nc $NC2chr -g 0.01 ${ampMapped} -t ${MAIN_TRANSCRIPT} -v ${SNPMANIAFILE} ${MUTATION_FLAGS}"
                        fi
                    fi
                else
                    ErrorLog "$SAMPLEID" "The annovar output file $ANNOVARFILE and/or pindel annovar output file $PINDELANNOVARFILE and/or SNPmania varition file $SNPMANIAFILE do NOT exist!!!";
                fi
            else
                    ErrorLog "$SAMPLEID" "Only implemented for NORMAL_SAMPLEID annovar so far!";
            fi
        else
            ErrorLog "$SAMPLEID" "Only supported for call_type h.sapiens so far!";
        fi
    else
        ErrorLog "$SAMPLEID" "READS has to be true to run the analysis!";
    fi
else
    ErrorLog "$SAMPLEID" "Normal_sampleid is false -> annovar is not run!";
fi

if [[ "$?" != "0" ]]; then
    ErrorLog "$SAMPLEID" "Failed in running Annovar";
else
    SuccessLog "$SAMPLEID" "Passed running Annovar";
fi