#!/bin/bash
#
# Script to export aligment archives as BAM-files.
#
#SBATCH -p core
#SBATCH -t 00:15:00
#SBATCH --qos=short

# Include functions
. $SERA_PATH/includes/logging.sh;

FILE="$ROOT_PATH/refFiles/plotFile.txt";

# if empty, create file
if [[ ! -e $FILE ]]; then
	echo -e "#Plotname type datafile1,title:datafile2,title\n#type=stenberg, stenberg or ericsson." > $FILE;
fi


# add reference and sample to file
function addReference {
	echo -e "${REFSEQ}.${1}.uniq.${2}\t" >> $FILE;
	# add if there are non-unique available
#	if [[ `ls -l $ROOT_PATH/plotValues/ | grep --exclude "uniq" | wc -l` == 0 ]]; then
#		echo -e "${REFSEQ}.${1}.${2}\t" >> $FILE;
#	fi;
}

# options to check for
types=(
"seqROI_seqRegion.uniq.nusbaum" "seqROI_seqRegion.nusbaum"
"seqROI_seqRegion.uniq.stenberg" "seqROI_seqRegion.stenberg"
"seqROI_seqROI.uniq.nusbaum" "seqROI_seqROI.nusbaum"
"seqROI_seqROI.uniq.stenberg" "seqROI_seqROI.stenberg"
);

# if there are no of this reference in file, append
if [[ `cat $FILE | grep $REFSEQ | wc -l` == 0 ]]; then

	for type in ${types[@]}; do

		PLOT=`echo $type | awk -F . '{print $NF}'`;

		# check if this type exist in folder
		if [[ `ls -l $ROOT_PATH/plotValues/ | grep ${type} | wc -l` != 0 ]]; then
			APPEND="$ROOT_PATH/plotValues/${SAMPLEID}.${type},${SAMPLEID}:";
			echo -e "${REFSEQ}.$type\t$PLOT\t$APPEND" >> $FILE;
		fi

	done;


# if there are reference already, just append sample
else

	# check so duplicates are not included
	if [[ `cat $FILE | grep $SAMPLEID | wc -l` == 0 ]]; then

		# check each line for each type
		for type in ${types[@]}; do

			CHECK="${REFSEQ}.${type}";
			PLOT=`echo $type | awk -F . '{print $NF}'`;
			APPEND="$ROOT_PATH/plotValues/${SAMPLEID}.${type},${SAMPLEID}:";

			# append with sed
			eval "sed -i 's|$REFSEQ.$type\t$PLOT\t|$REFSEQ.$type\t$PLOT\t$APPEND|' $FILE";

			# append line to correct lines by overwriting the old one
#			cat $FILE | awk '{if($1=="'$CHECK'"){print $0"'$APPEND'"}else{print $0}}' >> `echo ${FILE}`_cp;

			# overwrite
#			mv `echo $FILE`_cp $FILE;
#			rm `echo $FILE`_cp;

		done

	fi

fi

<<COMM

#look if the reference already exist in file
if [[ `cat $FILE | grep $REFSEQ | wc -l` == 0 ]]; then

	# append for sample
	APPEND="$ROOT_PATH/plotValues/${REFSEQ}.seqROI_seqRegion.uniq.nusbaum,${SAMPLEID}:";

	# add both nusbaum and stenberg for this reference
	echo -e "${REFSEQ}.seqROI_seqRegion.uniq.nusbaum\t" >> $FILE;
	echo -e "${REFSEQ}.seqROI_seqRegion.uniq.stenberg\t" >> $FILE;
	echo -e "${REFSEQ}.seqROI_seqROI.uniq.nusbaum\t" >> $FILE;
	echo -e "${REFSEQ}.seqROI_seqROI.uniq.stenberg\t" >> $FILE;

	# check if there are nonunique in path
	if [[ `ls -l $ROOT_PATH/plotValues/ | grep "uniq" | wc -l` != 0 ]]; then

		# add both nusbaum and stenberg for this reference
		echo -e "${REFSEQ}.seqROI_seqRegion.nusbaum\t" >> $FILE;
		echo -e "${REFSEQ}.seqROI_seqRegion.stenberg\t" >> $FILE;
		echo -e "${REFSEQ}.seqROI_seqROI.nusbaum\t" >> $FILE;
		echo -e "${REFSEQ}.seqROI_seqROI.stenberg\t" >> $FILE;

	fi

fi

# look if this sample already exists in plotfile-input
if [[ `cat $FILE | grep $SAMPLEID | wc -l` == 0 ]]; then

	# append line to correct lines by overwriting the old one
	APPEND="$ROOT_PATH/plotValues/${REFSEQ}.seqROI_seqRegion.uniq.nusbaum,${SAMPLEID}:";
	cat $FILE | awk '{split($0,a,"."); if(a[1]=="'$REFSEQ'"){print $0"'$APPEND'"}else{print $0}}' >> `echo ${FILE}`_cp;

	# overwrite
	rm $FILE;
	mv `echo $FILE`_cp $FILE;
	
fi

COMM


#check for control regions

#EGS38.ampregion.stenberg	stenberg	/data/seqdata/share/GATC_20110316/plotValues/EGS38_515_20110316.amproi_ampregion.stenberg,515,blue:/data/seqdata/share/GATC_20110316/plotValues/EGS38_516_20110316.amproi_ampregion.stenberg,516,black:/data/seqdata/share/GATC_20110316/plotValues/EGS38_517_20110316.amproi_ampregion.stenberg,517,red
#EGS38.ampregion.allSamples.stenberg	stenberg	/data/seqdata/share/GATC_20110316/plotValues/EGS38_515_20110316.amproi_ampregion.stenberg,515,blue:/data/seqdata/share/GATC_20110316/plotValues/EGS38_516_20110316.amproi_ampregion.stenberg,516,black:/data/seqdata/share/GATC_20110316/plotValues/EGS38_517_20110316.amproi_ampregion.stenberg,517,red:/data/seqdata/share/GATC_20110316/plotValues/EGS38_NA18508_1_20110316.amproi_ampregion.stenberg,NA18508 R1,blue:/data/seqdata/share/GATC_20110316/plotValues/EGS38_NA18508_2_20110316.amproi_ampregion.stenberg,NA18508 R2,black
#EGS38.amproi.ampregion.nusbaum	nusbaum	/data/seqdata/share/GATC_20110316/plotValues/EGS38_515_20110316.amproi_ampregion.nusbaum,515,blue:/data/seqdata/share/GATC_20110316/plotValues/EGS38_516_20110316.amproi_ampregion.nusbaum,516,black:/data/seqdata/share/GATC_20110316/plotValues/EGS38_517_20110316.amproi_ampregion.nusbaum,517,red
#EGS38.amproi.ampregion.allSamples.nusbaum	nusbaum	/data/seqdata/share/GATC_20110316/plotValues/EGS38_515_20110316.amproi_ampregion.nusbaum,515,blue:/data/seqdata/share/GATC_20110316/plotValues/EGS38_516_20110316.amproi_ampregion.nusbaum,516,black:/data/seqdata/share/GATC_20110316/plotValues/EGS38_517_20110316.amproi_ampregion.nusbaum,517,red:/data/seqdata/share/GATC_20110316/plotValues/EGS38_NA18508_1_20110316.amproi_ampregion.nusbaum,NA18508 R1,blue:/data/seqdata/share/GATC_20110316/plotValues/EGS38_NA18508_2_20110316.amproi_ampregion.nusbaum,NA18508 R2,black
#EGS38.amproi.amproi.nusbaum	nusbaum	/data/seqdata/share/GATC_20110316/plotValues/EGS38_515_20110316.amproi_amproi.nusbaum,515,blue:/data/seqdata/share/GATC_20110316/plotValues/EGS38_516_20110316.amproi_amproi.nusbaum,516,black:/data/seqdata/share/GATC_20110316/plotValues/EGS38_517_20110316.amproi_amproi.nusbaum,517,red
#EGS38.amproi.amproi.allSamples.nusbaum	nusbaum	/data/seqdata/share/GATC_20110316/plotValues/EGS38_515_20110316.amproi_amproi.nusbaum,515,blue:/data/seqdata/share/GATC_20110316/plotValues/EGS38_516_20110316.amproi_amproi.nusbaum,516,black:/data/seqdata/share/GATC_20110316/plotValues/EGS38_517_20110316.amproi_amproi.nusbaum,517,red:/data/seqdata/share/GATC_20110316/plotValues/EGS38_NA18508_1_20110316.amproi_amproi.nusbaum,NA18508 R1,blue:/data/seqdata/share/GATC_20110316/plotValues/EGS38_NA18508_2_20110316.amproi_amproi.nusbaum,NA18508 R2,black
