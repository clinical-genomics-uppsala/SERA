#!/bin/bash -l

#SBATCH -p core -n 1 
#SBATCH -t 2:00:00
##SBATCH --qos=short

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [ ! -d "${ROOT_PATH}/Molecules" ]; then
	mkdir ${ROOT_PATH}/Molecules;
fi

if [ ! -d "${ROOT_PATH}/SamBamFiles" ]; then
        mkdir ${ROOT_PATH}/SamBamFiles;
fi

SuccessLog "${SAMPLEID}" "Running BarD to add molecular tag to read";

# Check that molecules should be analysed
if [ ${MOLECULES} == "true" ]; then
	if [ ${PLATFORM} == "Illumina" ]; then
		if [[ "${MATE_PAIR}" == "true" && ${CALL_TYPE} == "h.sapiens" && ${DESIGN_TYPE} == "PCR" ]]; then
			# Change read name to fit in BarD
			if [[ ! -e  ${ROOT_PATH}/Molecules/${SAMPLEID}.index.fixed.fastq || ! -z $FORCE ]]; then 
				if [ -e ${ROOT_PATH}/filtereddata/${SAMPLEID}.index.fastq.gz ]; then
					zcat ${ROOT_PATH}/filtereddata/${SAMPLEID}.index.fastq.gz | awk -v name="${READ_ID_START}" 'BEGIN{FS=":"}{nm="^"name; if($1~nm){print $1"-"$3":"$4":"$5":"$6":"$7} else {print $0}}' /dev/stdin > ${ROOT_PATH}/Molecules/${SAMPLEID}.index.fixed.fastq;
				elif [ -e ${ROOT_PATH}/seqdata/${SAMPLEID}.index.fastq.gz ]; then
					zcat ${ROOT_PATH}/seqdata/${SAMPLEID}.index.fastq.gz | awk -v name="${READ_ID_START}" 'BEGIN{FS=":"}{nm="^"name; if($1~nm){print $1"AA114:"$2":"$5":"$6":"$7} else {print $0}}' /dev/stdin > ${ROOT_PATH}/Molecules/${SAMPLEID}.index.fixed.fastq;
				elif [ -e ${RAWDATA_INDEX} ]; then
					zcat ${RAWDATA_INDEX} | awk -v name="${READ_ID_START}" 'BEGIN{FS=":"}{nm="^"name; if($1~nm){print $1"A114:"$2":"$5":"$6":"$7} else {print $0}}' /dev/stdin > ${ROOT_PATH}/Molecules/${SAMPLEID}.index.fixed.fastq;
				else 
					ErrorLog "${SAMPLEID}" "Doesn't have an index read file!";
				fi
			else
				ErrorLog "${SAMPLEID}" "${ROOT_PATH}/Molecules/${SAMPLEID}.index.fixed.fastq already exists and force was NOT used!";
			fi
			
			# Check if the sorted file exists
			if [ -e ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.bam ]; then
				if [[ ! -e ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.fixed.bai || ! -z $FORCE ]]; then
					samtools view -h -F 4 ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.bam | awk -v name="${READ_ID_START}" 'BEGIN{OFS="\t"; sub(/@/,"",name); FS="\t"; }{nm="^"name; if($1~nm){split($1,a,":"); print a[1]"-"a[3]":"a[4]":"a[5]":"a[6]":"a[7]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12} else{print $0}}' /dev/stdin | samtools view -bS /dev/stdin > ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.fixed.bam;

#				samtools view -h -F 4 ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.bam | awk -v name="${READ_ID_START}" 'BEGIN{OFS="\t"; sub(/@/,"",name); FS="\t"; }{nm="^"name; if($1~nm){split($1,a,":"); print a[1]"AA114:"a[2]":"a[5]":"a[6]":"a[7]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12} else{print $0}}' /dev/stdin | samtools view -bS /dev/stdin > ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.fixed.bam;
				# Index bam-file
					
					if [ -e ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.fixed.bam ]; then 
						samtools index ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.fixed.bam ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.fixed.bai;
					else
						ErrorLog "${SAMPLEID}" "Couldn't generate the ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.fixed.bam file!";
					fi
				else
					ErrorLog "${SAMPLEID}" "The output file ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.fixed.bai already exist and force was NOT used!";
			else
				ErrorLog "${SAMPLEID}" "The sorted bam-file doesn't exists, run \"Create correct Alignment BAM-format\"...";
			fi
			
			# If all necessary files exist run AddBarcodeToBam - adding a field containing the molecular tag
			if [[ -e ${ROOT_PATH}/Molecules/${SAMPLEID}.index.fixed.fastq && -e ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.fixed.bam && -e ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.fixed.bai ]]; then
				java -classpath ${ROOT_PATH_BARD}/GenomeAnalysisTK.jar:${ROOT_PATH_BARD}/lib/Bard-20121008.jar org.broadinstitute.sting.gatk.CommandLineGATK -T AddBarcodeToBam -R ${GENOME_FASTA_REF} -I ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.fixed.bam -index ${ROOT_PATH}/Molecules/${SAMPLEID}.index.fixed.fastq -si 9 -ei 19 -o ${ROOT_PATH}/Molecules/${SAMPLEID}.h.sapiens.aligned.sorted.withIndex.bam;
				samtools index ${ROOT_PATH}/Molecules/${SAMPLEID}.h.sapiens.aligned.sorted.withIndex.bam ${ROOT_PATH}/Molecules/${SAMPLEID}.h.sapiens.aligned.sorted.withIndex.bai; 
				
				# Remove older files to save space
				#rm ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.fixed.bam;
				#rm ${ROOT_PATH}/SamBamFiles/${SAMPLEID}.h.sapiens.aligned.sorted.fixed.bai;
				#rm ${ROOT_PATH}/Molecules/${SAMPLEID}.index.fixed.fastq;
		
				# Change the bitwise flag to include 0x400 for reads with the same molecular tag
#				java -classpath ${ROOT_PATH_BARD}/GenomeAnalysisTK.jar:${ROOT_PATH_BARD}/lib/Bard-20121008.jar org.broadinstitute.sting.gatk.CommandLineGATK -T MarkDuplicatesByBarcode -R ${GENOME_FASTA_REF} -I ${ROOT_PATH}/Molecules/${SAMPLEID}.h.sapiens.aligned.sorted.withIndex.bam -o ${ROOT_PATH}/Molecules/${SAMPLEID}.h.sapiens.aligned.sorted.withIndex.duplicatesMarked.bam;
#				
#				# Remove older bam-files to save space
#				#rm ${ROOT_PATH}/Molecules/${SAMPLEID}.h.sapiens.aligned.sorted.withIndex.bam;
#				#rm ${ROOT_PATH}/Molecules/${SAMPLEID}.h.sapiens.aligned.sorted.withIndex.bai;
#	
#				# Remove duplicated reads
#				if [ -e ${ROOT_PATH}/Molecules/${SAMPLEID}.h.sapiens.aligned.sorted.withIndex.duplicatesMarked.bam ]; then
#					samtools view -h ${ROOT_PATH}/Molecules/${SAMPLEID}.h.sapiens.aligned.sorted.withIndex.duplicatesMarked.bam | awk '{if($1!~/^@/){if($2<1024){print $0}}else{print $0}}' | samtools view -bS /dev/stdin > ${ROOT_PATH}/Molecules/${SAMPLEID}.h.sapiens.aligned.sorted.molecules.bam
#					
#					# Remove older bam-files to save space
#					#rm ${ROOT_PATH}/Molecules/${SAMPLEID}.h.sapiens.aligned.sorted.withIndex.duplicatesMarked.bam;
#					#rm ${ROOT_PATH}/Molecules/${SAMPLEID}.h.sapiens.aligned.sorted.withIndex.duplicatesMarked.bai;
#	
#					if [ -e ${ROOT_PATH}/Molecules/${SAMPLEID}.h.sapiens.aligned.sorted.molecules.bam ]; then
#						 samtools index ${ROOT_PATH}/Molecules/${SAMPLEID}.h.sapiens.aligned.sorted.molecules.bam ${ROOT_PATH}/Molecules/${SAMPLEID}.h.sapiens.aligned.sorted.molecules.bai;
#					else
#						ErrorLog "${SAMPLEID}" "The bam-file with all duplicates removed could not be created...";
#					fi
#				else
#					ErrorLog "${SAMPLEID}" "The bam-file with duplicates marked could not be created...";
#				fi
#	
#
			else
				ErrorLog "${SAMPLEID}" "All files needed to run BarD was not present and/or could not be created...";
			fi
		fi
	fi
fi
if [ "$?" != "0" ]; then
	ErrorLog "${SAMPLEID}" "Failed in BarD...";
else
	SuccessLog "${SAMPLEID}" "Passed BarD";
fi		
