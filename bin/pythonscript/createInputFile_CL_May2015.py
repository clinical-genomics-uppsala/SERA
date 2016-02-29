"""

Script for creating input file for sera pipeline

Elin Falk Sorqvist 20141114

Edited by Claes Ladenvall 20150508
To handle files deposited by the SNP & SEQ platform for a Haloplex run

"""

import argparse
import re
import os
import shutil
import sys
import module_locator


# Parse commandline
#Eaxmple commandline
#/proj/a2013225/private/software/SERA_scripts_20150422/bin/pythonscript/createInputFile_CL.py
#createInputFile_CL.py -i INPUTFILE -p a2013223 -g PROJ -refDir /proj/a2013223/private/project/Manifest_etc/ -pindelFlags "-T 5 -x 2 -B 60 -j /proj/a2013223/private/project/Manifest_etc/Regions_for_Pindel_Haloplex_June2013.bed" -clinicalInfoFile false -n annovar
parser = argparse.ArgumentParser(description = "This script creates a SERA input file from a file with sample name and barcode")
parser.add_argument('-i', '--infile', help = 'Input file name', type = str, required = True)
parser.add_argument('-p', '--project', help = 'Name of the project', type = str, required = True)
parser.add_argument('-g', '--globals', help = 'Should the HOME or PROJ option be used', choices = ['HOME', 'PROJ'], type = str, required = True)
parser.add_argument('-refDir', '--refDir', help = 'If a directory for reference files is given, they will be copied to the analysis folder', type = str)
parser.add_argument('-n', '--normal', help = 'Name of normal sample to compare with.', type = str)
parser.add_argument('-s', '--software', help = 'Queing system for', type = str, default = "SLURM")
parser.add_argument('-f', '--fileFormat', help = 'File format of reference files. Default: bed', type = str, default = "bed")
parser.add_argument('-platform', '--platform', help = 'Platform for sequencing. Default: Illumina', type = str, default = "Illumina")
parser.add_argument('-readLength', '--readLength', help = 'Read length in sequencing. Default: 150', type = str, default = "150")
parser.add_argument('-matePair', '--matePair', help = 'If sequencing is run with mate pair set true otherwise false. Default: true', type = str, default = "true")
parser.add_argument('-designType', '--designType', help = 'Set the type of design, PCR or MDA. Default: PCR', type = str, default = "PCR")
parser.add_argument('-alignerFlags', '--alignerFlags', help = 'Set the alignment parameters. Default: -hs 15 -act 35 -mhp 100 -ms 7 -p 8 -minp 0.95 -mmp 0.05', type = str, default = "-hs 15 -act 35 -mhp 100 -ms 7 -p 8 -minp 0.95 -mmp 0.05")
parser.add_argument('-pindelFlags', '--pindelFlags', help = 'Set parameters for pindel analysis. Default: -T 5 -x 2 -B 25', type = str, default = "-T 5 -x 2 -B 25")
parser.add_argument('-snpmaniaFlags', '--snpmaniaFlags', help = 'Set the SNPmania parameters. Default: -a 5 -q 20 -e 0 -am', type = str, default = "-a 5 -q 20 -e 0 -am")
parser.add_argument('-tumorNormalFlags', '--tumorNormalFlags', help = 'Set parameters for pre-annovar filtering when both tumor and normal are present. Default: -tminRD 20 -nminRD 20 -tminVarRatio 0.10 -nHomoRefRatio 0.95 -am 5', type = str, default = "-tminRD 20 -nminRD 20 -tminVarRatio 0.10 -nHomoRefRatio 0.95 -am 5")
parser.add_argument('-annovarFlags', '--annovarFlags', help = 'Set parameters for pre-annovar filtering when only tumor is present. Default: -minRD 20 -minVarRatio 0.01 -am 5', type = str, default = " -minRD 20 -minVarRatio 0.01 -am 5")
parser.add_argument('-clinicalFlags', '--clinicalFlags', help = 'Set parameters for extracting clinical positions. Default: -minRD 30,100,300,500 -minVarRatio 0.01 -minAmpRD 5', type = str, default = "-minRD 30,100,300,500 -minVarRatio 0.01 -minAmpRD 5")
parser.add_argument('-pindelClinicalFlags', '--pindelClinicalFlags', help = 'Set parameters for extracting indels from pindel in clinical genes. Default: -minRD 30,100,300,500 -minVarRatio 0.01', type = str, default = "-minRD 30,100,300,500 -minVarRatio 0.01")
parser.add_argument('-pindelAnnovarFlags', '--pindelAnnovarFlags', help = 'Set parameters for pre-annovar filtering of pindel results. Default: -m 20 -v 0.01', type = str, default = "-m 20 -v 0.01")
parser.add_argument('-clinicalInfoFile', '--clinicalInfoFile', help = 'File with clinical hotspot and indel filenames per cancer type. If not wanted set to false! Default: clinicalCancerTypeFiles.txt', type = str, default = "clinicalCancerTypeFiles.txt")

args = parser.parse_args()
info = {}

clinicalInfo={}
# Go through file with info about different clinical files
if not args.clinicalInfoFile.lower() == "false":
    # Using module module_locator to determine the path of the pythonscript directory
    my_path = module_locator.module_path()
    my_pathParts = my_path.split("/") # Split path on /
    # Get the path of the SERA_script folder
    ciFilePath = ""
    for i in range(1,(len(my_pathParts)-2)):
        ciFilePath += "/"+my_pathParts[i]
    # Add the path to the file with info about files corresponding to different cancer types
    ciFilePath += "/res/"+args.clinicalInfoFile
    
    # Go through the cancer type file and add info about the corresponding files to a directory
    with open(ciFilePath, 'r') as cifile:
        for line in cifile:
            if not re.match('$', line): # If line is not empty start working
                line = line.strip()
                infoParts=line.split("=")
                cancer=infoParts[1].lower()
                clinicalInfo[cancer] = {}
               
                # Retrieve file info for each cancer type and add to dictionary
                for i in range (0,2):
                    line =next(cifile)
                    line = line.strip()
                    infoParts=line.split("=")
                    #Check if the file given i - or false, then add false as filename
                    if re.match('-', infoParts[1]) or re.match('false', infoParts[1].lower()): 
                        clinicalInfo[cancer][infoParts[0].lower()] = "false"
                    else: # Add the filename
                        clinicalInfo[cancer][infoParts[0].lower()] = infoParts[1]
            
        if not cifile.closed:
            cifile.close()
    
with open(args.infile, 'r') as infile:
    count = 0
    # Go through the file line by line
    for line in infile:
        # Check so the line isn't empty or starts with ##

        if not re.match('^##', line) and not re.match('$', line) and not re.match('^experiment', line.lower()):
            count += 1
            line = line.rstrip('\r\n')  # Remove new line character
            splitPattern = '\t'
            if re.search(";", line):
                splitPattern = ';'
            elif re.search(",", line):
                splitPattern = ','
            lineSplit = line.split(splitPattern)  # Split line by tab
            info[lineSplit[1]] = {}
            info[lineSplit[1]]['exp'] = lineSplit[0]
            info[lineSplit[1]]['index'] = lineSplit[2]
            info[lineSplit[1]]['barcode'] = lineSplit[3]
            info[lineSplit[1]]['sNummer'] = "S" + str(count)
            info[lineSplit[1]]['design'] = lineSplit[4]
            info[lineSplit[1]]['refseq'] = lineSplit[5]
            info[lineSplit[1]]['type'] = lineSplit[6].lower().strip()
            info[lineSplit[1]]['filename'] = lineSplit[7]
            
            if len(lineSplit) > 8 and args.normal:
                parser.print_usage()
                print("\nERROR: The normal has to be given either in the file or the commandline - NOT both!\n\n")
                sys.exit()
            else:
                if len(lineSplit) > 8:
                    info[lineSplit[1]]['normal'] = lineSplit[8]
                else:
                    if args.normal:
                        info[lineSplit[1]]['normal'] = args.normal
                    else:
                        parser.print_usage()
                        sys.exit("\nERROR: Either inputfile needs to have 9 columns with the normal given in the last column or the flag -normal has to be used!\n\n")
            # If clinical analysis is wanted add file names otherwise add false          
            if not args.clinicalInfoFile.lower() == "false":
                if info[lineSplit[1]]['type'] == "colon" or info[lineSplit[1]]['type'] == "kolon" :
                    # If the given hotspot file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo['colon']['hotspot']):
                        info[lineSplit[1]]['hotspot'] = "false"
                    else:
                        info[lineSplit[1]]['hotspot'] = "$FILE_PATH/refFiles/"+clinicalInfo["colon"]['hotspot']
                    # If the given indel file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["colon"]['indel']):
                        info[lineSplit[1]]['indel'] = "false"
                    else:
                        info[lineSplit[1]]['indel'] = "$FILE_PATH/refFiles/"+clinicalInfo["colon"]['indel']
                    
                elif info[lineSplit[1]]['type'] == "lung" or info[lineSplit[1]]['type'] == "lunga":
                    # If the given hotspot file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["lung"]['hotspot']):
                        info[lineSplit[1]]['hotspot'] = "false"
                    else:
                        info[lineSplit[1]]['hotspot'] = "$FILE_PATH/refFiles/"+clinicalInfo["lung"]['hotspot']
                    # If the given indel file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["lung"]['indel']):
                        info[lineSplit[1]]['indel'] = "false"
                    else:
                        info[lineSplit[1]]['indel'] = "$FILE_PATH/refFiles/"+clinicalInfo["lung"]['indel']
                    
                else:
                    print ("\nERROR: Unknown cancer type given in input file"+info[lineSplit[1]]['type']+", so far only Lung and Colon are supported!\n\n")
                    sys.exit()
            else:
                info[lineSplit[1]]['hotspot'] = "false"
                info[lineSplit[1]]['indel'] = "false"
                
    if not infile.closed:
        infile.close()

expFolder = ""
count = 1
infoSort = sorted(info)
folder = ""
for sample in infoSort:
    if count == 1:
        expFolder = info[sample]['exp'] + "_" + sample
        folder = info[sample]['exp']
        count += 1
    else:
        expFolder += "-" + sample
        folder = info[sample]['exp']

#Set raw data folder
rawPath = "/proj/" + args.project + "/private/project/2015/" + folder + "_rawdata"
filePath = "/proj/" + args.project + "/nobackup/private/project/2015/" + folder

if not os.path.exists(rawPath):
    os.makedirs(rawPath)
if not os.path.exists(filePath):
    os.makedirs(filePath)

# Check if a reference file dir is given, if so copy the files to refFiles in the analysis folder
if args.refDir:
    # List files in the reference file directory
    refDir_files = os.listdir(args.refDir)
    refFilePath = filePath + "/refFiles"  # Setting output folder path
    # Check if the output folder exists, if not create
    if not os.path.exists(refFilePath):
        os.makedirs(refFilePath)
    # Go through all files in the reference file dir
    for file_name in refDir_files:
        full_file_name = os.path.join(args.refDir, file_name)
        # Check that it is a file
        if (os.path.isfile(full_file_name)):
             destFile = refFilePath + "/" + file_name  # Set the destination file name
             shutil.copy(full_file_name, destFile)  # Copy



output = filePath + "/inputFile"
with (open(output, mode = 'w'))as outfile:
    outfile.write('################################\n### SEQUENCING-SPECIFIC SETTINGS\n################################\n\n# These settings are shared for all samples included in this analysis\n\n')
    outfile.write("SOFTWARE=\"" + args.software + "\";\n")
    outfile.write("UPPNEX_PROJECT_ID=\"" + args.project + "\";\n")
    outfile.write("GLOBALS=\"" + args.globals + "\";\n")
    outfile.write("RAW_PATH=\"" + rawPath + "\";\n")
    outfile.write("FILE_PATH=\"" + filePath + "\";\n")
    outfile.write("\nexport FILE_FORMAT=\"" + args.fileFormat + "\";\n")
    outfile.write("export READS=\"true\";\n")
    outfile.write("export MOLECULES=\"false\";\n")
    outfile.write("\n## Alignment\n")
    outfile.write("export PLATFORM=\"" + args.platform + "\" # SOLiD, Illumina or IonTorrent\n")
    outfile.write("export READ_LENGTH=\"" + args.readLength + "\";\n")
    outfile.write("export MATE_PAIR=\"" + str(args.matePair) + "\"; # true or false\n")
    outfile.write("export DESIGN_TYPE=\"" + args.designType + "\"; # PCR or MDA\n")
    outfile.write("export ALIGNERFLAGS_GENOME=\"" + args.alignerFlags + "\";\n")
    outfile.write("export PINDEL_FLAGS=\"" + args.pindelFlags + "\";\n")
    outfile.write("\n## Filtering\n")
    outfile.write("export CALL_TYPE=\"h.sapiens\";\n")
    outfile.write("export SNPMANIAFLAGS=\"" + args.snpmaniaFlags + "\";\n")
    outfile.write("export TUMOR_NORMAL_FLAGS=\"" + args.tumorNormalFlags + "\";\n")
    outfile.write("export ANNOVAR_FLAGS=\"" + args.annovarFlags + "\";\n")
    outfile.write("export PINDEL_ANNOVAR_FLAGS=\"" + args.pindelAnnovarFlags + "\";\n")
    
    outfile.write("\n## Clinical filtering\n")
    outfile.write("export CLINICAL_FLAGS=\"" + args.clinicalFlags + "\";\n")
    outfile.write("export PINDEL_CLINICAL_FLAGS=\"" + args.pindelClinicalFlags + "\";\n")
   
    outfile.write('\n################################\n### SAMPLE SETTINGS\n################################\n\n')

    infoSort = sorted(info)
    for sample in infoSort:
        outfile.write("let COUNT=COUNT+1;\n")
        outfile.write ("SAMPLEID_ARR_[${COUNT}]=\"" + sample + "\";\n")
        outfile.write ("SEQUENCING_TAG_ARR_[${COUNT}]=\"" + info[sample]['barcode'] + "\";\n")
        #Make sure to get the correct lane information!
        outfile.write ("RAWDATA_PE1_ARR_[${COUNT}]=\"$RAW_PATH/" + info[sample]['index'] + "_" + info[sample]['barcode'] + "_L007_R1_001.fastq.gz" + "\";\n")
        outfile.write ("RAWDATA_PE2_ARR_[${COUNT}]=\"$RAW_PATH/" + info[sample]['index'] + "_" + info[sample]['barcode'] + "_L007_R2_001.fastq.gz" + "\";\n")
        outfile.write ("RAWDATA_INDEX_ARR_[${COUNT}]=\"false\";\n")
        outfile.write ("REFSEQ_ARR_[${COUNT}]=\"" + info[sample]['refseq'] + "\";\n")
        outfile.write ("ROIFILE_ARR_[${COUNT}]=\"$FILE_PATH/refFiles/" + info[sample]['design'] + "_Regions.bed" + "\";\n")
        outfile.write ("SELECTIONFILE_ARR_[${COUNT}]=\"$FILE_PATH/refFiles/" + info[sample]['design'] + "_Amplicons.bed" + "\";\n")
        outfile.write ("NORMAL_SAMPLEID_ARR_[${COUNT}]=\"" + info[sample]['normal'] + "\";\n")
        outfile.write("HOTSPOTFILE_ARR_[${COUNT}]=\"" +  info[sample]['hotspot'] + "\";\n")
        outfile.write("INDELFILE_ARR_[${COUNT}]=\"" +  info[sample]['indel'] + "\";\n")
        outfile.write("\n")
    if not outfile.closed:
        outfile.close()
