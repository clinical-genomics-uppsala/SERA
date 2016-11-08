"""

Script for creating input file for sera pipeline

Elin Falk Sorqvist 20141114

"""

import argparse
import re
import os
import shutil
import sys
import module_locator


# Parse commandline
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
parser.add_argument('-pindelFlags', '--pindelFlags', help = 'Set parameters for pindel analysis. Default: -T 5 -x 2 -B 100', type = str, default = "-T 5 -x 2 -B 100")
parser.add_argument('-snpmaniaFlags', '--snpmaniaFlags', help = 'Set the SNPmania parameters. Default: -a 5 -q 20 -e 0 -am', type = str, default = "-a 5 -q 20 -e 0 -am")
parser.add_argument('-tumorNormalFlags', '--tumorNormalFlags', help = 'Set parameters for pre-annovar filtering when both tumor and normal are present. Default: -tminRD 20 -nminRD 20 -tminVarRatio 0.10 -nHomoRefRatio 0.95 -am 5', type = str, default = "-tminRD 20 -nminRD 20 -tminVarRatio 0.10 -nHomoRefRatio 0.95 -am 5")
parser.add_argument('-annovarFlags', '--annovarFlags', help = 'Set parameters for pre-annovar filtering when only tumor is present. Default: -minRD 20 -minVarRatio 0.01 -am 5', type = str, default = " -minRD 20 -minVarRatio 0.01 -am 5")
parser.add_argument('-annovarPlasmaFlags', '--annovarPlasmaFlags', help = 'Set parameters for pre-annovar filtering when only tumor is present. Default: -minRD 20 -minVarRatio 0.001', type = str, default = " -minRD 20 -minVarRatio 0.001")
parser.add_argument('-pindelAnnovarFlags', '--pindelAnnovarPlasmaFlags', help = 'Set parameters for pre-annovar filtering of pindel results. Default: -m 20 -v 0.01', type = str, default = "-m 20 -v 0.01")
parser.add_argument('-pindelAnnovarPlasmaFlags', '--pindelAnnovarFlags', help = 'Set parameters for pre-annovar filtering of pindel results. Default: -m 20 -v 0.001', type = str, default = "-m 20 -v 0.001")
parser.add_argument('-clinicalFlags', '--clinicalFlags', help = 'Set parameters for extracting clinical positions. Default: -minRD 30,300 -minVarRatio 0.01 -minAmpRD 5', type = str, default = "-minRD 30,300 -minVarRatio 0.01 -minAmpRD 5")
parser.add_argument('-clinicalPlasmaFlags', '--clinicalPlasmaFlags', help = 'Set parameters for extracting clinical positions. Default: -minRD 30,300 -minVarRatio 0.001', type = str, default = "-minRD 30,300 -minVarRatio 0.001")
parser.add_argument('-pindelClinicalFlags', '--pindelClinicalFlags', help = 'Set parameters for extracting indels from pindel in clinical genes. Default: -minRD 30,300 -minVarRatio 0.01', type = str, default = "-minRD 30,300 -minVarRatio 0.01")
parser.add_argument('-pindelClinicalPlasmaFlags', '--pindelClinicalPlasmaFlags', help = 'Set parameters for extracting indels from pindel in clinical genes. Default: -minRD 30,300 -minVarRatio 0.001', type = str, default = "-minRD 30,300 -minVarRatio 0.001")
parser.add_argument('-clinicalInfoFile', '--clinicalInfoFile', help = 'File with clinical hotspot and indel filenames per cancer type. If not wanted set to false! Default: clinicalCancerTypeFiles.txt', type = str, default = "clinicalCancerTypeFiles.txt")
parser.add_argument('-regionClinicalFlags', '--regionClinicalFlags', help = 'Set parameters for reporting all positions in whole regions. Default: -minRD 30,300', type = str, default = "-minRD 30,300")

args = parser.parse_args()
info = {}

clinicalInfo = {}
# Go through file with info about different clinical files
if not args.clinicalInfoFile.lower() == "false":
    # Using module module_locator to determine the path of the pythonscript directory
    my_path = module_locator.module_path()
    my_pathParts = my_path.split("/")  # Split path on /
    # Get the path of the SERA_script folder
    ciFilePath = ""
    for i in range(1, (len(my_pathParts) - 2)):
        ciFilePath += "/" + my_pathParts[i]
    # Add the path to the file with info about files corresponding to different cancer types
    ciFilePath += "/res/" + args.clinicalInfoFile

    # Go through the cancer type file and add info about the corresponding files to a directory
    with open(ciFilePath, 'r') as cifile:
        for line in cifile:
            if not re.match('$', line):  # If line is not empty start working
                line = line.strip()
                infoParts = line.split("=")
                cancer = infoParts[1].lower()
                clinicalInfo[cancer] = {}

                # Retrieve file info for each cancer type and add to dictionary
                for i in range (0, 6):
                    line = next(cifile)
                    line = line.strip()
                    infoParts = line.split("=")
                    # Check if the file given i - or false, then add false as filename
                    if re.match('-', infoParts[1]) or re.match('false', infoParts[1].lower()):
                        clinicalInfo[cancer][infoParts[0].lower()] = "false"
                    else:  # Add the filename
                        clinicalInfo[cancer][infoParts[0].lower()] = infoParts[1]

        if not cifile.closed:
            cifile.close()

with open(args.infile, 'r') as infile:
    count = 0
    # Go through the file line by line
    for line in infile:
        # Check so the line isn't empty or starts with ##

        if not re.match('^##', line) and not re.match('$', line) and not re.match('^Experiment', line):
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
            info[lineSplit[1]]['sNummer'] = "S" + str(count)
            info[lineSplit[1]]['design'] = lineSplit[2]
            info[lineSplit[1]]['refseq'] = lineSplit[3]
            info[lineSplit[1]]['cutadapt'] = lineSplit[4]
            info[lineSplit[1]]['type'] = lineSplit[5].lower().strip()
            info[lineSplit[1]]['tissue'] = lineSplit[6].lower().strip()
            info[lineSplit[1]]['barcodeI7'] = lineSplit[7]
            info[lineSplit[1]]['barcodeI5'] = lineSplit[8]

            if info[lineSplit[1]]['barcodeI7'] == "":
                info[lineSplit[1]]['barcodeI7'] = "false"
            if info[lineSplit[1]]['barcodeI5'] == "":
                info[lineSplit[1]]['barcodeI5'] = "false"


            if len(lineSplit) > 9 and args.normal:
                parser.print_usage()
                print("\nERROR: The normal has to be given either in the file or the commandline - NOT both!\n\n")
                sys.exit()
            else:
                if len(lineSplit) > 9:
                    info[lineSplit[1]]['normal'] = lineSplit[6]
                else:
                    if args.normal:
                        info[lineSplit[1]]['normal'] = args.normal
                    else:
                        parser.print_usage()
                        sys.exit("\nERROR: Either inputfile needs to have 7 columns with the normal given in the last column or the flag -normal has to be used!\n\n")
            # If clinical analysis is wanted add file names otherwise add false
            if not args.clinicalInfoFile.lower() == "false":
                # ## COLON
                if info[lineSplit[1]]['tissue'] == "colon" or info[lineSplit[1]]['tissue'] == "kolon" :
                    info[lineSplit[1]]['tissue'] = "colon"
                    # If the given hotspot file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo['colon']['hotspot']):
                        info[lineSplit[1]]['hotspot'] = "false"
                    else:
                        info[lineSplit[1]]['hotspot'] = "$FILE_PATH/refFiles/" + clinicalInfo["colon"]['hotspot']
                    # If the given indel file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["colon"]['indel']):
                        info[lineSplit[1]]['indel'] = "false"
                    else:
                        info[lineSplit[1]]['indel'] = "$FILE_PATH/refFiles/" + clinicalInfo["colon"]['indel']

                    # If the given region file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["colon"]['region']):
                        info[lineSplit[1]]['region'] = "false"
                    else:
                        info[lineSplit[1]]['region'] = "$FILE_PATH/refFiles/" + clinicalInfo["colon"]['region']

                    # If the given keep file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["colon"]['keep']):
                        info[lineSplit[1]]['keep'] = "false"
                    else:
                        info[lineSplit[1]]['keep'] = "$FILE_PATH/refFiles/" + clinicalInfo["colon"]['keep']

                    # If the given amplification file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["colon"]['amplification']):
                        info[lineSplit[1]]['amplification'] = "false"
                    else:
                        info[lineSplit[1]]['amplification'] = "$FILE_PATH/refFiles/" + clinicalInfo["colon"]['amplification']

                    # If the given background file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["colon"]['background']):
                        info[lineSplit[1]]['background'] = "false"
                    else:
                        info[lineSplit[1]]['background'] = "$FILE_PATH/refFiles/" + clinicalInfo["colon"]['background']

                # ## LUNG
                elif info[lineSplit[1]]['tissue'] == "lung" or info[lineSplit[1]]['tissue'] == "lunga":
                    info[lineSplit[1]]['tissue'] = "lung"
                    # If the given hotspot file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["lung"]['hotspot']):
                        info[lineSplit[1]]['hotspot'] = "false"
                    else:
                        info[lineSplit[1]]['hotspot'] = "$FILE_PATH/refFiles/" + clinicalInfo["lung"]['hotspot']
                    # If the given indel file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["lung"]['indel']):
                        info[lineSplit[1]]['indel'] = "false"
                    else:
                        info[lineSplit[1]]['indel'] = "$FILE_PATH/refFiles/" + clinicalInfo["lung"]['indel']

                    # If the given region file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["lung"]['region']):
                        info[lineSplit[1]]['region'] = "false"
                    else:
                        info[lineSplit[1]]['region'] = "$FILE_PATH/refFiles/" + clinicalInfo["lung"]['region']

                     # If the given keep file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["lung"]['keep']):
                        info[lineSplit[1]]['keep'] = "false"
                    else:
                        info[lineSplit[1]]['keep'] = "$FILE_PATH/refFiles/" + clinicalInfo["lung"]['keep']

                    # If the given amplification file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["lung"]['amplification']):
                        info[lineSplit[1]]['amplification'] = "false"
                    else:
                        info[lineSplit[1]]['amplification'] = "$FILE_PATH/refFiles/" + clinicalInfo["lung"]['amplification']

                    # If the given background file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["lung"]['background']):
                        info[lineSplit[1]]['background'] = "false"
                    else:
                        info[lineSplit[1]]['background'] = "$FILE_PATH/refFiles/" + clinicalInfo["lung"]['background']

                # ## GIST
                elif info[lineSplit[1]]['tissue'] == "gist":
                    info[lineSplit[1]]['tissue'] = "gist"
                    # If the given hotspot file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["gist"]['hotspot']):
                        info[lineSplit[1]]['hotspot'] = "false"
                    else:
                        info[lineSplit[1]]['hotspot'] = "$FILE_PATH/refFiles/" + clinicalInfo["gist"]['hotspot']
                    # If the given indel file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["gist"]['indel']):
                        info[lineSplit[1]]['indel'] = "false"
                    else:
                        info[lineSplit[1]]['indel'] = "$FILE_PATH/refFiles/" + clinicalInfo["gist"]['indel']

                    # If the given region file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["gist"]['region']):
                        info[lineSplit[1]]['region'] = "false"
                    else:
                        info[lineSplit[1]]['region'] = "$FILE_PATH/refFiles/" + clinicalInfo["gist"]['region']

                     # If the given keep file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["gist"]['keep']):
                        info[lineSplit[1]]['keep'] = "false"
                    else:
                        info[lineSplit[1]]['keep'] = "$FILE_PATH/refFiles/" + clinicalInfo["gist"]['keep']

                    # If the given amplification file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["gist"]['amplification']):
                        info[lineSplit[1]]['amplification'] = "false"
                    else:
                        info[lineSplit[1]]['amplification'] = "$FILE_PATH/refFiles/" + clinicalInfo["gist"]['amplification']

                    # If the given background file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["gist"]['background']):
                        info[lineSplit[1]]['background'] = "false"
                    else:
                        info[lineSplit[1]]['background'] = "$FILE_PATH/refFiles/" + clinicalInfo["gist"]['background']

                # ## OVARIAL
                elif info[lineSplit[1]]['tissue'] == "ovarial":
                    info[lineSplit[1]]['tissue'] = "ovarial"
                    # If the given hotspot file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["ovarial"]['hotspot']):
                        info[lineSplit[1]]['hotspot'] = "false"
                    else:
                        info[lineSplit[1]]['hotspot'] = "$FILE_PATH/refFiles/" + clinicalInfo["ovarial"]['hotspot']
                    # If the given indel file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["ovarial"]['indel']):
                        info[lineSplit[1]]['indel'] = "false"
                    else:
                        info[lineSplit[1]]['indel'] = "$FILE_PATH/refFiles/" + clinicalInfo["ovarial"]['indel']

                    # If the given region file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["ovarial"]['region']):
                        info[lineSplit[1]]['region'] = "false"
                    else:
                        info[lineSplit[1]]['region'] = "$FILE_PATH/refFiles/" + clinicalInfo["ovarial"]['region']

                     # If the given keep file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["ovarial"]['keep']):
                        info[lineSplit[1]]['keep'] = "false"
                    else:
                        info[lineSplit[1]]['keep'] = "$FILE_PATH/refFiles/" + clinicalInfo["ovarial"]['keep']

                    # If the given amplification file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["ovarial"]['amplification']):
                        info[lineSplit[1]]['amplification'] = "false"
                    else:
                        info[lineSplit[1]]['amplification'] = "$FILE_PATH/refFiles/" + clinicalInfo["ovarial"]['amplification']

                    # If the given background file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["ovarial"]['background']):
                        info[lineSplit[1]]['background'] = "false"
                    else:
                        info[lineSplit[1]]['background'] = "$FILE_PATH/refFiles/" + clinicalInfo["ovarial"]['background']

                # ## MELANOM
                elif info[lineSplit[1]]['tissue'] == "melanom" or info[lineSplit[1]]['tissue'] == "melanoma":
                    info[lineSplit[1]]['tissue'] = "melanom"
                    # If the given hotspot file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["melanom"]['hotspot']):
                        info[lineSplit[1]]['hotspot'] = "false"
                    else:
                        info[lineSplit[1]]['hotspot'] = "$FILE_PATH/refFiles/" + clinicalInfo["melanom"]['hotspot']
                    # If the given indel file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["melanom"]['indel']):
                        info[lineSplit[1]]['indel'] = "false"
                    else:
                        info[lineSplit[1]]['indel'] = "$FILE_PATH/refFiles/" + clinicalInfo["melanom"]['indel']

                    # If the given region file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["melanom"]['region']):
                        info[lineSplit[1]]['region'] = "false"
                    else:
                        info[lineSplit[1]]['region'] = "$FILE_PATH/refFiles/" + clinicalInfo["ovarial"]['region']

                     # If the given keep file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["melanom"]['keep']):
                        info[lineSplit[1]]['keep'] = "false"
                    else:
                        info[lineSplit[1]]['keep'] = "$FILE_PATH/refFiles/" + clinicalInfo["melanom"]['keep']

                    # If the given amplification file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["melanom"]['amplification']):
                        info[lineSplit[1]]['amplification'] = "false"
                    else:
                        info[lineSplit[1]]['amplification'] = "$FILE_PATH/refFiles/" + clinicalInfo["melanom"]['amplification']

                    # If the given background file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["melanom"]['background']):
                        info[lineSplit[1]]['background'] = "false"
                    else:
                        info[lineSplit[1]]['background'] = "$FILE_PATH/refFiles/" + clinicalInfo["melanom"]['background']

                else:
                    print ("\nERROR: Unknown cancer type given in input file" + info[lineSplit[1]]['tissue'] + ", so far only lung, colon, GIST, melanoma and ovarial are supported!\n\n")
                    sys.exit()
            else:
                info[lineSplit[1]]['hotspot'] = "false"
                info[lineSplit[1]]['indel'] = "false"
                info[lineSplit[1]]['region'] = "false"
                info[lineSplit[1]]['keep'] = "false"
                info[lineSplit[1]]['amplification'] = "false"
                info[lineSplit[1]]['background'] = "false"

    if not infile.closed:
        infile.close()

expFolder = ""
count = 1
infoSort = sorted(info)
for sample in infoSort:
    if count == 1:
        expFolder = info[sample]['exp'] + "_" + sample
        count += 1
    else:
        expFolder += "-" + sample

rawPath = "/proj/" + args.project + "/private/" + info[sample]['exp'] + "_rawdata"
filePath = "/proj/" + args.project + "/nobackup/private/" + info[sample]['exp']

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
    outfile.write("export ANNOVAR_PLASMA_FLAGS=\"" + args.annovarPlasmaFlags + "\";\n")
    outfile.write("export PINDEL_ANNOVAR_PLASMA_FLAGS=\"" + args.pindelAnnovarPlasmaFlags + "\";\n")

    outfile.write("\n## Clinical filtering\n")
    outfile.write("export CLINICAL_FLAGS=\"" + args.clinicalFlags + "\";\n")
    outfile.write("export PINDEL_CLINICAL_FLAGS=\"" + args.pindelClinicalFlags + "\";\n")
    outfile.write("export CLINICAL_PLASMA_FLAGS=\"" + args.clinicalPlasmaFlags + "\";\n")
    outfile.write("export PINDEL_CLINICAL__PLASMA_FLAGS=\"" + args.pindelClinicalPlasmaFlags + "\";\n")
    outfile.write("export REGION_CLINICAL_FLAGS=\"" + args.regionClinicalFlags + "\";\n")



    outfile.write('\n################################\n### SAMPLE SETTINGS\n################################\n\n')

    infoSort = sorted(info)
    for sample in infoSort:
        outfile.write("let COUNT=COUNT+1;\n")
        outfile.write ("SAMPLEID_ARR_[${COUNT}]=\"" + sample + "\";\n")
        outfile.write ("BARCODE_I7_ARR_[${COUNT}]=\"" + info[sample]['barcodeI7'] + "\";\n")
        outfile.write ("BARCODE_I5_ARR_[${COUNT}]=\"" + info[sample]['barcodeI5'] + "\";\n")
        outfile.write ("CUTADAPT_PREFIX_ARR_[${COUNT}]=\"" + info[sample]['cutadapt'] + "\";\n")
        outfile.write ("RAWDATA_PE1_ARR_[${COUNT}]=\"$RAW_PATH/" + sample + "_" + info[sample]['sNummer'] + "_L001_R1_001.fastq.gz" + "\";\n")
        outfile.write ("RAWDATA_PE2_ARR_[${COUNT}]=\"$RAW_PATH/" + sample + "_" + info[sample]['sNummer'] + "_L001_R2_001.fastq.gz" + "\";\n")
        outfile.write ("RAWDATA_INDEX_ARR_[${COUNT}]=\"false\";\n")
        outfile.write ("REFSEQ_ARR_[${COUNT}]=\"" + info[sample]['refseq'] + "\";\n")
        outfile.write ("ROIFILE_ARR_[${COUNT}]=\"$FILE_PATH/refFiles/" + info[sample]['design'] + "_Regions.bed" + "\";\n")
        outfile.write ("SELECTIONFILE_ARR_[${COUNT}]=\"$FILE_PATH/refFiles/" + info[sample]['design'] + "_Amplicons.bed" + "\";\n")
        outfile.write ("NORMAL_SAMPLEID_ARR_[${COUNT}]=\"" + info[sample]['normal'] + "\";\n")
        outfile.write ("HOTSPOTFILE_ARR_[${COUNT}]=\"" + info[sample]['hotspot'] + "\";\n")
        outfile.write ("INDELFILE_ARR_[${COUNT}]=\"" + info[sample]['indel'] + "\";\n")
        outfile.write ("REGIONFILE_ARR_[${COUNT}]=\"" + info[sample]['region'] + "\";\n")
        outfile.write ("KEEPFILE_ARR_[${COUNT}]=\"" + info[sample]['keep'] + "\";\n")
        outfile.write ("AMPLIFICATIONFILE_ARR_[${COUNT}]=\"" + info[sample]['amplification'] + "\";\n")
        outfile.write ("BACKGROUNDFILE_ARR_[${COUNT}]=\"" + info[sample]['background'] + "\";\n")
        outfile.write ("TYPE_ARR_[${COUNT}]=\"" + info[sample]['type'] + "\";\n")
        outfile.write ("TISSUE_ARR_[${COUNT}]=\"" + info[sample]['tissue'] + "\";\n")
        outfile.write ("\n")
    if not outfile.closed:
        outfile.close()
