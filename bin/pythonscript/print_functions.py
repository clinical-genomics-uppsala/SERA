from variant_functions import *
import re

# Print variants, both bwa and pindel on hotspot positions. If now variant found print base and it's total read depth
def printHotspots(hotspots, minRDs, sample, transcripts):
    for chrom in hotspots['hotspot']:
            for start in hotspots['hotspot'][chrom]:
                for end in hotspots['hotspot'][chrom][start]:

                    # If there are mutations from bwa print them
                    if not re.match("no", str(hotspots['hotspot'][chrom][start][end]['bwa'])):
                        pin = "no"
                        # check if there are also pindel variants present at this position
                        if not re.match("no", str(hotspots['hotspot'][chrom][start][end]['pindel'])):
                            pin = "yes"
                        # Go through all bwa variants
                        for bwaVar in hotspots['hotspot'][chrom][start][end]['bwa']:

                            if re.match("yes", pin):  # If there are pindel variants available extract them
                                for pindelVar in hotspots['hotspot'][chrom][start][end]['pindel']:
                                    # Check that there is no pindel variant which is identical to the bwa variant
                                    if not (bwaVar[2] == pindelVar[2] and bwaVar[3] == pindelVar[3] and bwaVar[4] == pindelVar[4] and bwaVar[5] == pindelVar[5]):
                                        aa, cds, accNum, exonicType = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information
                                        found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                                        if re.match("-", found):
                                            found = "yes"
                                        print (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + hotspots['hotspot'][chrom][start][end]['comment'] + "\thotspot\t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(bwaVar[28]) + "\t" + str(bwaVar[29]) + "\t" + str(bwaVar[26]) + "\t" + str(bwaVar[27]) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(bwaVar[31]) + "\t" + str(bwaVar[30]) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]))
                            else:  # If no pindel variants are available print the variant
                                aa, cds, accNum, exonicType = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information
                                found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                                if re.match("-", found):
                                    found = "yes"

                                print (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + hotspots['hotspot'][chrom][start][end]['comment'] + "\thotspot\t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(bwaVar[28]) + "\t" + str(bwaVar[29]) + "\t" + str(bwaVar[26]) + "\t" + str(bwaVar[27]) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(bwaVar[31]) + "\t" + str(bwaVar[30]) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]))


                    # Check if there are pindel variants available
                    if not re.match("no", str(hotspots['hotspot'][chrom][start][end]['pindel'])):

                        for pindelVar in hotspots['hotspot'][chrom][start][end]['pindel']:

                            aa, cds, accNum, exonicType = getTranscriptInfo(pindelVar, transcripts)

                            # Check the level of reads. Below the lowest => -, between first and second => low, above highest => ok
                            found, readLevel = getReadLevel(minRDs, hotspots['hotspot'][chrom][start][end]['rd'][0])
                            if re.match("-", found):
                                found = "yes"

                            print (str(pindelVar[0]) + "\t" + str(pindelVar[6]) + "\t" + exonicType + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + hotspots['hotspot'][chrom][start][end]['comment'] + "\thotspot\t" + found + "\t" + readLevel + "\t" + str(pindelVar[12]) + "\t" + str(pindelVar[10]) + "\t" + str(pindelVar[11]) + "\t" + str(pindelVar[9]) + "\t" + str(pindelVar[14]) + "\t" + str(pindelVar[13]) + "\t" + str(pindelVar[16]) + "\t" + str(pindelVar[15]) + "\t" + str(pindelVar[17]) + "\t" + str(pindelVar[18]) + "\t" + str(pindelVar[19]) + "\t" + str(pindelVar[28]) + "\t" + str(pindelVar[29]) + "\t" + str(pindelVar[26]) + "\t" + str(pindelVar[27]) + "\t" + str(pindelVar[20]) + "\t" + str(pindelVar[21]) + "\t" + str(pindelVar[22]) + "\t" + str(pindelVar[23]) + "\t" + str(pindelVar[24]) + "\t" + str(pindelVar[25]) + "\t" + str(pindelVar[31]) + "\t" + str(pindelVar[30]) + "\t" + str(pindelVar[1]) + "\t" + str(pindelVar[2]) + "\t" + str(pindelVar[3]) + "\t" + str(pindelVar[4]) + "\t" + str(pindelVar[5]) + "\t" + str(pindelVar[32]))


                    if re.match("no", str(hotspots['hotspot'][chrom][start][end]['bwa'])) and re.match("no", str(hotspots['hotspot'][chrom][start][end]['pindel'])):  # If no variant is added print position
                        for c in range(0, int(end) - int(start) + 1):
                            found = readLevel = rd = "NA"
                            if isinstance(hotspots['hotspot'][chrom][start][end]['rd'], list):  # check that the read depths are within an array, otherwise it's not in design

                                if re.match("-", str(hotspots['hotspot'][chrom][start][end]['rd'][c])):  # if read depth is - the position is not within design
                                    found = "not in design"
                                    readLevel = "-"
                                    rd = "-"
                                else:  # read depth exists
                                    found, readLevel = getReadLevel (minRDs, int(hotspots['hotspot'][chrom][start][end]['rd'][c]))
                                    if re.match("-", found):
                                        found = "no"
                                    rd = getReadLevel (minRDs, int(hotspots['hotspot'][chrom][start][end]['rd'][c]))
                            else:  # position is not within design
                                found = "not in design"
                                readLevel = "-"
                                rd = "-"

                            print (sample + "\t" + str(hotspots['hotspot'][chrom][start][end]['gene']) + "\t-\t" + str(hotspots['hotspot'][chrom][start][end]['aa']) + "\t" + str(hotspots['hotspot'][chrom][start][end]['cds']) + "\t" + str(hotspots['hotspot'][chrom][start][end]['accNum']) + "\t" + str(hotspots['hotspot'][chrom][start][end]['comment']) + "\thotspot\t" + found + "\t" + readLevel + "\t" + str(hotspots['hotspot'][chrom][start][end]['rd'][c]) + "\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t" + str(chrom) + "\t" + str(start) + "\t" + str(end) + "\t-\t-\t-")

# print variants within region, both bwa pindel and the total read depth for positions without variant.
def printRegionAll(hotspots, minRDs, sample, transcripts):
    for chrom in hotspots['region_all']:
        for start in hotspots['region_all'][chrom]:
            for end in hotspots['region_all'][chrom][start]:
                # If there are mutations from bwa print them
                if not re.match("no", str(hotspots['region_all'][chrom][start][end]['bwa'])):
                    pin = "no"
                    # check if there are also pindel variants present at this position
                    if not re.match("no", str(hotspots['region_all'][chrom][start][end]['pindel'])):
                        pin = "yes"
                    # Go through all bwa variants
                    for bwaVar in hotspots['region_all'][chrom][start][end]['bwa']:

                        if re.match("yes", pin):  # If there are pindel variants available extract them
                            for pindelVar in hotspots['region_all'][chrom][start][end]['pindel']:
                                # Check that there is no pindel variant which is identical to the bwa variant
                                if not (bwaVar[2] == pindelVar[2] and bwaVar[3] == pindelVar[3] and bwaVar[4] == pindelVar[4] and bwaVar[5] == pindelVar[5]):
                                    aa, cds, accNum, exonicType = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information
                                    found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                                    if re.match("-", found):
                                        found = "yes"
                                    print (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + hotspots['region_all'][chrom][start][end]['comment'] + "\tregion_all\t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(bwaVar[28]) + "\t" + str(bwaVar[29]) + "\t" + str(bwaVar[26]) + "\t" + str(bwaVar[27]) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(bwaVar[31]) + "\t" + str(bwaVar[30]) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]))
                        else:  # If no pindel variants are available print the variant
                            aa, cds, accNum, exonicType = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information
                            found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                            if re.match("-", found):
                                found = "yes"

                            print (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + hotspots['region_all'][chrom][start][end]['comment'] + "\tregion_all\t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(bwaVar[28]) + "\t" + str(bwaVar[29]) + "\t" + str(bwaVar[26]) + "\t" + str(bwaVar[27]) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(bwaVar[31]) + "\t" + str(bwaVar[30]) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]))


                # Check if there are pindel variants available
                if not re.match("no", str(hotspots['region_all'][chrom][start][end]['pindel'])):

                    for pindelVar in hotspots['region_all'][chrom][start][end]['pindel']:

                        aa, cds, accNum, exonicType = getTranscriptInfo(pindelVar, transcripts)

                        # Check the level of reads. Below the lowest => -, between first and second => low, above highest => ok
                        found, readLevel = getReadLevel(minRDs, hotspots['region_all'][chrom][start][end]['rd'][0])
                        if re.match("-", found):
                            found = "yes"

                        print (str(pindelVar[0]) + "\t" + str(pindelVar[6]) + "\t" + exonicType + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + hotspots['region_all'][chrom][start][end]['comment'] + "\tregion_allt\t" + found + "\t" + readLevel + "\t" + str(pindelVar[12]) + "\t" + str(pindelVar[10]) + "\t" + str(pindelVar[11]) + "\t" + str(pindelVar[9]) + "\t" + str(pindelVar[14]) + "\t" + str(pindelVar[13]) + "\t" + str(pindelVar[16]) + "\t" + str(pindelVar[15]) + "\t" + str(pindelVar[17]) + "\t" + str(pindelVar[18]) + "\t" + str(pindelVar[19]) + "\t" + str(pindelVar[28]) + "\t" + str(pindelVar[29]) + "\t" + str(pindelVar[26]) + "\t" + str(pindelVar[27]) + "\t" + str(pindelVar[20]) + "\t" + str(pindelVar[21]) + "\t" + str(pindelVar[22]) + "\t" + str(pindelVar[23]) + "\t" + str(pindelVar[24]) + "\t" + str(pindelVar[25]) + "\t" + str(pindelVar[31]) + "\t" + str(pindelVar[30]) + "\t" + str(pindelVar[1]) + "\t" + str(pindelVar[2]) + "\t" + str(pindelVar[3]) + "\t" + str(pindelVar[4]) + "\t" + str(pindelVar[5]) + "\t" + str(pindelVar[32]))

                # Check if neither bwa nor pindel variant exists
                if re.match("no", str(hotspots['region_all'][chrom][start][end]['bwa'])) and re.match("no", str(hotspots['region_all'][chrom][start][end]['pindel'])):  # If no variant is added print position
                    test = int(end) - int(start) + 1
                    for c in range(0, int(end) - int(start) + 1):
                        found = readLevel = rd = "NA"
                        if isinstance(hotspots['region_all'][chrom][start][end]['rd'], list):  # check that the read depths are within an array, otherwise it's not in design

                            if re.match("-", str(hotspots['region_all'][chrom][start][end]['rd'][c])):  # if read depth is - the position is not within design
                                found = "not in design"
                                readLevel = "-"
                                rd = "-"
                            else:  # read depth exists
                                found, readLevel = getReadLevel (minRDs, int(hotspots['region_all'][chrom][start][end]['rd'][c]))
                                if re.match("-", found):
                                    found = "no"
                                rd = getReadLevel (minRDs, int(hotspots['region_all'][chrom][start][end]['rd'][c]))
                        else:  # position is not within design
                            found = "not in design"
                            readLevel = "-"
                            rd = "-"

                        print (sample + "\t" + str(hotspots['region_all'][chrom][start][end]['gene']) + "\t-\t" + str(hotspots['region_all'][chrom][start][end]['aa']) + "\t" + str(hotspots['region_all'][chrom][start][end]['cds']) + "\t" + str(hotspots['region_all'][chrom][start][end]['accNum']) + "\t" + str(hotspots['region_all'][chrom][start][end]['comment']) + "\tregion_all\t" + found + "\t" + readLevel + "\t" + str(hotspots['region_all'][chrom][start][end]['rd'][c]) + "\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t" + str(chrom) + "\t" + str(start) + "\t" + str(end) + "\t-\t-\t-")

# print variants within region, both bwa and pindel. However it does not print total read depth for all positions without variant
def printRegion(hotspots, minRDs, sample, transcripts):
    for chrom in hotspots['region']:
        for start in hotspots['region'][chrom]:
            for end in hotspots['region'][chrom][start]:
                # If there are mutations from bwa print them
                if not re.match("no", str(hotspots['region'][chrom][start][end]['bwa'])):
                    pin = "no"
                    # check if there are also pindel variants present at this position
                    if not re.match("no", str(hotspots['region'][chrom][start][end]['pindel'])):
                        pin = "yes"
                    # Go through all bwa variants
                    for bwaVar in hotspots['region'][chrom][start][end]['bwa']:

                        if re.match("yes", pin):  # If there are pindel variants available extract them
                            for pindelVar in hotspots['region'][chrom][start][end]['pindel']:
                                # Check that there is no pindel variant which is identical to the bwa variant
                                if not (bwaVar[2] == pindelVar[2] and bwaVar[3] == pindelVar[3] and bwaVar[4] == pindelVar[4] and bwaVar[5] == pindelVar[5]):
                                    aa, cds, accNum, exonicType = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information
                                    found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                                    if re.match("-", found):
                                        found = "yes"
                                    print (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + hotspots['region'][chrom][start][end]['comment'] + "\tregion\t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(bwaVar[28]) + "\t" + str(bwaVar[29]) + "\t" + str(bwaVar[26]) + "\t" + str(bwaVar[27]) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(bwaVar[31]) + "\t" + str(bwaVar[30]) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]))
                        else:  # If no pindel variants are available print the variant
                            aa, cds, accNum, exonicType = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information
                            found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                            if re.match("-", found):
                                found = "yes"

                            print (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + hotspots['region'][chrom][start][end]['comment'] + "\tregion\t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(bwaVar[28]) + "\t" + str(bwaVar[29]) + "\t" + str(bwaVar[26]) + "\t" + str(bwaVar[27]) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(bwaVar[31]) + "\t" + str(bwaVar[30]) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]))


                # Check if there are pindel variants available
                if not re.match("no", str(hotspots['region'][chrom][start][end]['pindel'])):

                    for pindelVar in hotspots['region'][chrom][start][end]['pindel']:

                        aa, cds, accNum, exonicType = getTranscriptInfo(pindelVar, transcripts)

                        # Check the level of reads. Below the lowest => -, between first and second => low, above highest => ok
                        found, readLevel = getReadLevel(minRDs, hotspots['region'][chrom][start][end]['rd'][0])
                        if re.match("-", found):
                            found = "yes"

                        print (str(pindelVar[0]) + "\t" + str(pindelVar[6]) + "\t" + exonicType + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + hotspots['region'][chrom][start][end]['comment'] + "\tregion\t" + found + "\t" + readLevel + "\t" + str(pindelVar[12]) + "\t" + str(pindelVar[10]) + "\t" + str(pindelVar[11]) + "\t" + str(pindelVar[9]) + "\t" + str(pindelVar[14]) + "\t" + str(pindelVar[13]) + "\t" + str(pindelVar[16]) + "\t" + str(pindelVar[15]) + "\t" + str(pindelVar[17]) + "\t" + str(pindelVar[18]) + "\t" + str(pindelVar[19]) + "\t" + str(pindelVar[28]) + "\t" + str(pindelVar[29]) + "\t" + str(pindelVar[26]) + "\t" + str(pindelVar[27]) + "\t" + str(pindelVar[20]) + "\t" + str(pindelVar[21]) + "\t" + str(pindelVar[22]) + "\t" + str(pindelVar[23]) + "\t" + str(pindelVar[24]) + "\t" + str(pindelVar[25]) + "\t" + str(pindelVar[31]) + "\t" + str(pindelVar[30]) + "\t" + str(pindelVar[1]) + "\t" + str(pindelVar[2]) + "\t" + str(pindelVar[3]) + "\t" + str(pindelVar[4]) + "\t" + str(pindelVar[5]) + "\t" + str(pindelVar[32]))


# print indels within indel region, both bwa and pindel. However it does not print total read depth for all positions without variant
def printIndel(hotspots, minRDs, sample, transcripts):
    for chrom in hotspots['indel']:
        for start in hotspots['indel'][chrom]:
            for end in hotspots['indel'][chrom][start]:
                # If there are mutations from bwa print them
                if not re.match("no", str(hotspots['indel'][chrom][start][end]['bwa'])):
                    pin = "no"
                    # check if there are also pindel variants present at this position
                    if not re.match("no", str(hotspots['indel'][chrom][start][end]['pindel'])):
                        pin = "yes"
                    # Go through all bwa variants
                    for bwaVar in hotspots['indel'][chrom][start][end]['bwa']:

                        if re.match("yes", pin):  # If there are pindel variants available extract them
                            for pindelVar in hotspots['indel'][chrom][start][end]['pindel']:
                                # Check that there is no pindel variant which is identical to the bwa variant
                                if not (bwaVar[2] == pindelVar[2] and bwaVar[3] == pindelVar[3] and bwaVar[4] == pindelVar[4] and bwaVar[5] == pindelVar[5]):
                                    aa, cds, accNum, exonicType = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information
                                    found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                                    if re.match("-", found):
                                        found = "yes"
                                    print (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + hotspots['indel'][chrom][start][end]['comment'] + "\tindel\t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(bwaVar[28]) + "\t" + str(bwaVar[29]) + "\t" + str(bwaVar[26]) + "\t" + str(bwaVar[27]) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(bwaVar[31]) + "\t" + str(bwaVar[30]) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]))
                        else:  # If no pindel variants are available print the variant
                            aa, cds, accNum, exonicType = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information
                            found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                            if re.match("-", found):
                                found = "yes"

                            print (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + hotspots['indel'][chrom][start][end]['comment'] + "\tindel\t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(bwaVar[28]) + "\t" + str(bwaVar[29]) + "\t" + str(bwaVar[26]) + "\t" + str(bwaVar[27]) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(bwaVar[31]) + "\t" + str(bwaVar[30]) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]))


                # Check if there are pindel variants available
                if not re.match("no", str(hotspots['indel'][chrom][start][end]['pindel'])):

                    for pindelVar in hotspots['indel'][chrom][start][end]['pindel']:

                        aa, cds, accNum, exonicType = getTranscriptInfo(pindelVar, transcripts)

                        # Check the level of reads. Below the lowest => -, between first and second => low, above highest => ok
                        found, readLevel = getReadLevel(minRDs, hotspots['indel'][chrom][start][end]['rd'][0])
                        if re.match("-", found):
                            found = "yes"

                        print (str(pindelVar[0]) + "\t" + str(pindelVar[6]) + "\t" + exonicType + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + hotspots['indel'][chrom][start][end]['comment'] + "\tindel\t" + found + "\t" + readLevel + "\t" + str(pindelVar[12]) + "\t" + str(pindelVar[10]) + "\t" + str(pindelVar[11]) + "\t" + str(pindelVar[9]) + "\t" + str(pindelVar[14]) + "\t" + str(pindelVar[13]) + "\t" + str(pindelVar[16]) + "\t" + str(pindelVar[15]) + "\t" + str(pindelVar[17]) + "\t" + str(pindelVar[18]) + "\t" + str(pindelVar[19]) + "\t" + str(pindelVar[28]) + "\t" + str(pindelVar[29]) + "\t" + str(pindelVar[26]) + "\t" + str(pindelVar[27]) + "\t" + str(pindelVar[20]) + "\t" + str(pindelVar[21]) + "\t" + str(pindelVar[22]) + "\t" + str(pindelVar[23]) + "\t" + str(pindelVar[24]) + "\t" + str(pindelVar[25]) + "\t" + str(pindelVar[31]) + "\t" + str(pindelVar[30]) + "\t" + str(pindelVar[1]) + "\t" + str(pindelVar[2]) + "\t" + str(pindelVar[3]) + "\t" + str(pindelVar[4]) + "\t" + str(pindelVar[5]) + "\t" + str(pindelVar[32]))



