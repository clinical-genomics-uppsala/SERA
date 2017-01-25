import re

def createHotspotHash(lineSplit, hotspots, intronic):
    chrom = lineSplit[0]
    start = int(lineSplit[1])
    end = int(lineSplit[2])
    gene = lineSplit[3]
    cds = lineSplit[4]
    aa = lineSplit[5]
    report = lineSplit[6]
    comment = lineSplit[7]
    exon = lineSplit[8]
    accNum = lineSplit[9].split(".")[0]
    varInfo = "no"

    if not report in hotspots:
        hotspots[report] = {}
        hotspots[report][chrom] = {}
        hotspots[report][chrom][start] = {}
        hotspots[report][chrom][start][end] = {}
        hotspots[report][chrom][start][end]['gene'] = gene
        hotspots[report][chrom][start][end]['cds'] = cds
        hotspots[report][chrom][start][end]['aa'] = aa
        hotspots[report][chrom][start][end]['comment'] = comment
        hotspots[report][chrom][start][end]['exon'] = exon
        hotspots[report][chrom][start][end]['accNum'] = accNum
        hotspots[report][chrom][start][end]['bwa'] = varInfo
        hotspots[report][chrom][start][end]['pindel'] = varInfo
        hotspots[report][chrom][start][end]['rd'] = "-"

    elif not chrom in hotspots[report]:
        hotspots[report][chrom] = {}
        hotspots[report][chrom][start] = {}
        hotspots[report][chrom][start][end] = {}
        hotspots[report][chrom][start][end]['gene'] = gene
        hotspots[report][chrom][start][end]['cds'] = cds
        hotspots[report][chrom][start][end]['aa'] = aa
        hotspots[report][chrom][start][end]['comment'] = comment
        hotspots[report][chrom][start][end]['exon'] = exon
        hotspots[report][chrom][start][end]['accNum'] = accNum
        hotspots[report][chrom][start][end]['bwa'] = varInfo
        hotspots[report][chrom][start][end]['pindel'] = varInfo
        hotspots[report][chrom][start][end]['rd'] = "-"

    elif not start in hotspots[report][chrom]:
        hotspots[report][chrom][start] = {}
        hotspots[report][chrom][start][end] = {}
        hotspots[report][chrom][start][end]['gene'] = gene
        hotspots[report][chrom][start][end]['cds'] = cds
        hotspots[report][chrom][start][end]['aa'] = aa
        hotspots[report][chrom][start][end]['comment'] = comment
        hotspots[report][chrom][start][end]['exon'] = exon
        hotspots[report][chrom][start][end]['accNum'] = accNum
        hotspots[report][chrom][start][end]['bwa'] = varInfo
        hotspots[report][chrom][start][end]['pindel'] = varInfo
        hotspots[report][chrom][start][end]['rd'] = "-"

    elif not end in hotspots[report][chrom][start]:
        hotspots[report][chrom][start][end] = {}
        hotspots[report][chrom][start][end]['gene'] = gene
        hotspots[report][chrom][start][end]['cds'] = cds
        hotspots[report][chrom][start][end]['aa'] = aa
        hotspots[report][chrom][start][end]['comment'] = comment
        hotspots[report][chrom][start][end]['exon'] = exon
        hotspots[report][chrom][start][end]['accNum'] = accNum
        hotspots[report][chrom][start][end]['bwa'] = varInfo
        hotspots[report][chrom][start][end]['pindel'] = varInfo
        hotspots[report][chrom][start][end]['rd'] = "-"

    else:
        print ("Position already exists: " + chrom + " " + str(start) + " " + str(end))

    # If report is region_all add list for rdata_address
    if re.match("hotspot", report) or re.match("region_all", report):
        hotspots[report][chrom][start][end]['rd'] = ["-"] * (end - start + 1)

    # If this is intronic region add to intronic hash
    if re.match("intronic", exon):
        if not chrom in intronic:
            intronic[chrom] = {}
            intronic[chrom][start] = {}
            intronic[chrom][start][end] = "intronic"
        elif not start in intronic:
            intronic[chrom][start] = {}
            intronic[chrom][start][end] = "intronic"


# Check if the variant is a hotspot
def hotspotVariant(lineSplit, hotspots):
    chrom = lineSplit[1]
    start = int(lineSplit[2])
    end = int(lineSplit[3])

    # This is based on that an hotspot in the hotspot file is only 1bp (the variant can differ between one and more bp)
    if "hotspot" in hotspots:
        # start with checking if the chromosome exists in the hotspot hash
        if chrom in hotspots['hotspot']:
            # If the variant is a one bp variant check if that position exists in the hotspot hash
            if end - start == 0:
                if start in hotspots['hotspot'][chrom]:
                    if not ("+" in lineSplit[24]) and not ("+" in lineSplit[25]):  # if bwa results
                        if re.match("no", hotspots['hotspot'][chrom][start][start]['bwa']):  # If no bwa variant is added before
                            hotspots['hotspot'][chrom][start][end]['bwa'] = []  # Create list for the variant

                        hotspots['hotspot'][chrom][start][end]['bwa'].append(lineSplit)  # Add the variant to the list
                        return True  # return True that the variant is added
                    else:  # Pindel result
                        if re.match("no", hotspots['hotspot'][chrom][start][start]['pindel']):  # If no pindel variant was added before
                            hotspots['hotspot'][chrom][start][end]['pindel'] = []  # Create list for the variant

                        hotspots['hotspot'][chrom][start][end]['pindel'].append(lineSplit)  # Add the variant to the list
                        return True  # return True that the variant is added
            # If it's a multiple bp variant check if the start position exists in the hash
            else:
                if not ("+" in lineSplit[24]) and not ("+" in lineSplit[25]):  # bwa results
                    # start of the variant exist in hotspot hash add variant to start position
                    if start in hotspots['hotspot'][chrom]:
                        if re.match("no", hotspots['hotspot'][chrom][start][start]['bwa']):
                            hotspots['hotspot'][chrom][start][start]['bwa'] = []

                        hotspots['hotspot'][chrom][start][start]['bwa'].append(lineSplit)
                        return True
                    # If the start position is not within the hash check if the variant overlaps a hotspot position
                    else:
                        for s in hotspots['hotspot'][chrom]:
                            for e in hotspots['hotspot'][chrom][s]:
                                if start <= e and end > s:
                                    if re.match("no", hotspots['hotspot'][chrom][s][e]['bwa']):
                                        hotspots['hotspot'][chrom][s][e]['bwa'] = []
                                    hotspots['hotspot'][chrom][s][e]['bwa'].append(lineSplit)
                                    return True
                else:  # pindel result
                    if start in hotspots['hotspot'][chrom]:
                        if re.match("no", hotspots['hotspot'][chrom][start][start]['pindel']):
                            hotspots['hotspot'][chrom][start][start]['pindel'] = []
                        hotspots['hotspot'][chrom][start][start]['pindel'].append(lineSplit)
                        return True
                    # If the start position is not within the hash check if the variant overlaps a hotspot position
                    else:
                        for s in hotspots['hotspot'][chrom]:
                            for e in hotspots['hotspot'][chrom][s]:
                                if start <= e and end > s:
                                    if re.match("no", hotspots['hotspot'][chrom][s][e]['pindel']):
                                        hotspots['hotspot'][chrom][s][e]['pindel'] = []
                                    hotspots['hotspot'][chrom][s][e]['pindel'].append(lineSplit)
                                    return True

    return False

def regionVariant(lineSplit, hotspots):
    chrom = lineSplit[1]
    start = lineSplit[2]
    end = lineSplit[3]

    if "region_all" in hotspots:
        # start with checking if the chromosome exists in the region_all part of the hash
        if chrom in hotspots['region_all']:
            if not ("+" in lineSplit[24]) and not ("+" in lineSplit[25]):  # bwa results
                if start in hotspots['region_all'][chrom]:
                    added = False
                    for e in hotspots['region_all'][chrom][start]:
                        if not added:  # if the variant is not added already
                            if re.match("no", hotspots['region_all'][chrom][start][e]['bwa']):
                                hotspots['region_all'][chrom][start][e]['bwa'] = []
                            hotspots['region_all'][chrom][start][e]['bwa'].append(lineSplit)
                            added = True  # set the variant to added
                            return True

                         # If the start position is not within the hash check if the variant overlaps a hotspot position
                else:
                    added = False
                    for s in hotspots['region_all'][chrom]:
                        for e in hotspots['region_all'][chrom][s]:
                            if start <= e and end > s:  # if the variant overlaps at least 1 bp (start position excluded since that is handled above)
                                if re.match("no", hotspots['region_all'][chrom][s][e]['bwa']):
                                    hotspots['region_all'][chrom][s][e]['bwa'] = []
                                hotspots['region_all'][chrom][s][e]['bwa'].append(lineSplit)
                                added = True
                                return True

            else:  # pindel results
                if start in hotspots['region_all'][chrom]:
                    added = False
                    for e in hotspots['region_all'][chrom][start]:
                        if not added:  # if the variant is not added already
                            if re.match("no", hotspots['region_all'][chrom][start][e]['pindel']):
                                hotspots['region_all'][chrom][start][e]['pindel'] = []
                            hotspots['region_all'][chrom][start][e]['pindel'].append(lineSplit)
                            added = True  # set the variant to added
                            return True

                         # If the start position is not within the hash check if the variant overlaps a hotspot position
                else:
                    added = False
                    for s in hotspots['region_all'][chrom]:
                        for e in hotspots['region_all'][chrom][s]:
                            if start <= e and end > s and not added:  # if the variant overlaps at least 1 bp (start position excluded since that is handled above)
                                if re.match("no", hotspots['region_all'][chrom][s][e]['pindel']):
                                    hotspots['region_all'][chrom][s][e]['pindel'] = []
                                hotspots['region_all'][chrom][s][e]['pindel'].append(lineSplit)
                                added = True
                                return True

    if "region" in hotspots:
        # start with checking if the chromosome exists in the region part of the hash
        if chrom in hotspots['region']:
            if not ("+" in lineSplit[24]) and not ("+" in lineSplit[25]):  # bwa results
                if start in hotspots['region'][chrom]:
                    added = False
                    for e in hotspots['region'][chrom][start]:
                        if not added:  # if the variant is not added already
                            if re.match("no", hotspots['region'][chrom][start][e]['bwa']):
                                hotspots['region'][chrom][start][e]['bwa'] = []
                            hotspots['region'][chrom][start][e]['bwa'].append(lineSplit)
                            added = True  # set the variant to added
                            return True

                         # If the start position is not within the hash check if the variant overlaps a hotspot position
                else:
                    added = False
                    for s in hotspots['region'][chrom]:
                        for e in hotspots['region'][chrom][s]:
                            if start <= e and end > s:  # if the variant overlaps at least 1 bp (start position excluded since that is handled above)
                                if re.match("no", hotspots['region'][chrom][s][e]['bwa']):
                                    hotspots['region'][chrom][s][e]['bwa'] = []
                                hotspots['region'][chrom][s][e]['bwa'].append(lineSplit)
                                added = True
                                return True

            else:  # pindel results
                if start in hotspots['region'][chrom]:
                    added = False
                    for e in hotspots['region'][chrom][start]:
                        if not added:  # if the variant is not added already
                            if re.match("no", hotspots['region'][chrom][start][e]['pindel']):
                                hotspots['region'][chrom][start][e]['pindel'] = []
                            hotspots['region'][chrom][start][e]['pindel'].append(lineSplit)
                            added = True  # set the variant to added
                            return True

                         # If the start position is not within the hash check if the variant overlaps a hotspot position
                else:
                    added = False
                    for s in hotspots['region'][chrom]:
                        for e in hotspots['region'][chrom][s]:
                            if start <= e and end > s and not added:  # if the variant overlaps at least 1 bp (start position excluded since that is handled above)
                                if re.match("no", hotspots['region'][chrom][s][e]['pindel']):
                                    hotspots['region'][chrom][s][e]['pindel'] = []
                                hotspots['region'][chrom][s][e]['pindel'].append(lineSplit)
                                added = True
                                return True

    return False


def indelVariant(lineSplit, hotspots):
    chrom = lineSplit[1]
    start = lineSplit[2]
    end = lineSplit[3]

    if "indel" in hotspots:
        # start with checking if the chromosome exists in the indel part of the hash
        if chrom in hotspots['indel']:
            if not ("+" in lineSplit[24]) and not ("+" in lineSplit[25]):  # bwa results
                if start in hotspots['indel'][chrom]:
                    added = False
                    for e in hotspots['indel'][chrom][start]:
                        if not added:  # if the variant is not added already
                            if re.match("no", hotspots['indel'][chrom][start][e]['bwa']):
                                hotspots['indel'][chrom][start][e]['bwa'] = []
                            hotspots['indel'][chrom][start][e]['bwa'].append(lineSplit)
                            added = True  # set the variant to added
                            return True

                         # If the start position is not within the hash check if the variant overlaps a hotspot position
                else:
                    added = False
                    for s in hotspots['indel'][chrom]:
                        for e in hotspots['indel'][chrom][s]:
                            if start <= e and end > s:  # if the variant overlaps at least 1 bp (start position excluded since that is handled above)
                                if re.match("no", hotspots['indel'][chrom][s][e]['bwa']):
                                    hotspots['indel'][chrom][s][e]['bwa'] = []
                                hotspots['indel'][chrom][s][e]['bwa'].append(lineSplit)
                                added = True
                                return True

            else:  # pindel results
                if start in hotspots['indel'][chrom]:
                    added = False
                    for e in hotspots['indel'][chrom][start]:
                        if not added:  # if the variant is not added already
                            if re.match("no", hotspots['indel'][chrom][start][e]['pindel']):
                                hotspots['indel'][chrom][start][e]['pindel'] = []
                            hotspots['indel'][chrom][start][e]['pindel'].append(lineSplit)
                            added = True  # set the variant to added
                            return True

                         # If the start position is not within the hash check if the variant overlaps a hotspot position
                else:
                    added = False
                    for s in hotspots['indel'][chrom]:
                        for e in hotspots['indel'][chrom][s]:
                            if start <= e and end > s and not added:  # if the variant overlaps at least 1 bp (start position excluded since that is handled above)
                                if re.match("no", hotspots['indel'][chrom][s][e]['pindel']):
                                    hotspots['indel'][chrom][s][e]['pindel'] = []
                                hotspots['indel'][chrom][s][e]['pindel'].append(lineSplit)
                                added = True
                                return True

    return False

def filterAnnovar(lineSplit, minRD, blacklist, g1000, ampliconmapped, intronic):
    chrom = lineSplit[1]
    start = lineSplit[2]
    end = lineSplit[3]
    ref = lineSplit[4]
    var = lineSplit[5]
    gene = lineSplit[6]

    intr = False  # This is not an intronic position to save
    if chrom in intronic:
        for c in intronic:
            for s in intronic[c]:
                for e in intronic[c][s]:
                    if start <= e and end >= s:  # Check if the variant overlaps an intronic region to keep
                        intr = True

    # Keep if it's exonic or splicing or intronic to keep
    if re.match('^exonic', lineSplit[7]) or re.match('splicing', lineSplit[7]) or intr:

        # Remove all synonymous variants
        if (not re.match('^synonymous', lineSplit[8])) or intr:

            # Keep those without a value in 1000G or a value that is lower than the given input
            if re.match('-', str(lineSplit[13])) or float(lineSplit[13]) <= g1000:

                varOK = True
                # Check if ampliconmapped is set and if so check amplicon status
                if ampliconmapped:
                    varOK = False
                    # Set amplicon info

                    ref_plus = ref_minus = var_plus = var_minus = 0
                    if not re.match('-', lineSplit[26]):
                        var_plus = int(lineSplit[26])
                    if not re.match('-', lineSplit[27]):
                        var_minus = int(lineSplit[27])
                    if not re.match('-', lineSplit[28]):
                        ref_plus = int(lineSplit[28])
                    if not re.match('-', lineSplit[29]):
                        ref_minus = int(lineSplit[29])

                    if (ref_plus + ref_minus) >= 3:
                        if (var_plus + var_minus) >= 2:
                            varOK = True
                    elif (ref_plus + ref_minus) >= 1:
                        if (var_plus + var_minus) >= 1:
                            varOK = True
                    elif (ref_plus + ref_minus) >= 0:
                        if (var_plus + var_minus) >= 0:
                            varOK = True

                # Check if the amplicon info is okey
                if varOK:
                    # Remove blacklist variants
                    if not chrom in blacklist:
                        return True
                    else:
                        if not start in blacklist[chrom]:
                            return True
                        else:
                            if not end in blacklist[chrom][start]:
                                return True
                            else:
                                if not ref in blacklist[chrom][start][end]:
                                    return True
                                else:
                                    if not var in blacklist[chrom][start][end][ref]:
                                        return True
                                    else:
                                        return False

    return False

def createBlacklist(line, blacklist):
    # Remove new line characters and split string on tab
    line = line.rstrip('\r\n').split("\t")

    chrom = line[0]
    start = line[1]
    end = line[2]
    ref = line[3]
    var = line[4]
    gene = line[5]

    # Create hash with blacklist variants
    if not chrom in blacklist:
        blacklist[chrom] = {}
        blacklist[chrom][start] = {}
        blacklist[chrom][start][end] = {}
        blacklist[chrom][start][end][ref] = {}
        blacklist[chrom][start][end][ref][var] = gene
    elif not start in blacklist[chrom]:
        blacklist[chrom][start] = {}
        blacklist[chrom][start][end] = {}
        blacklist[chrom][start][end][ref] = {}
        blacklist[chrom][start][end][ref][var] = gene
    elif not end in blacklist[chrom][start]:
        blacklist[chrom][start][end] = {}
        blacklist[chrom][start][end][ref] = {}
        blacklist[chrom][start][end][ref][var] = gene
    elif not ref in blacklist[chrom][start][end]:
        blacklist[chrom][start][end][ref] = {}
        blacklist[chrom][start][end][ref][var] = gene
    elif not var in blacklist[chrom][start][end][ref]:
        blacklist[chrom][start][end][ref][var] = gene

# Add info about a variant existing in the hotspot hash
def addVariantInfo(lineSplit, minRD, blacklist, g1000, ampliconmapped, hotspots, intronic):
    chrom = lineSplit[1]
    start = int(lineSplit[2])
    ref = lineSplit[4]
    var = lineSplit[5]

    # If the variant is okey add to hash
    varOK = filterAnnovar(lineSplit, minRD, blacklist, g1000, ampliconmapped, intronic)

    variantAdded = False
    if varOK:

        variantAdded = hotspotVariant(lineSplit, hotspots)
        if not variantAdded:
            variantAdded = regionVariant(lineSplit, hotspots)

        if not variantAdded:
            if re.match("-", ref) or re.match("-", var):  # Add indel info
                variantAdded = indelVariant(lineSplit, hotspots)

    else:
        variantAdded = True  # Variant should not be added set to true
    return variantAdded

# Method to return the read level of the position
def getReadLevel(minRDs, rd):
    # Check the level of reads. Below the lowest => -, between first and second => low, above highest => ok
    readLevel = "ok"
    found = "-"
    # If the read depth is lower than minRDs[0] report found as not analyzable
    if rd < int(minRDs[0]):
        found = "not analyzable"
    # If the read depth is lower than the higher rd report rd as low
    if rd < int(minRDs[1]):
        readLevel = "low"

    return found, readLevel

# Get the aa, cds and accessionnumber information extracted
def getTranscriptInfo(lineSplit, transcripts):
    exonicType = lineSplit[8]
    gene = lineSplit[6]
    if re.match("splicing", lineSplit[7]):
        exonicType = lineSplit[7]
    aa = cds = accNum = exon = "-"  # set the parameters to - as default
    comm = "-"  # If variable to add comment in
    found = False;
    if not re.match("-", lineSplit[32]):  # Check that there are any transcript
        allTranscript = lineSplit[32].split(",")  # Split all transcripts on ,
        # If the gene exist in preferd transcript hash use that transcript, otherwise take the first
        if gene in transcripts:  # check if the genename exist in the transcript hash
            for tr in allTranscript:  # if so go through and see if any transcript overlaps with the main transcript in hash
                if transcripts[gene] in tr:
                    transcriptInfo = tr.split(":")
                    found = True  # The main transript is found in the list of transcript for this mutation
                    if len(transcriptInfo) >= 5:
                        aa = transcriptInfo[4]
                        cds = transcriptInfo[3]
                        exon = transcriptInfo[2]
                        accNum = transcriptInfo[1]
                    elif len(transcriptInfo) >= 4:
                        cds = transcriptInfo[3]
                        exon = transcriptInfo[2]
                        accNum = transcriptInfo[1]
                    elif len(transcriptInfo) >= 3:
                        exon = transcriptInfo[2]
                        accNum = transcriptInfo[1]
                    elif len(transcriptInfo) >= 2:
                        accNum = transcriptInfo[1]

            if not found:
                comm = "altTranscript"
        else:
            transcriptInfo = allTranscript[0].split(":")

            if len(transcriptInfo) >= 5:
                aa = transcriptInfo[4]
                cds = transcriptInfo[3]
                exon = transcriptInfo[2]
                accNum = transcriptInfo[1]
            elif len(transcriptInfo) >= 4:
                cds = transcriptInfo[3]
                exon = transcriptInfo[2]
                accNum = transcriptInfo[1]
            elif len(transcriptInfo) >= 3:
                exon = transcriptInfo[2]
                accNum = transcriptInfo[1]
            elif len(transcriptInfo) >= 2:
                accNum = transcriptInfo[1]

    return aa, cds, accNum, exon, exonicType, comm

# Give the comment to report
def setComment(transcriptComm, inputComm):
    comm = "-"
    if not re.match("-", transcriptComm):
        if not re.match("-", inputComm):
            comm = transcriptComm + " " + inputComm
        else:
            comm = transcriptComm
    else:
        comm = inputComm



