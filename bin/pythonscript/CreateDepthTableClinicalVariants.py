import argparse
import re
from os import walk

from collections import OrderedDict

# Parse commandline
parser = argparse.ArgumentParser(description = "")
parser.add_argument('-f', '--file',
                    help = 'File(s) that will be processed',
                    action='append',
                    default = [],
                    required = True)
parser.add_argument('-c', '--codons',
                    help = 'File with codons that should be extracted (Gene name\tcondon nr',
                    type = str,
                    required = False)

args = parser.parse_args()
# outfile = open(args.output, 'w')

# Go through all files in the list

data = {}

def columnsToDict(headerMap, columns):
    return {'#Gene': columns[headerMap['#Gene']],
            'Cosmic_id': columns[headerMap['Cosmic_id']],
            'Found': columns[headerMap['Found']],
            'Min_read_depth30': columns[headerMap['Min_read_depth30']],
            'Min_read_depth100': columns[headerMap['Min_read_depth100']],
            'Min_read_depth300': columns[headerMap['Min_read_depth300']],
            'Min_read_depth500': columns[headerMap['Min_read_depth500']],
            'Total_read_depth': columns[headerMap['Total_read_depth']],
            'Chr': columns[headerMap['Chr']],
            'Start': columns[headerMap['Start']],
            'End': columns[headerMap['End']],
            'Reference': columns[headerMap['Reference']],
            'Variant': columns[headerMap['Variant']],
            'Variant_read_depth': columns[headerMap['Variant_read_depth']],
            'Reference_read_depth': columns[headerMap['Reference_read_depth']],
            'Variant_allele_ratio': columns[headerMap['Variant_allele_ratio']],
            'CDS_change': columns[headerMap['CDS_change']],
            'AA_change': columns[headerMap['AA_change']],
            'Reference_plus_amplicons': columns[headerMap['Reference_plus_amplicons']],
            'Reference_minus_amplicons': columns[headerMap['Reference_minus_amplicons']],
            'Variant_plus_amplicons': columns[headerMap['Variant_plus_amplicons']],
            'Variant_minus_amplicons': columns[headerMap['Variant_minus_amplicons']],
            'Transcript': columns[headerMap['Transcript']],
            'Ref_amplicons': columns[headerMap['Ref_amplicons']],
            'Var_amplicons': columns[headerMap['Var_amplicons']]}

codons = {}
if args.codons is not None:
    with open(args.codons, 'r') as c:
        for line in c:
            line = line.rstrip()
            columns = line.split('\t')
            if columns[0] in codons:
                codons[columns[0]].append(columns[1])
            else:
                codons[columns[0]] = [columns[1]]

for f in args.file:
    sample = re.search('\/*([A-Za-z0-9]+)\.h\.sapiens\.clinicalPositions\.txt$',f).group(1)
    data[sample] = {}
    with open(f, 'r') as f:
        headerMap = {}
        for line in f:
            line = line.rstrip()
            columns = line.split('\t')
            if line.startswith('#Gene'):
                for i, c in enumerate(columns):
                    headerMap[c] = i
            else:
                codon = ""
                if columns[headerMap['AA_change']] == 'p.?':
                    codon = "Unknown"
                else:
                    codon = re.search('^[a-zA-Z.*]*([0-9]+|[0-9]+_[A-Z*]*[0-9]+)[0-9A-Za-z*>]*$',columns[headerMap['AA_change']]).group(1)
                if (not codons or
                       (columns[headerMap['#Gene']] in codons and
                        codon in codons[columns[headerMap['#Gene']]])):
                    if columns[headerMap['#Gene']] in data[sample]:
                        if codon in data[sample][columns[headerMap['#Gene']]]:
                            data[sample][columns[headerMap['#Gene']]][codon].append(columnsToDict(headerMap,columns))
                        else:
                            data[sample][columns[headerMap['#Gene']]][codon] = [columnsToDict(headerMap,columns)]
                    else:
                        data[sample][columns[headerMap['#Gene']]] = {codon: [columnsToDict(headerMap,columns)]}

pattern = re.compile('^.*(del|ins).*$')

printData = OrderedDict()

for sample in data:
    for gene in data[sample]:
        for codon in data[sample][gene]:
            row_info = ""
            depth = 0
            counter = 0
            if len(data[sample][gene][codon]) > 1:
                for row in data[sample][gene][codon]:
                    row_info = gene + "\t" + codon + "\t"
                    if not pattern.match(row['AA_change']) and not row['Total_read_depth'] == "-" :
                        counter = counter + 1
                        depth = depth + int(row['Total_read_depth'])
            else:
                row_info = data[sample][gene][codon][0]['#Gene'] + "\t" + codon + "\t"
                depth = int(data[sample][gene][codon][0]['Total_read_depth'])
                counter = 1
            if row_info in  printData:
                printData[row_info][sample] = str(depth/counter)
            else:
                printData[row_info] = OrderedDict()
                printData[row_info][sample] = str(depth / counter)

import sys
sys.stdout.write("#Gere\tcodon\t")
for sample in data:
    sys.stdout.write("\t" + sample)
sys.stdout.write("\n")

for codon in printData:
    sys.stdout.write(codon)
    for sample in data:
        sys.stdout.write("\t" + printData[codon][sample])
    sys.stdout.write("\n")
