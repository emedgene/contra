#!/usr/bin/python
# Usage: python 2_process.py [_renorm.txt] [_wg.txt] [_exon.txt] [bed.annotated.txt] [debug]
import copy
import sys

#-----------------------------------------------------------------------------------
# Based on the 27/08 version (Ca2015_Histograms/renorm.py)
# Modified to just be run on one file
# More descriptive for the exon level (coordinates in columns)
# (And eventually, command line input!)
# 17/09 - removed "Name" column
#       - edited _exon.txt output to use matched/"right" coordinates
# 28/10 - no longer outputs NA
#-----------------------------------------------------------------------------------

# Global dictionary - {<Gene>: {<Exon> : [List of log ratios]}}
# where Exon = startcoord_endcoord (Should be chromosome_startcoord_endcoord but thats not in renorm)
gdict = {}

if len(sys.argv) < 4:
    print("2_process.py: not enough input arguments")
    sys.exit(1)

filename = sys.argv[1]
wholegene_output = sys.argv[2]
exon_output = sys.argv[3]
debug='F'
if len(sys.argv) == 6:
    debug=sys.argv[5]


# Create gdict entries based on myTSCA.bed.annotated.txt
# Format: chromosome startcoord endcoord gene
if len(sys.argv) >= 5:
    myTSCA_in = sys.argv[4]
else:
    myTSCA_in = "/mnt/Storage/NextGenSequencing/Data/Yeap/myTSCA.bed.annotated.txt"

myTSCA = open(myTSCA_in, 'r')
for line in myTSCA.readlines():
    entries = line.split('\t')
    entries[3] = entries[3].rstrip('\n')
    if entries[3] not in gdict:
        gdict[entries[3]] = {}
#    gdict[entries[3]]['_'.join(entries[1:3])] = []
    gdict[entries[3]][entries[0]+":"+'_'.join(entries[1:3])] = []
myTSCA.close()

# Helper function
def calculate_average(numblist):
    total = 0.0
    for item in numblist:
        total += float(item)
    return total / float(len(numblist))

# Helper function - list may contain "NA" entries
def calculate_average_withmissing(numblist):
    total = 0.0
    count = 0
    for item in numblist:
        if item != "NA":
            total += float(item)
            count += 1
    if count > 0:
        return total / float(count)
    else:
        return "NA"

# Given a CONTRA output file, generate a list of exons with average log ratio per exon.
def process_file(f_in):
    GENE_IND=0
    CHR_IND=1
    START_IND=2
    END_IND=3
    LR_IND=4

    f = open(f_in, 'r')
    # 0 - Gene Sym
    # 1 - OriStCoord, 2 - OriEndCoord
    # 3 - Mean.of.LogRatio (renormalised)

    # Exon Data Format:
    # XXX<GeneName>: { <ExonNumber>: [StartCoord, EndCoord, MeanLogRatio] }
    # <GeneName>: { <ExonNumber>: [StartCoord, EndCoord, MeanLogRatio,Chr] }
    exons = {}

    # Temp list of all genes in gdict
    all_genes = list(gdict.keys())
    
    # Preprocessing
    contra = []
    for line in f.readlines()[1:]:
        newline = line.split('\t')
        name = newline[GENE_IND].split('(')[0]
        newline[GENE_IND] = name
        newline[-1] = newline[-1].rstrip('\n')
        contra.append(newline)

    f.close()

    currgene = contra[0][GENE_IND]
    temp = [contra[0][START_IND], contra[0][END_IND], [contra[0][LR_IND]], contra[0][CHR_IND]] # start end list of averages
    exons[currgene] = {}
    exon_n = 1
    for line in contra[1:]:
        # If its a new gene...
        if line[GENE_IND] != currgene:
            # Add old to exons dictionary, after calculating average log ratio over the exon.
            temp[2] = calculate_average(temp[2])
            exons[currgene][exon_n] = temp
            currgene = line[GENE_IND]
            temp = [line[START_IND], line[END_IND], [line[LR_IND]], line[CHR_IND]]
            if currgene not in exons: # if exon data gets broken up...
                exon_n = 1
                exons[currgene] = {}
            else:
                exon_n = len(list(exons[currgene].keys())) + 1
                
        # If its a new exon...
        elif line[START_IND] != temp[1]: # ie coordinates don't match up
            # Add old to exons dictionary, after calculating average log ratio over the exon.
            temp[2] = calculate_average(temp[2])
            exons[currgene][exon_n] = temp
            exon_n += 1
            temp = [line[START_IND], line[END_IND], [line[LR_IND]], line[CHR_IND]]
        else:
            # Part of the current exon. change the current end coord and add average to list of averages.
            temp[2].append(line[LR_IND])
            temp[1] = line[END_IND]
    # Finish off current gene - BUGFIX
    temp[2] = calculate_average(temp[2])
    exons[currgene][exon_n] = temp

    # Add appropriate data to global dictionary gdict
    # If gdict not yet populated (empty):
    for gene in list(exons.keys()):
        exon_list = list(gdict[gene].keys())
        for exon in list(exons[gene].keys()):
            start = int(exons[gene][exon][0])
            end = int(exons[gene][exon][1])
            ch = exons[gene][exon][3].replace("chr","")
            exon_name = ch+":"+str(start)+'_'+str(end)
            if exon_name in exon_list:
                gdict[gene][exon_name].append(exons[gene][exon][2])
                exon_list.remove(exon_name)
            else:
                # Attempt to match it to something in exon_list: up to 20 bases missing from each side is OK
                matched = False
                for possible_match in exon_list:
                    components = possible_match.split(":")[1].split('_')
                    possible_start = int(components[0])
                    possible_end = int(components[1])
                    start_diff = start - possible_start
                    end_diff = possible_end - end
                    if 0 <= start_diff < 21 and 0 <= end_diff < 21:
                        # Can match these two together
                        gdict[gene][possible_match].append(exons[gene][exon][2])
                        exon_list.remove(possible_match)
                        matched = True
                    # Couldn't find a match for it: output error
        # Take care of missing data: remaining exons in exonlist

        for missing_exon in exon_list:
            gdict[gene][missing_exon].append("NA")
        all_genes.remove(gene)

    for gene in all_genes:
        # These genes had nothing in the input data. So, add NA to each exon of the gene.
        for exon in gdict[gene]:
            gdict[gene][exon].append("NA")


process_file(filename)


f2 = open(exon_output, 'w')
headers = ['Gene','Chr', 'Start Coord', 'End Coord', 'Mean Log Ratio\n']
f2.write('\t'.join(headers))
for gene in sorted(gdict.keys()):
    for exon in sorted(gdict[gene].keys()):
        log_ratio = str(gdict[gene][exon][0])
        if (log_ratio != "NA"):
            ch,coord=exon.split(":")
            out = [gene, ch, coord.split('_')[0], coord.split('_')[1], str(gdict[gene][exon][0])]
            f2.write('\t'.join(out))
            f2.write('\n')
f2.close()


# Duplicate gdict, just in case...

gdict2 = copy.deepcopy(gdict)

for gene in list(gdict2.keys()):
    if ';' in gene: # If gene needs to be split
        components = gene.split(';')
        for component_gene in components:
            if component_gene not in gdict2:
                gdict2[component_gene] = {}
        for exon in gdict[gene]:
            gdict2[components[0]][exon] = gdict[gene][exon][:]
            gdict2[components[1]][exon] = gdict[gene][exon][:]
        del gdict2[gene]

wholegene = open(wholegene_output, 'w')
# write headers
wholegene.write("Gene\tRatio\n")
for gene in sorted(gdict2.keys()):
    sample_temp = []
    for exon in list(gdict2[gene].keys()): # For each exon (in that gene):
        sample_temp.append(gdict2[gene][exon][0]) # collect the log ratio of that exon of that sample (of that gene)
    log_ratio = calculate_average_withmissing(sample_temp)
    if log_ratio != "NA":
        wholegene.write(gene)
        wholegene.write('\t'+str(log_ratio)) # whole gene log ratio = average of log ratios across all exons of that gene
        wholegene.write('\n')
wholegene.close()


if debug == 'T':
    print(("2_process.py: created "+wholegene_output))
    print(("2_process.py: created "+exon_output))


