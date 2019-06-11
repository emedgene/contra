# ----------------------------------------------------------------------#
# Copyright (c) 2019, Richard Lupat, Jason Li, Adam Lev-Libfeld
#
# > Source License <
# This file is part of CONTRA.
#
#    CONTRA is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    CONTRA is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with CONTRA.  If not, see <http://www.gnu.org/licenses/>.
#
# 
#-----------------------------------------------------------------------#
# Last Updated : 2019-04-17 15:59

import os

def splitByChromosome(destFolder, inputfile, skip_headers=True):
    try:
        os.mkdir(destFolder + "chr/")
    except:
        print("folder exist")

    outputfile = destFolder + "chr/chr1.txt"
    file = open(inputfile,"r")
    output = open(outputfile,"w")
    check = "1"

    for row in file:
        if skip_headers and row[0] == '#':
            continue

        cols = row.split()
        chr = cols[0].strip("chr")
        if (chr != check):
            output.close()
            check = chr
            output = open(destFolder+ "chr/chr"+check+".txt","w")
        output.write(row)

    output.close()

def splitFolderByChromosome(destFolder):
    infile = destFolder + "sample.BEDGRAPH"
    splitByChromosome(destFolder, infile, skip_headers=True)


def splitFileByChromosome(infile):
    destFolder = os.path.dirname(infile)+"/"
    splitByChromosome(destFolder, infile, skip_headers=False)