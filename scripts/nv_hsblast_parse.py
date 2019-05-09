"""
This script parses the output from HS-BLASTN.

Copyright (C) 2019 Tham Cheng Yong, Roberto Tirado Magallanes, Touati Benoukraf.

This file is part of NanoVar.

NanoVar is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

NanoVar is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with NanoVar.  If not, see <https://www.gnu.org/licenses/>.
"""

from sys import argv
import sys
import os
import logging

#Ensure correct number of inputs
if len(argv)<4 or len(argv)>=5:
    sys.exit("Usage: python nv_hsblast_parse.py file.hsblast file.fa log_file > output.tsv")

#Assign variable to inputs
file1 = argv[1]
file2 = argv[2]
log_path = argv[3]

#Calculating query lengths
def measureQlen(fasta):
    global qlendict
    qlendict = {}
    with open(fasta) as f:
        line1 = [next(f) for x in xrange(2)]
        if line1[0][0] == '>': #fasta
            qlendict[line1[0].split()[0][1:].strip()] = len(line1[1].strip())
            lc = 0
            for line in f:
                qlendict[line.split()[0][1:].strip()] = len(next(f).strip())
        elif line1[0][0] == '@': #fastq
            qlendict[line1[0].split()[0][1:].strip()] = len(line1[1].strip())
            line2 = [next(f) for x in xrange(2)]
            lc = 0
            for line in f:
                if lc == 0:
                    qlendict[line.split()[0][1:].strip()] = len(next(f).strip())
                    lc = 1
                else:
                    lc += 1
                    if lc == 3:
                        lc = 0
    return qlendict

#Detects strandedness of an alignment
def strander(line):
    if line.split('\t')[8] < line.split('\t')[9]:
        return "+"
    elif line.split('\t')[8] > line.split('\t')[9]:
        return "-"

#Gathering alignment information and parsing
def info(line):
    qid = line.split('\t')[0]
    sid = line.split('\t')[1]
    piden = line.split('\t')[2]
    nmismatch = line.split('\t')[4]
    ngapopen = line.split('\t')[5]
    qstart = line.split('\t')[6]
    qstretch = int(line.split('\t')[7]) - int(line.split('\t')[6])
    sstart = min(int(line.split('\t')[8]), int(line.split('\t')[9]))
    sstretch = abs(int(line.split('\t')[9]) - int(line.split('\t')[8]))
    evalue = line.split('\t')[10]
    bitscore = line.split('\t')[11].strip('\n')
    strand = strander(line)
    qlen = qlendict[qid]
    entry = sid + '\t' + str(sstart) + '\t' + str(sstretch) + '\t+\t' + qid + '\t' + qstart + '\t' + str(qstretch) + '\t' + strand + '\t' + str(qlen) + '\t' + evalue + '\t' + bitscore + '\t' + piden + '\t' + nmismatch + '\t' + ngapopen
    return entry

def main():
    qlendict = measureQlen(file2)
    tmp = []
    with open(file1, 'r') as hsblast:
        for line in hsblast:
            tmp.append(info(line))
    tmp.append('null\tnull\tnull\tnull\tnull\tnull')
    l = len(tmp) - 1
    overlap_tolerance = 20
    temp = []
    temp2 = []
    output = []
    finalout = []
    #For Collecting chromosome number
    chromocollect = []
    for i in range(l):
        if tmp[i].split('\t')[4] == tmp[i+1].split('\t')[4]: #Grouping alignments by readname
            temp.append(tmp[i])
            if tmp[i].split('\t')[0].strip() not in chromocollect:
                chromocollect.append(tmp[i].split('\t')[0].strip())
        else:
            temp.append(tmp[i])
            if tmp[i].split('\t')[0].strip() not in chromocollect:
                chromocollect.append(tmp[i].split('\t')[0].strip())
            nchr = len(chromocollect)
            output = []
            j = len(temp)
            refrange = []
            refrange.append([int(temp[0].split('\t')[5]), int(temp[0].split('\t')[5]) + int(temp[0].split('\t')[6])])
            temp2.append(temp[0] + '\tn=' + str(j) + '\t' + str(nchr) + 'chr')
            for p in range(j)[1:]:
                queryx = [int(temp[p].split('\t')[5]), int(temp[p].split('\t')[5]) + int(temp[p].split('\t')[6])]
                queryrange = range(int(temp[p].split('\t')[5]), int(temp[p].split('\t')[5]) + int(temp[p].split('\t')[6]) + 1)
                clean = 1
                Trim = 0
                for r in refrange:
                    intersect = xrange(max(r[0], queryx[0]), min(r[1], queryx[1]) + 1)
                    if len(intersect) == 0: #Check if theres no alignment overlap
                        continue
                    else: #If there is alignment overlap
                        if len(intersect) <= overlap_tolerance: #Check if alignment overlap length is within limit
                            for t in intersect:#Trim alignment to remove overlapping
                                queryrange.remove(t)
                                Trim = 1
                        else: #If alignment overlap length exceed limit
                            clean = 0
                            break
                if clean == 1:
                    if Trim == 0:
                        refrange.append(queryx)
                        temp2.append(temp[p] + '\tn=' + str(j) + '\t' + str(nchr) + 'chr')
                    elif Trim == 1:
                        if len(queryrange) > 1: #Check if the trimmed query sequence is still acceptable in length
                            #Redefine alignment infomation
                            qstart = min(queryrange)
                            qend = max(queryrange)
                            qstretch = len(queryrange) - 1
                            refrange.append([qstart, qend])
                            temp2.append('\t'.join(temp[p].split('\t')[0:5]) + '\t' + str(qstart) + '\t' + str(qstretch) + '\t' + '\t'.join(temp[p].split('\t')[7:]) + '\tn=' + str(j) + '\t' + str(nchr) + 'chr')
            sortdict = {}
            for k in temp2:
                sortdict[int(k.split('\t')[5])] = k
            output = [value for (key, value) in sorted(sortdict.items())]
            for n in output:
                print n
            temp = []
            temp2 = []
            chromocollect = []

if __name__ == '__main__':
    main()
