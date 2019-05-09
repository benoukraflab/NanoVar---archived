"""
This script verifies long-read FASTA files.

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

if len(argv)<2 or len(argv)>=3:
    sys.exit("Usage: python nv_fasta_validator.py longread.fasta")

longr = argv[1]

c = 0
lc = 0
tc = 0
corrupt = 0
qtype = ''
totalbases = 0

with open(longr) as f:
    try:
        line1 = [next(f) for x in xrange(2)]
        line2 = [next(f) for x in xrange(2)]
    except:
        pass
    try:
        if line2:
            if line1[0][0] == '>':
                if line2[0][0] == '>': #fasta format
                    qtype = 'fasta'
                    for line in f:
                        c += 1
                        if lc == 0:
                            if line[0] != '>':
                                corrupt = 1
                                break
                            lc = 1
                        else:
                            #totalbases += len(line)
                            lc = 0
                            if c == 100000:
                                tc = tc + c
                                c = 0
                else:
                    c = 3
                    corrupt = 1
            elif line1[0][0] == '@':
                if line2[0][0] == '+': #fastq format
                    qtype = 'fastq'
                    for line in f:
                        c += 1
                        if lc == 0:
                            if line[0] != '@':
                                corrupt = 1
                                break
                            lc = 1
                        else:
                            #if lc == 1:
                                #totalbases += len(line)
                            lc += 1
                            if lc == 4:
                                lc = 0
                            if c == 100000:
                                tc = tc + c
                                c = 0
                else:
                    c = 3
                    corrupt = 1
            else:
                c = 1
                corrupt = 1
    except:
        try:
            if line1:
                if line1[0][0] == '>':
                    c = 2
                elif line1[0][0] == '@':
                    c = 2
                    corrupt = 1
                else:
                    c = 1
                    corrupt = 1
        except:
            c = 0
            corrupt = 2

if corrupt == 0:
    if qtype == 'fasta':
        if c%2 == 0:
            print str(((tc+c)/2)+2)# + ' ' + str(totalbases)
        else:
            print str(tc+c) + ' corrupt'
    elif qtype == 'fastq':
        if c%4 == 0:
            print str(((tc+c)/4)+1)# + ' ' + str(totalbases)
        else:
            print str(tc+c) + ' corrupt'
elif corrupt == 1:
    print str(tc+c) + ' corrupt'
else:
    print 'empty'
