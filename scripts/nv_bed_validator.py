"""
This script verifies input BED files.

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

if len(argv)<3 or len(argv)>=4:
    sys.exit("Usage: python nv_bed_validator.py bed.file ref_prefix.genomesizes")

bed_file = argv[1]
ref_file = argv[2]

ref = open(ref_file, 'r').read().splitlines()

chromsizedict = {}
chromdict = {}
p = 1
c = 0

for i in ref:
    chromsizedict[i.split('\t')[0]] = [0, i.split('\t')[1]]
    chromdict[i.split('\t')[0]] = []

with open(bed_file) as bed:
    for line in bed:
        c += 1
        try:
            if chromdict[line.split('\t')[0]] == []:
                try:
                    if int(chromsizedict[line.split('\t')[0]][0]) > int(line.split('\t')[1]) or int(line.split('\t')[2]) > int(chromsizedict[line.split('\t')[0]][1]):
                        p = 0
                        print str(c) + " invalidbedregion"
                        break
                    else:
                        if int(line.split('\t')[1]) > int(line.split('\t')[2]):
                            p = 0
                            print str(c) + " invalidbedrange"
                            break
                except:
                    p = 0
                    print str(c) + " invalidbedvalue"
                    break
        except:
            p = 0
            print str(c) + " invalidbedchrom"
            break

if p == 1:
    print "bedpass"
