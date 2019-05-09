"""
This script verifies short-read FASTA files.

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
    sys.exit("Usage: python nv_fastq-a_validator.py shortread.fasta")

shortr = argv[1]
c = 0

with open(shortr) as s:
    try:
        line1 = [next(s) for x in xrange(2)]
        line2 = [next(s) for x in xrange(2)]
    except:
        pass
    try:
        if line2:
            if line1[0][0] == '>':
                if line2[0][0] == '>':
                    print 'fasta'
                else:
                    print 'corrupt'
            elif line1[0][0] == '@':
                if line2[0][0] == '+':
                    print 'fastq'
                else:
                    print 'corrupt'
            else:
                print 'corrupt'
    except:
        if line1[0][0] == '>':
            print 'fasta'
        elif line1[0][0] == '@': 
            print 'fastq'
        else:
            print 'corrupt'