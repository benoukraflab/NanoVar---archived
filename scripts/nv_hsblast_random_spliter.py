"""
This script splits aligned reads randomly into groups.

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
import random

#Ensure correct number of inputs
if len(argv)<4 or len(argv)>=5:
    sys.exit("Usage: python nv_hsblast_random_spliter.py file.parse.tsv split-number seed")

#Assign variable to inputs
file1 = argv[1]
splitno = int(argv[2])
seed = int(argv[3])

#Splitting parse.tsv data according to fasta name random sampling
def splitter(parse):
    random.seed(seed)
    pdata = open(parse, 'r').read().splitlines()
    n = len(pdata)
    pdata.append('dum\tdum\tdum\tdum\tdum\tdum\tdum\tdum\tdum\tdum\tdum')
    temp = []
    outputdict = {}
    for i in range(splitno):
        outputdict[i + 1] = open('parse_' + str(i+1) + '.tsv', 'w')
    for i in xrange(n):
        if pdata[i].split('\t')[0] == pdata[i+1].split('\t')[0]:
            temp.append(pdata[i])
        else:
            temp.append(pdata[i])
            nrand = random.randint(1,splitno)
            outputdict[nrand].write('\n'.join(temp) + '\n')
            temp = []
    for i in range(splitno):
        outputdict[i+1].close()

def main():
    splitter(file1)

if __name__ == '__main__':
    main()
