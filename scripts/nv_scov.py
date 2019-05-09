"""
This script clusters short reads around each SV breakend.

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

#Ensure correct number of inputs
if len(argv)<5 or len(argv)>=6:
    sys.exit("Usage: python nv_scov.py ANN.tsv short-long.bam buffer bedtoolpath > ANN.cov.tsv")

#Assign variable to inputs
file1 = argv[1]
bam = argv[2]
buf = int(argv[3])
bedpath = argv[4]

#Read file
rdata = open(file1, 'r').read().splitlines()

#Open temporary file 
tmpread = []

#Calculating number of queries
l = len(rdata)

#Set range according to number of queries
k = range(l)

#Create last line dummy
rdata.append('dum\tdum\tdum\tdum\tdum\tdum\tdum\tdum\tdum')

out = []

for i in k:
    #Sensing consecutive SV breakpoint unique identifier
    if rdata[i].split('\t')[8] == rdata[i+1].split('\t')[8]:
        tmpread.append(rdata[i])
    else:
        tmpread.append(rdata[i])
        if len(tmpread) == 2:
            if int(tmpread[0].split('\t')[1]) >= buf:
                out.append(tmpread[0].split('\t')[0] + '\t' + str(int(tmpread[0].split('\t')[1]) - buf) + '\t' + str(int(tmpread[0].split('\t')[1]) + buf) + '\t' + tmpread[0].split('\t')[8])
            else:
                out.append(tmpread[0].split('\t')[0] + '\t0\t' + str(int(tmpread[0].split('\t')[1]) + buf) + '\t' + tmpread[0].split('\t')[8])
            if int(tmpread[1].split('\t')[2]) >= buf:
                out.append(tmpread[1].split('\t')[0] + '\t' + str(int(tmpread[1].split('\t')[2]) - buf) + '\t' + str(int(tmpread[1].split('\t')[2]) + buf) + '\t' + tmpread[1].split('\t')[8])
            else:
                out.append(tmpread[1].split('\t')[0] + '\t0\t' + str(int(tmpread[1].split('\t')[2]) + buf) + '\t' + tmpread[1].split('\t')[8]) 
        elif len(tmpread) == 1:
            if int(tmpread[0].split('\t')[1]) >= buf:
                out.append(tmpread[0].split('\t')[0] + '\t' + str(int(tmpread[0].split('\t')[1]) - buf) + '\t' + str(int(tmpread[0].split('\t')[1]) + buf) + '\t' + tmpread[0].split('\t')[8])
            else:
                out.append(tmpread[0].split('\t')[0] + '\t0\t' + str(int(tmpread[0].split('\t')[1]) + buf) + '\t' + tmpread[0].split('\t')[8])
        else:
            print 'SCOV.py error'
            break
            sys.exit()
        tmpread = []

sbed = open('sbed', 'w')
sbed.write('\n'.join(out))
sbed.close()
os.system(bedpath + ' sort -i sbed > sbed.sort')
os.system(bedpath + " coverage -a sbed.sort -b " + bam + " > scov.bed")
scoverage = open('scov.bed', 'r').read().splitlines()
scovdict = {}
for i in scoverage:
    if scovdict.has_key(i.split('\t')[3]):
        scovdict[i.split('\t')[3]] += float(i.split('\t')[4])
    else:
        scovdict[i.split('\t')[3]] = float(i.split('\t')[4])

for i in k:
    #Sensing consecutive SV breakpoint unique identifier
    if rdata[i].split('\t')[8] == rdata[i+1].split('\t')[8]:
        tmpread.append(rdata[i])
    else:
        tmpread.append(rdata[i])
        if len(tmpread) == 2:
            print tmpread[0] + '\t' + str(float(scovdict[tmpread[0].split('\t')[8]])/2)
            print tmpread[1] + '\t' + str(float(scovdict[tmpread[0].split('\t')[8]])/2)
        elif len(tmpread) == 1:
            print tmpread[0] + '\t' + str(float(scovdict[tmpread[0].split('\t')[8]]))
        tmpread = []
