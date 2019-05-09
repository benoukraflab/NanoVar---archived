"""
This script crudely corrects for Nanopore sequencing errors within each read.

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

import os
from sys import argv
import sys
import re
import logging

if len(argv)<6 or len(argv)>=7:
    sys.exit("Usage: python nv_nanocorrector.py [SV+.rlen.tsv] [longreads.fa] [reference.fa] [overlap.DNN.tsv] [path_to_log_file] > output.fa")

#Define paths
file1 = argv[1]
longfa = argv[2]
reffa = argv[3]
overlapDNN = argv[4]
log_path=argv[5]

logging.basicConfig(filename=log_path, level=logging.DEBUG, format='%(asctime)s - %(message)s')

#Create dictionary for query name to query bed
def querymapper(x):
    bedoutq = {}
    g = 0
    f = 0
    readlength = int(seqinfo[x].split('\t')[1])
    qname = str(seqinfo[x].split('\t')[0])
    #Start coord
    if int(eval(seqinfo[x].split('\t')[4])[0][0]) > 1: # actually it is -1
        start1 = 0
        start2 = int(eval(seqinfo[x].split('\t')[4])[0][0]) - 1
        bedoutq[qname + '*1'] = [str(start1), str(start2), '1']
        g = 1
        f = 1
    #Middle coord
    j = len(eval(seqinfo[x].split('\t')[4])) - 1
    if j != 0:
        for i in range(j):
            if int(eval(seqinfo[x].split('\t')[4])[i+1][0]) - int(eval(seqinfo[x].split('\t')[4])[i][1]) > -1:
                bedoutq[qname + '*' + str(g + 2)] = [str(eval(seqinfo[x].split('\t')[4])[i][1]), str(eval(seqinfo[x].split('\t')[4])[i+1][0] - 1), str(g + 2)]
                g = g + 2
    #End coord
    if readlength - int(eval(seqinfo[x].split('\t')[4])[-1][1]) > 0: # actually it is -1
        end1 = int(eval(seqinfo[x].split('\t')[4])[-1][1])
        end2 = readlength
        bedoutq[qname + '*' + str(g + 2)] = [str(end1), str(end2), str(g + 2)]
        g = g + 2
    return f,g,bedoutq

#Create a dictionary of reference chromosome names to its sequence
def fastadict(x):
    ref = open(x, 'r').read()
    refn = ref.split('\n')
    refh = ref.replace('\n','').split('>')
    ref = ''
    refnlen = []
    l = len(refn) -1
    for i in range(l):
        if refn[i][0] == '>':
            refnlen.append(len(refn[i]))
    refn = ''
    l = len(refnlen)
    refdict = {}
    for i in range(l):
        refdict['>' + refh[i+1][0:refnlen[i]-1]] = refh[i+1][refnlen[i]-1:]
    return refdict

#Reverse Complement
def rc(seq):
    s = seq[::-1]
    seq = []
    for i in s:
        if i == 'A':
            seq.append('T')
        elif i == 'T':
            seq.append('A')
        elif i == 'C':
            seq.append('G')
        elif i == 'G':
            seq.append('C')
        elif i == 'a':
            seq.append('t')
        elif i == 't':
            seq.append('a')
        elif i == 'c':
            seq.append('g')
        elif i == 'g':
            seq.append('c')    
        elif i == 'N':
            seq.append('N')
    seq = ''.join(seq)
    return seq

#Create dictionary for subject name to subject bed
def subjmapper(x,y):
    bedouts = {}
    g = y -1
    j = len(eval(seqinfo[x].split('\t')[2]))
    for i in range(j):
        if int(eval(seqinfo[x].split('\t')[3])[i][0]) < int(eval(seqinfo[x].split('\t')[3])[i][1]): # Pos strand
            chrom = str(eval(seqinfo[x].split('\t')[2])[i])
            begin = int(eval(seqinfo[x].split('\t')[3])[i][0]) - 1
            end = int(eval(seqinfo[x].split('\t')[3])[i][1])
            strand = '+'
            bedouts[chrom + '*' + str(g + 2)] = [str(begin), str(end), str(g + 2), str(1), str(strand), str(chrom)]
            g = g + 2
        elif int(eval(seqinfo[x].split('\t')[3])[i][0]) > int(eval(seqinfo[x].split('\t')[3])[i][1]): # Neg strand
            chrom = str(eval(seqinfo[x].split('\t')[2])[i])
            begin = int(eval(seqinfo[x].split('\t')[3])[i][1]) -1
            end = int(eval(seqinfo[x].split('\t')[3])[i][0])
            strand = '-'
            bedouts[chrom + '*' + str(g + 2)] = [str(begin), str(end), str(g + 2), str(1), str(strand), str(chrom)]
            g = g + 2
    return bedouts

def fastsearch(x,t,e):
    qfasta = ''
    sfasta = ''
    qname = str(seqinfo[x].split('\t')[0])
    longfasta = longread[e+1]
    for i in bedoutq:
        qfasta += '>' + str(bedoutq[i][2]) + '\n' + longfasta[int(bedoutq[i][0]):int(bedoutq[i][1])] + '\n'
    for i in bedouts:
        if bedouts[i][4] == '+':
            sfasta += '>' + str(bedouts[i][2]) + '\n' + refdict['>' + bedouts[i][5]][int(bedouts[i][0]):int(bedouts[i][1])] + '\n'
        elif bedouts[i][4] == '-':
            sfasta += '>' + str(bedouts[i][2]) + '\n' + rc(str(refdict['>' + bedouts[i][5]][int(bedouts[i][0]):int(bedouts[i][1])])) + '\n'
    qfasta = list(qfasta.split('\n'))
    sfasta = list(sfasta.split('\n'))
    w = len(qfasta) - 1
    v = len(sfasta) - 1
    fast = '>' + qname + '\n'
    for i in range(t):
        c = i + 1
        for i in range(w):
            if qfasta[i][0:2] == '>' + str(c):
                fast += qfasta[i+1]
                break
        for i in range(v):
            if sfasta[i][0:2] == '>' + str(c):
                fast += sfasta[i+1]
                break 
    return fast

logging.info('Start read correcting...')

sequenceinfo = os.popen("cut -f 1,2,23-27 " + file1).read()
seqinfo = list(sequenceinfo.split('\n'))

svdict = {}
l = len(seqinfo)
#Create dictionary of selected SV reads to their index
for i in range(l-1):
    svdict[seqinfo[i].split('\t')[0]] = i

DNNsvdict = {}
with open(overlapDNN) as f:
    for line in f:
        DNNsvdict[line.split('\t')[0]] = svdict[line.split('\t')[0]]

longread = open(longfa, 'r').read().split('\n')
lr = range(len(longread) - 1)

if longread[0][0] == '>':
    lre = [k for k in lr if k%2 == 0]
elif longread[0][0] == '@':
    lre = [k for k in lr if k%4 == 0]
refdict = fastadict(reffa)

for i in lre:
    try:
        if DNNsvdict[longread[i].split()[0][1:]] >= 0:
            q = DNNsvdict[longread[i][1:]]
            d,s,bedoutq = querymapper(q)
            bedouts = subjmapper(q,d)
            print fastsearch(q,s,i)
    except:
        pass

logging.info("Read correction Finished")
