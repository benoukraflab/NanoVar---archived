"""
This script creates the output VCF file.

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
from time import gmtime, strftime
import math

#Define Config of logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(message)s')
if len(argv)<10 or len(argv)>=11:
    sys.exit("Usage: python nv_parser_VCF.py overlap.DNN.cov.tsv Reference_genome.fa Long_read.fa Short_read_1.fa Short_read_2.fa bowcmd hsblastcmd genomesizes bedtools > output.tsv")

path_file = argv[1]
ref = argv[2]
longr = argv[3]
shortr1 = argv[4]
shortr2 = argv[5]
bedpath = argv[9]

bowcmd = open(argv[6], 'r').read().split('\n')[0]
hsblastcmd = open(argv[7], 'r').read().split('\n')[0]
gensize = open(argv[8], 'r').read().split('\n')

#Phred quality calculator
def phredC(p):
    probfalse = 1.000000001 - float(p)
    return math.log10(probfalse)*(-10)

def intersection(): #Tethering orphan Ins_bp SVs to parent without changing total score
    connectdict = {}
    weakread = {}
    for i in bed1:
        connectdict[i.split('\t')[3]] = []
    for i in intersect:
        if insbpdict.has_key(i.split('\t')[3]): #only Ins_bp SVs can be fellows
            if not weakread.has_key(i.split('\t')[8]):
                if float(i.split('\t')[9]) >= float(i.split('\t')[4]):
                    connectdict[i.split('\t')[8]].append(i.split('\t')[3])
                    weakread[i.split('\t')[3]] = 1
    return connectdict,weakread

#Read file
rdata = open(path_file, 'r').read()

#Open temporary file 
tmpread = []

#Coverages
covl = 0
covs = 0

#Run stopper
stop = 0

#Overall query counter
n = 0

#Temporary query counter
c = 0

#Calculating number of queries
l = len(rdata.split('\n'))

#Set range according to number of queries
k = range(l-1)

#Convert string to list object
rdata = list(rdata.split('\n'))

#Create last line dummy
rdata[l-1] = 'dum' + '\t' + 'dum' + '\t' + 'dum' + '\t' + 'dum' + '\t' + 'dum' + '\t' + 'dum' + '\t' + 'dum' + '\t' + 'dum' + '\t' + 'dum'

print '##fileformat=VCFv4.2'
print '##fileDate=' + strftime("%Y-%m-%d", gmtime())
print '##source_longread=' + longr 
print '##source_shortreadpair1=' + shortr1
print '##source_shortreadpair2=' + shortr2
print '##reference=' + ref
print '##longread-to-ref-alignment=HS-BLAST=' + hsblastcmd
print '##shortread-to-longread-alignment=Bowtie2=' + bowcmd
print '##phasing=none'
#print '##ALT=<ID=NOV-INS,Description="Novel-sequence insertion">'
print '##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">'
print '##ALT=<ID=DEL,Description="Sequence deletion relative to the reference">'
print '##ALT=<ID=INV,Description="Sequence inversion relative to the reference">'
print '##ALT=<ID=DUP,Description="Sequence tandem duplication relative to the reference">'
print '##ALT=<ID=BND,Description="Sequence breakend relative to the reference (Possibly translocations/genomic-sequence insertions)">'
print '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
print '##INFO=<ID=CHR2,Number=1,Type=String,Description="End chromosome number of the described variant">'
print '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the described variant">'
print '##INFO=<ID=SCOV,Number=1,Type=Float,Description="Average number of short-reads covering the SV breakpoint">'
print '##INFO=<ID=LCOV,Number=1,Type=Float,Description="Number of long-reads supporting the SV breakpoint">'
print '##INFO=<ID=NCOV,Number=1,Type=Integer,Description="Number of long-reads not supporting the SV breakpoint">'
print '##INFO=<ID=SVRATIO,Number=1,Type=Float,Description="Ratio of SV reads to total reads at SV breakpoint">'
print '##INFO=<ID=PROB,Number=1,Type=Float,Description="Neural network probability score of variant being true">'
print '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Estimated length of the SV">'
print '##INFO=<ID=SVLEN_UNKN,Number=0,Type=Flag,Description="Unknown length of SV">'
#print '##INFO=<ID=NOVEL,Number=0,Type=Flag,Description="Indicates a novel sequence insertional structural variant">'

lg = len(gensize) -1
for i in range(lg):
    print '##contig=<ID=' + gensize[i].split('\t')[0] + ',length=' + gensize[i].split('\t')[1] + '>'

print '#CHROM' + '\t' + 'POS' + '\t' + 'ID' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'QUAL' + '\t' + 'FILTER' + '\t' + 'INFO'

def negtozero(x):
    if x < 0:
        return 0
    else:
        return x

def num(x):
    try:
        return float(x)
    except:
        return float(0)

out = []

insbpdict = {}

for i in k:
    n = n + 1
    #Sensing consecutive SV breakpoint unique identifier
    if rdata[i].split('\t')[8] == rdata[i+1].split('\t')[8]:
        tmpread.append(rdata[i])
        c = c + 1
    else:
        tmpread.append(rdata[i])
        c = c + 1
        covs = 0
        covl = float(tmpread[0].split('\t')[10])
        normcov = float(tmpread[0].split('\t')[12])
        DNNscore = float(tmpread[0].split('\t')[13])
        for j in range(c):
            covs += num(tmpread[j].split('\t')[14])
        covl = str(round(float(covl), 2))
        covs = str(round(float(covs)/c, 2))
        DNN = str(round(float(DNNscore), 4))
        phre = phredC(DNNscore)
        phred = str(round(float(phre), 1))
        sv_id = tmpread[0].split('\t')[6].split('~')[0] + '~' + tmpread[0].split('\t')[0]
        bp_name = tmpread[0].split('\t')[3].split(' ')[0]
        chrm1 = tmpread[0].split('\t')[6].split('~')[1].split(':')[0]
        if bp_name == 'Nov_Ins':
            sv = 'INS'
            sv_len = tmpread[0].split('\t')[3].split(' ')[1]
            chrm2 = '.' #tmpread[0].split('\t')[6].split('~')[1].split(':')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            coord2 = '.' #int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(sv_id) + '\t' + '.' + '\t' + str(sv) + '\t' + str(phred) + '\t' + 'PASS' + '\t' + 'SVTYPE=' + str(sv) + ';' + 'CHR2=' + str(chrm2) + ';' + 'END=' + str(coord2) + ';' + 'SCOV=' + str(covs) + ';' + 'LCOV=' + str(covl) + ';' + 'NCOV=' + str(normcov) + ';' + 'PROB=' + str(DNN) + ';' + 'SVLEN=' + str(sv_len) + ';')
        elif bp_name == 'E-Nov_Ins_bp' or bp_name == 'S-Nov_Ins_bp':
            insbpdict[sv_id] = 1
            sv = 'BND' 
            sv_len = tmpread[0].split('\t')[3].split(' ')[1]
            chrm2 = '.'
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            coord2 = '.'
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(sv_id) + '\t' + '.' + '\t' + str(sv) + '\t' + str(phred) + '\t' + 'PASS' + '\t' + 'SVTYPE=' + str(sv) + ';' + 'CHR2=' + str(chrm2) + ';' + 'END=' + str(coord2) + ';' + 'SCOV=' + str(covs) + ';' + 'LCOV=' + str(covl) + ';' + 'NCOV=' + str(normcov) + ';' + 'PROB=' + str(DNN) + ';' + 'SVLEN=' + str(sv_len) + ';')
        elif bp_name == 'Del':
            sv = 'DEL'
            sv_len = tmpread[0].split('\t')[3].split(' ')[1]
            chrm2 = tmpread[0].split('\t')[6].split('~')[1].split(':')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(sv_id) + '\t' + '.' + '\t' + str(sv) + '\t' + str(phred) + '\t' + 'PASS' + '\t' + 'SVTYPE=' + str(sv) + ';' + 'CHR2=' + str(chrm2) + ';' + 'END=' + str(coord2) + ';' + 'SCOV=' + str(covs) + ';' + 'LCOV=' + str(covl) + ';' + 'NCOV=' + str(normcov) + ';' + 'PROB=' + str(DNN) + ';' + 'SVLEN=' + str(sv_len) + ';')
        elif bp_name == 'Inv':
            sv = 'INV'
            sv_len = '.'
            chrm2 = tmpread[0].split('\t')[6].split('~')[1].split(':')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(sv_id) + '\t' + '.' + '\t' + str(sv) + '\t' + str(phred) + '\t' + 'PASS' + '\t' + 'SVTYPE=' + str(sv) + ';' + 'CHR2=' + str(chrm2) + ';' + 'END=' + str(coord2) + ';' + 'SCOV=' + str(covs) + ';' + 'LCOV=' + str(covl) + ';' + 'NCOV=' + str(normcov) + ';' + 'PROB=' + str(DNN) + ';' + 'SVLEN_UNKN' + ';')
        elif bp_name == 'Inv(1)' or bp_name == 'Inv(2)':
            sv = 'INV'
            sv_len = tmpread[0].split('\t')[3].split(' ')[1]
            chrm2 = tmpread[0].split('\t')[6].split('~')[1].split(':')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(sv_id) + '\t' + '.' + '\t' + str(sv) + '\t' + str(phred) + '\t' + 'PASS' + '\t' + 'SVTYPE=' + str(sv) + ';' + 'CHR2=' + str(chrm2) + ';' + 'END=' + str(coord2) + ';' + 'SCOV=' + str(covs) + ';' + 'LCOV=' + str(covl) + ';' + 'NCOV=' + str(normcov) + ';' + 'PROB=' + str(DNN) + ';' + 'SVLEN=' + str(sv_len) + ';')
        elif bp_name == 'TDupl':
            sv = 'DUP'
            sv_len = '.'
            chrm2 = tmpread[0].split('\t')[6].split('~')[1].split(':')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(sv_id) + '\t' + '.' + '\t' + str(sv) + '\t' + str(phred) + '\t' + 'PASS' + '\t' + 'SVTYPE=' + str(sv) + ';' + 'CHR2=' + str(chrm2) + ';' + 'END=' + str(coord2) + ';' + 'SCOV=' + str(covs) + ';' + 'LCOV=' + str(covl) + ';' + 'NCOV=' + str(normcov) + ';' + 'PROB=' + str(DNN) + ';' + 'SVLEN_UNKN' + ';')
        elif bp_name == 'Intra-Ins':
            sv = 'BND'
            sv_len = '.'
            chrm2 = tmpread[0].split('\t')[6].split('~')[1].split(':')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(sv_id) + '\t' + '.' + '\t' + str(sv) + '\t' + str(phred) + '\t' + 'PASS' + '\t' + 'SVTYPE=' + str(sv) + ';' + 'CHR2=' + str(chrm2) + ';' + 'END=' + str(coord2) + ';' + 'SCOV=' + str(covs) + ';' + 'LCOV=' + str(covl) + ';' + 'NCOV=' + str(normcov) + ';' + 'PROB=' + str(DNN) + ';' + 'SVLEN_UNKN' + ';')
        elif bp_name == 'Intra-Ins(1)' or bp_name == 'Intra-Ins(2)':
            sv = 'BND'
            sv_len = tmpread[0].split('\t')[3].split(' ')[1]
            chrm2 = tmpread[0].split('\t')[6].split('~')[1].split(':')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(sv_id) + '\t' + '.' + '\t' + str(sv) + '\t' + str(phred) + '\t' + 'PASS' + '\t' + 'SVTYPE=' + str('Ins') + ';' + 'CHR2=' + str(chrm2) + ';' + 'END=' + str(coord2) + ';' + 'SCOV=' + str(covs) + ';' + 'LCOV=' + str(covl) + ';' + 'NCOV=' + str(normcov) + ';' + 'PROB=' + str(DNN) + ';' + 'SVLEN=' + str(sv_len) + ';')
        elif bp_name == 'Inter-Ins(1)' or bp_name == 'Inter-Ins(2)':
            sv = 'BND'
            sv_len = tmpread[0].split('\t')[3].split(' ')[1]
            chrm2 = tmpread[0].split('\t')[6].split('~')[2].split(':')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[2].split(':')[1])
            #chrmdict = {chrm1:coord1, chrm2:coord2}
            #chrmlist = sorted([chrm1, chrm2])
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(sv_id) + '\t' + '.' + '\t' + str(sv) + '\t' + str(phred) + '\t' + 'PASS' + '\t' + 'SVTYPE=' + str('Ins') + ';' + 'CHR2=' + str(chrm2) + ';' + 'END=' + str(coord2) + ';' + 'SCOV=' + str(covs) + ';' + 'LCOV=' + str(covl) + ';' + 'NCOV=' + str(normcov) + ';' + 'PROB=' + str(DNN) + ';' + 'SVLEN=' + str(sv_len) + ';')
        elif bp_name == 'InterTx':
            sv = 'BND'
            sv_len = '.'
            chrm2 = tmpread[0].split('\t')[6].split('~')[2].split(':')[0]
            coord1 = int(tmpread[0].split('\t')[6].split('~')[1].split(':')[1])
            coord2 = int(tmpread[0].split('\t')[6].split('~')[2].split(':')[1])
            #chrmdict = {chrm1:coord1, chrm2:coord2}
            #chrmlist = sorted([chrm1, chrm2])
            out.append(str(chrm1) + '\t' + str(coord1) + '\t' + str(sv_id) + '\t' + '.' + '\t' + str(sv) + '\t' + str(phred) + '\t' + 'PASS' + '\t' + 'SVTYPE=' + str(sv) + ';' + 'CHR2=' + str(chrm2) + ';' + 'END=' + str(coord2) + ';' + 'SCOV=' + str(covs) + ';' + 'LCOV=' + str(covl) + ';' + 'NCOV=' + str(normcov) + ';' + 'PROB=' + str(DNN) + ';' + 'SVLEN_UNKN' + ';')
        else:
            print 'ERROR in breakpoint name'
            sys.exit()
            break
        tmpread = []
        c=0

################################Redundant step, overlap.py was updated, edited to removed step (23/9/18)#######
bed1 = []
bed2 = []
outdict = {}
for i in out:
    bed1.append('\t'.join(i.split('\t')[0:2]) + '\t' + str(int(i.split('\t')[1]) + 1) + '\t' + i.split('\t')[2] + '\t' + i.split('\t')[5])
    if int(i.split('\t')[1]) - 200 < 0:
        coord1 = 0
    else:
        coord1 = int(i.split('\t')[1]) - 200
    bed2.append(i.split('\t')[0] + '\t' + str(coord1) + '\t' + str(int(i.split('\t')[1]) + 200) + '\t' + i.split('\t')[2] + '\t' + i.split('\t')[5])
    outdict[i.split('\t')[2]] = i

bedout1 = open('hsblast_longreads/bed1', 'w')
bedout1.write('\n'.join(bed1))
bedout1.close()
bedout2 = open('hsblast_longreads/bed2', 'w')
#bedout2.write('\n'.join(bed2))
bedout2.close()
os.system(bedpath + ' sort -i hsblast_longreads/bed1 > hsblast_longreads/bed1.sort')
os.system(bedpath + ' sort -i hsblast_longreads/bed2 > hsblast_longreads/bed2.sort')
os.system(bedpath + " intersect -wa -wb -a hsblast_longreads/bed1.sort -b hsblast_longreads/bed2.sort | awk -F'\t' '{if ($4 != $9) print $0}' > hsblast_longreads/intersect.txt")
intersect = open('hsblast_longreads/intersect.txt', 'r').read().splitlines()
connectdict,weakread = intersection()
for key in connectdict:
    if not weakread.has_key(key):
        lcov = 0
        for i in connectdict[key]:
            lcov += float(outdict[i].split('\t')[7].split(';')[4].split('=')[1])
        lcov += float(outdict[key].split('\t')[7].split(';')[4].split('=')[1])
        print '\t'.join(outdict[key].split('\t')[0:7]) + '\t' + ';'.join(outdict[key].split('\t')[7].split(';')[0:4]) + ';LCOV=' + str(lcov) + ';' + outdict[key].split('\t')[7].split(';')[5] + ';SVRATIO=' + str(round(float(lcov)/(float(lcov) + float(outdict[key].split('\t')[7].split(';')[5].split('=')[1])),3)) + ';' + ';'.join(outdict[key].split('\t')[7].split(';')[6:])
