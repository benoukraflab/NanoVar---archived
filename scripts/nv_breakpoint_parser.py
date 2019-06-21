"""
This script parses the output from nv_hsblast_SV_detector.py.

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
import random
import ast

#Seed
random.seed(1)

#Ensure correct number of inputs
if len(argv)<4 or len(argv)>=5:
    sys.exit("Usage: python nv_breakpoint_parser.py file.SV.tsv SV_length_cutoff[INT] path_to_log_file > output.tsv")

#Define path to SV file
file1 = argv[1]
cutoff = int(argv[2]) - 1
log_path=argv[3]

#Define Config of logging
logging.basicConfig(filename=log_path, level=logging.DEBUG, format='%(asctime)s - %(message)s')
breakpointfile = os.popen("cut -f 1,3-18,20-22,25-28 " + file1).read()
#Ouput columns 0=readname, 1=total_chromosome(s), 2=No._of_maps, 3=query_signature, 4=evalue, 5=bitscore, 6=deletion, 7=Isolated_insertion, 8=Insertion, 9=complex_insertion, 10=tandemDupl, 11=inversion, 12=intra_insertion, 13=inter_insertion, 14=inter_tx, 15=majstrand, 16=short_frag_perc, 17=complex_SV, 18=Total_SV_complex, 19=sv_range, 20=query_string, 21=%identity, 22=mismatch, 23=gap_ratio

breakpointfile = list(breakpointfile.split('\n'))
ran = "01234567890ABCDEFGHIJKLMNOPQRSTUVWXYZ"
ran_len = 5
l = len(breakpointfile) -1
c = range(l)
g = 0

#Function to convert chromosome names to numbers for ranking (Not used)
def chrnum(chrm):
    return chrm

'''
    if chrm[3:] == 'X':
        rank = 23
    elif chrm[3:] == 'Y':
        rank = 24
    elif chrm[3:] == 'M':
        rank = 25
    else:
        try:
            rank = int(chrm[3:])
        except:
            rank = str(chrm)
    return rank
'''

#Function to convert numbers back to respective chromosome names (Not used)
def revchrnum(rank):
    return rank

'''
    if rank == 23:
        chrm = 'chrX'
    elif rank == 24:
        chrm = 'chrY'
    elif rank == 25:
        chrm = 'chrM'
    else:
        try:
            chrm = 'chr' + str(int(rank))
        except:
            chrm = rank
    return chrm
'''

index = [0, 2, 4]
def getSignature(bp1, name, querymap, sign, realnmaps, elist, bitlist, bps_no, complex_sv, nchrom, short, piden, mismatch, gap_ratio):
    signlist = sign.split(',')
    if name == 'S-Nov_Ins':
        signature = ','.join(signlist[0:3]) + ',0,0,' + str(elist[0]) + ',0,' + str(bitlist[0]) + ',0,' + str(complex_sv) + ',' + str(realnmaps) + ',' + str(bps_no) + ',' + str(nchrom) + ',' + str(short)  + ',' + str(piden[0]) + ',0,' + str(mismatch[0]) + ',0,' + str(gap_ratio[0]) + ',0'
    elif name == 'E-Nov_Ins':
        signature = ','.join(signlist[-3:]) + ',0,0,' + str(elist[-1]) + ',0,' + str(bitlist[-1]) + ',0,' + str(complex_sv) + ',' + str(realnmaps) + ',' + str(bps_no) + ',' + str(nchrom) + ',' + str(short)  + ',' + str(piden[-1]) + ',0,' + str(mismatch[-1]) + ',0,' + str(gap_ratio[-1]) + ',0'
    else:
        maps = ast.literal_eval(querymap)
        nmaps = len(maps)
        for i in range(nmaps):
            if int(maps[i][1]) == bp1:
                signature = ','.join(signlist[index[i]:index[i]+5]) + ',' + str(elist[i]) + ',' + str(elist[i+1]) + ',' + str(bitlist[i]) + ',' + str(bitlist[i+1]) + ',' + str(complex_sv) + ',' + str(realnmaps) + ',' + str(bps_no) + ',' + str(nchrom) + ',' + str(short)  + ',' + str(piden[i])  + ',' + str(piden[i+1]) + ',' + str(mismatch[i]) + ',' + str(mismatch[i+1]) + ',' + str(gap_ratio[i]) + ',' + str(gap_ratio[i+1])
                break
    return signature

def trysort(a,b):
    try:
        return sorted([a,b])
    except:
        return [a,b]

logging.info('Start SV Parsing...')
for i in c:
    g = g + 1
    read_name = breakpointfile[i].split('\t')[0]
    if breakpointfile[i].split('\t')[17] == '':
        complex_sv = 0
    else:
        complex_sv = len(breakpointfile[i].split('\t')[17].split(','))
    if int(breakpointfile[i].split('\t')[17].count(',')) >= 1:
        multi_sv = len(breakpointfile[i].split('\t')[18].split(','))
    else:
        multi_sv = 0
    del_count = del_counter = int(breakpointfile[i].split('\t')[19].count('Del'))
    iso_ins_count = iso_ins_counter = int(breakpointfile[i].split('\t')[19].count('S-Nov_Ins')) + int(breakpointfile[i].split('\t')[19].count('E-Nov_Ins'))
    ins_count = ins_counter = int(breakpointfile[i].split('\t')[19].count('Nov_Ins'))
    tdup_count = tdup_counter = int(breakpointfile[i].split('\t')[19].count('TDupl'))
    inv_count = inv_counter = int(breakpointfile[i].split('\t')[19].count('Inv'))
    intra_ins_count = intra_ins_counter = int(breakpointfile[i].split('\t')[19].count('Intra-Ins'))
    inter_ins_count = inter_ins_counter = int(breakpointfile[i].split('\t')[19].count('Inter-Ins'))
    inter_tx_count = inter_tx_counter = int(breakpointfile[i].split('\t')[19].count('InterTx'))
    bps = breakpointfile[i].split('\t')[19].split(',')
    bps_no = len(bps)
    querymap = breakpointfile[i].split('\t')[20]
    nchrom = breakpointfile[i].split('\t')[1]
    evaluelist = breakpointfile[i].split('\t')[4].split(',')
    bitscorelist = breakpointfile[i].split('\t')[5].split(',')
    short = breakpointfile[i].split('\t')[16]
    piden = breakpointfile[i].split('\t')[21].split(',')
    mismatch = breakpointfile[i].split('\t')[22].split(',')
    gap_ratio = breakpointfile[i].split('\t')[23].split(',')
    sign = breakpointfile[i].split('\t')[3]
    realnmaps = int(breakpointfile[i].split('\t')[2].split(' ')[0])
    for j in range(bps_no):
        bp_name = bps[j].split(':')[0]
        bp_range = re.sub('^.*:', '', str(bps[j]))
        bp1 = int(bp_range.split('-')[0])
        bp2 = int(bp_range.split('-')[1])
        signature = getSignature(bp1, bp_name, querymap, sign, realnmaps, evaluelist, bitscorelist, bps_no, complex_sv, nchrom, short, piden, mismatch, gap_ratio)
        if bp_name == 'S-Nov_Ins':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = breakpointfile[i].split('\t')[7].split(' ')[1].split(',')[int(iso_ins_count - iso_ins_counter)].split(':')[0]
            coord = breakpointfile[i].split('\t')[7].split(' ')[1].split(',')[int(iso_ins_count - iso_ins_counter)].split(':')[1]
            s_nov_ins_size = breakpointfile[i].split('\t')[7].split(' ')[0].split(',')[int(iso_ins_count - iso_ins_counter)].split('ns')[1]
            pair = '.'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(coord)
            uniqname = read_name + '~' + bp_uname
            if int(s_nov_ins_size) > cutoff:
                print str(read_name) + '\t' + str(bp2-10) + '\t' + str(bp2) + '\t' + str(bp_name + '_bp') + ' ' + str(s_nov_ins_size) + '\t' + str(chrom) + '\t' + str(coord) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            iso_ins_counter = iso_ins_counter - 1
        elif bp_name == 'E-Nov_Ins':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = breakpointfile[i].split('\t')[7].split(' ')[1].split(',')[int(iso_ins_count - iso_ins_counter)].split(':')[0]
            coord = breakpointfile[i].split('\t')[7].split(' ')[1].split(',')[int(iso_ins_count - iso_ins_counter)].split(':')[1]
            e_nov_ins_size = breakpointfile[i].split('\t')[7].split(' ')[0].split(',')[int(iso_ins_count - iso_ins_counter)].split('ns')[1]
            pair = '.'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(coord)
            uniqname = read_name + '~' + bp_uname
            if int(e_nov_ins_size) > cutoff:
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp1+10) + '\t' + str(bp_name + '_bp') + ' ' + str(e_nov_ins_size) + '\t' + str(chrom) + '\t' + str(coord) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            iso_ins_counter = iso_ins_counter - 1
        elif bp_name == 'Del':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = breakpointfile[i].split('\t')[6].split(' ')[1].split(',')[int(del_count - del_counter)].split(':')[0]
            coord1 = int(breakpointfile[i].split('\t')[6].split(' ')[1].split(',')[int(del_count - del_counter)].split(':')[1].split('-')[0])
            coord2 = int(breakpointfile[i].split('\t')[6].split(' ')[1].split(',')[int(del_count - del_counter)].split(':')[1].split('-')[1])
            del_size = breakpointfile[i].split('\t')[6].split(' ')[0].split(',')[int(del_count - del_counter)].split('l')[1]
            pair = '.'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(min(coord1,coord2)) + '-' + str(max(coord1,coord2))
            uniqname = read_name + '~' + bp_uname
            if int(del_size) > cutoff:
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(del_size) + '\t' + str(chrom) + '\t' + str(coord1) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(del_size) + '\t' + str(chrom) + '\t' + str(coord2) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            del_counter = del_counter - 1
        elif bp_name == 'Nov_Ins':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = breakpointfile[i].split('\t')[8].split(' ')[1].split(',')[int(ins_count - ins_counter)].split(':')[0]
            coord1 = int(breakpointfile[i].split('\t')[8].split(' ')[1].split(',')[int(ins_count - ins_counter)].split(':')[1].split('-')[0])
            coord2 = int(breakpointfile[i].split('\t')[8].split(' ')[1].split(',')[int(ins_count - ins_counter)].split(':')[1].split('-')[1])
            nov_ins_size = breakpointfile[i].split('\t')[8].split(' ')[0].split(',')[int(ins_count - ins_counter)].split('ns')[1]
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(min(coord1,coord2)) + '-' + str(max(coord1,coord2))
            uniqname = read_name + '~' + bp_uname
            if int(nov_ins_size) > cutoff:
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp1+10) + '\t' + str(bp_name) + ' ' + str(nov_ins_size) + '\t' + str(chrom) + '\t' + str(coord1) + '\t' + str(uniqidx) + '\t' + str('Nov1') + '\t' + str(uniqname) + '\t' + str(signature)
                print str(read_name) + '\t' + str(bp2-10) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(nov_ins_size) + '\t' + str(chrom) + '\t' + str(coord2) + '\t' + str(uniqidx) + '\t' + str('Nov2') + '\t' + str(uniqname) + '\t' + str(signature)
            ins_counter = ins_counter - 1
        elif bp_name == 'TDupl':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = breakpointfile[i].split('\t')[10].split(' ')[1].split(',')[int(tdup_count - tdup_counter)].split(':')[0]
            coord1 = int(breakpointfile[i].split('\t')[10].split(' ')[1].split(',')[int(tdup_count - tdup_counter)].split(':')[1].split('-')[0])
            coord2 = int(breakpointfile[i].split('\t')[10].split(' ')[1].split(',')[int(tdup_count - tdup_counter)].split(':')[1].split('-')[1])
            pair = '.'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(min(coord1,coord2)) + '-' + str(max(coord1,coord2))
            uniqname = read_name + '~' + bp_uname
            print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' 99.99' + '\t' + str(chrom) + '\t' + str(coord1) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' 99.99' + '\t' + str(chrom) + '\t' + str(coord2) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            tdup_counter = tdup_counter - 1
        elif bp_name == 'Inv':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = breakpointfile[i].split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[0]
            coord1 = int(breakpointfile[i].split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[1].split('-')[0])
            coord2 = int(breakpointfile[i].split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[1].split('-')[1])
            pair = '.'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(min(coord1,coord2)) + '-' + str(max(coord1,coord2))
            uniqname = read_name + '~' + bp_uname
            print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' 99.99' + '\t' + str(chrom) + '\t' + str(coord1) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' 99.99' + '\t' + str(chrom) + '\t' + str(coord2) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            inv_counter = inv_counter - 1
        elif bp_name == 'Inv(1)':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = breakpointfile[i].split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[0]
            coord1 = int(breakpointfile[i].split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[1].split('-')[0])
            coord2 = int(breakpointfile[i].split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[1].split('-')[1])
            inv_size = breakpointfile[i].split('\t')[11].split(' ')[0].split(',')[int(inv_count - inv_counter)].split('nv')[1]
            pair = 'Inv1'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(min(coord1,coord2)) + '-' + str(max(coord1,coord2))
            uniqname = read_name + '~' + bp_uname
            if int(inv_size) > cutoff:
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(inv_size) + '\t' + str(chrom) + '\t' + str(coord1) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(inv_size) + '\t' + str(chrom) + '\t' + str(coord2) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            inv_counter = inv_counter - 1
        elif bp_name == 'Inv(2)':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = breakpointfile[i].split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[0]
            coord1 = int(breakpointfile[i].split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[1].split('-')[0])
            coord2 = int(breakpointfile[i].split('\t')[11].split(' ')[1].split(',')[int(inv_count - inv_counter)].split(':')[1].split('-')[1])
            inv_size = breakpointfile[i].split('\t')[11].split(' ')[0].split(',')[int(inv_count - inv_counter)].split('nv')[1]
            pair = 'Inv2'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(min(coord1,coord2)) + '-' + str(max(coord1,coord2))
            uniqname = read_name + '~' + bp_uname
            if int(inv_size) > cutoff:
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(inv_size) + '\t' + str(chrom) + '\t' + str(coord1) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(inv_size) + '\t' + str(chrom) + '\t' + str(coord2) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            inv_counter = inv_counter - 1
        elif bp_name == 'Intra-Ins':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = breakpointfile[i].split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count - intra_ins_counter)].split(':')[0]
            coord1 = int(breakpointfile[i].split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count - intra_ins_counter)].split(':')[1].split('-')[0])
            coord2 = int(breakpointfile[i].split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count - intra_ins_counter)].split(':')[1].split('-')[1])
            pair = '.'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(min(coord1,coord2)) + '-' + str(max(coord1,coord2))
            uniqname = read_name + '~' + bp_uname
            print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' 99.99' + '\t' + str(chrom) + '\t' + str(coord1) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' 99.99' + '\t' + str(chrom) + '\t' + str(coord2) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            intra_ins_counter = intra_ins_counter - 1
        elif bp_name == 'Intra-Ins(1)':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = breakpointfile[i].split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count - intra_ins_counter)].split(':')[0]
            coord1 = int(breakpointfile[i].split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count - intra_ins_counter)].split(':')[1].split('-')[0])
            coord2 = int(breakpointfile[i].split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count - intra_ins_counter)].split(':')[1].split('-')[1])
            intrains_size = breakpointfile[i].split('\t')[12].split(' ')[0].split(',')[int(intra_ins_count - intra_ins_counter)].split(')')[1]
            pair = 'Intra1'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(min(coord1,coord2)) + '-' + str(max(coord1,coord2))
            uniqname = read_name + '~' + bp_uname
            if int(intrains_size) > cutoff:
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(intrains_size) + '\t' + str(chrom) + '\t' + str(coord1) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(intrains_size) + '\t' + str(chrom) + '\t' + str(coord2) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            intra_ins_counter = intra_ins_counter - 1
        elif bp_name == 'Intra-Ins(2)':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom = breakpointfile[i].split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count - intra_ins_counter)].split(':')[0]
            coord1 = int(breakpointfile[i].split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count - intra_ins_counter)].split(':')[1].split('-')[0])
            coord2 = int(breakpointfile[i].split('\t')[12].split(' ')[1].split(',')[int(intra_ins_count - intra_ins_counter)].split(':')[1].split('-')[1])
            intrains_size = breakpointfile[i].split('\t')[12].split(' ')[0].split(',')[int(intra_ins_count - intra_ins_counter)].split(')')[1]
            pair = 'Intra2'
            uniqidx = bp_uname + '~' + str(chrom) + ':' + str(min(coord1,coord2)) + '-' + str(max(coord1,coord2))
            uniqname = read_name + '~' + bp_uname
            if int(intrains_size) > cutoff:
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(intrains_size) + '\t' + str(chrom) + '\t' + str(coord1) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(intrains_size) + '\t' + str(chrom) + '\t' + str(coord2) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            intra_ins_counter = intra_ins_counter - 1
        elif bp_name == 'Inter-Ins(1)':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom1= breakpointfile[i].split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count - inter_ins_counter)].split('~')[0].split(':')[0]
            coord1_1 = breakpointfile[i].split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count - inter_ins_counter)].split('~')[0].split(':')[1]
            chrom2 = breakpointfile[i].split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count - inter_ins_counter)].split('~')[1].split(':')[0]
            coord2_1 = breakpointfile[i].split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count - inter_ins_counter)].split('~')[1].split(':')[1]
            interins_size = breakpointfile[i].split('\t')[13].split(' ')[0].split(',')[int(inter_ins_count - inter_ins_counter)].split(')')[1]
            chrmdict = {chrnum(chrom1):coord1_1, chrnum(chrom2):coord2_1}
            chrmlist = trysort(chrnum(chrom1),chrnum(chrom2)) #chrmlist = sorted([chrnum(chrom1), chrnum(chrom2)])
            pair = 'Inter1'
            uniqidx = bp_uname + '~' + str(revchrnum(chrmlist[0])) + ':' + str(chrmdict[chrmlist[0]]) + '~' + str(revchrnum(chrmlist[1])) + ':' + str(chrmdict[chrmlist[1]])
            uniqname = read_name + '~' + bp_uname
            if int(interins_size) > cutoff:
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(interins_size) + '\t' + str(revchrnum(chrmlist[0])) + '\t' + str(chrmdict[chrmlist[0]]) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(interins_size) + '\t' + str(revchrnum(chrmlist[1])) + '\t' + str(chrmdict[chrmlist[1]]) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            inter_ins_counter = inter_ins_counter - 1
        elif bp_name == 'Inter-Ins(2)':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom1= breakpointfile[i].split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count - inter_ins_counter)].split('~')[0].split(':')[0]
            coord1_1 = breakpointfile[i].split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count - inter_ins_counter)].split('~')[0].split(':')[1]
            chrom2 = breakpointfile[i].split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count - inter_ins_counter)].split('~')[1].split(':')[0]
            coord2_1 = breakpointfile[i].split('\t')[13].split(' ')[1].split(',')[int(inter_ins_count - inter_ins_counter)].split('~')[1].split(':')[1]
            interins_size = breakpointfile[i].split('\t')[13].split(' ')[0].split(',')[int(inter_ins_count - inter_ins_counter)].split(')')[1]
            chrmdict = {chrnum(chrom1):coord1_1, chrnum(chrom2):coord2_1}
            chrmlist = trysort(chrnum(chrom1),chrnum(chrom2)) #chrmlist = sorted([chrnum(chrom1), chrnum(chrom2)])
            pair = 'Inter2'
            uniqidx = bp_uname + '~' + str(revchrnum(chrmlist[0])) + ':' + str(chrmdict[chrmlist[0]]) + '~' + str(revchrnum(chrmlist[1])) + ':' + str(chrmdict[chrmlist[1]])
            uniqname = read_name + '~' + bp_uname
            if int(interins_size) > cutoff:
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(interins_size) + '\t' + str(revchrnum(chrmlist[0])) + '\t' + str(chrmdict[chrmlist[0]]) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
                print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' ' + str(interins_size) + '\t' + str(revchrnum(chrmlist[1])) + '\t' + str(chrmdict[chrmlist[1]]) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            inter_ins_counter = inter_ins_counter - 1
        elif bp_name == 'InterTx':
            bp_uname = "".join(random.sample(ran, ran_len))
            chrom1= breakpointfile[i].split('\t')[14].split(' ')[1].split(',')[int(inter_tx_count - inter_tx_counter)].split('~')[0].split(':')[0]
            coord1_1 = breakpointfile[i].split('\t')[14].split(' ')[1].split(',')[int(inter_tx_count - inter_tx_counter)].split('~')[0].split(':')[1]
            chrom2 = breakpointfile[i].split('\t')[14].split(' ')[1].split(',')[int(inter_tx_count - inter_tx_counter)].split('~')[1].split(':')[0]
            coord2_1 = breakpointfile[i].split('\t')[14].split(' ')[1].split(',')[int(inter_tx_count - inter_tx_counter)].split('~')[1].split(':')[1]
            chrmdict = {chrnum(chrom1):coord1_1, chrnum(chrom2):coord2_1}
            chrmlist = trysort(chrnum(chrom1),chrnum(chrom2)) #chrmlist = sorted([chrnum(chrom1), chrnum(chrom2)])
            pair = '.'
            uniqidx = bp_uname + '~' + str(revchrnum(chrmlist[0])) + ':' + str(chrmdict[chrmlist[0]]) + '~' + str(revchrnum(chrmlist[1])) + ':' + str(chrmdict[chrmlist[1]])
            uniqname = read_name + '~' + bp_uname
            print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' 99.99' + '\t' + str(revchrnum(chrmlist[0])) + '\t' + str(chrmdict[chrmlist[0]]) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            print str(read_name) + '\t' + str(bp1) + '\t' + str(bp2) + '\t' + str(bp_name) + ' 99.99' + '\t' + str(revchrnum(chrmlist[1])) + '\t' + str(chrmdict[chrmlist[1]]) + '\t' + str(uniqidx) + '\t' + str(pair) + '\t' + str(uniqname) + '\t' + str(signature)
            inter_tx_counter = inter_tx_counter - 1
        else:
            sys.exit("ERROR @ line " + str(g) + ' Unidentifiable SV')

logging.info("Total SV reads " + str(g))
logging.info("Job Finished")