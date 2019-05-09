"""
This script detects and characterizes SVs.

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
if len(argv)<5 or len(argv)>=6:
    sys.exit("Usage: python nv_hsblast_SV_detector.py hsblast.tsv SV%threshold[float] path_to_gap_file path_to_log_file > output.tsv")

#Assign variable to inputs
path_file = argv[1]
SVthres = float(argv[2])
hg38gap_path = argv[3]
log_path = argv[4]

#Set SV size cutoff to 19bp
ncutoff = 19

#Define Config of logging
logging.basicConfig(filename=log_path, level=logging.DEBUG, format='%(asctime)s - %(message)s')

#Function to create hg38 gap dictionary
def makegapdict(hg38gap_path):
    global gapdict
    gapdict = {}
    if hg38gap_path != '0':
        rgapdata = open(hg38gap_path, 'r').read().split('\n')
        l = len(rgapdata)
        k = range(l-1)
        for i in k: #Define gap dict list
            gapdict[rgapdata[i].split('\t')[0]] = []
        for i in k: #Making gap dictionary
            gapdict[rgapdata[i].split('\t')[0]].append([int(rgapdata[i].split('\t')[1]), int(rgapdata[i].split('\t')[2])])
        logging.info('Gap dictionary successfully loaded')
        logging.info('Start SV Detection..')
        return gapdict
    else:
        logging.info('Gap dictionary not loaded and will not be used')

def main():
    #Read file and split data into read alignment queries
    rdata = open(path_file, 'r').read().split('\n')
    #Calculating number of queries
    l = len(rdata)
    #Set range according to number of queries
    k = range(l-1)
    #Temporary query counter, Overall query counter, Total Read counter, Passed Read Counter
    c = n = o = r = 0
    #Create last line dummy
    rdata[l-1] = 'dum\tdum\tdum\tdum\tdum\tdum\tdum\tdum\tdum\tdum\tdum'
    #Reset false counter
    false_count = 0
    tmpread = []
    gapdict = makegapdict(hg38gap_path)
    #Loop through each event in file
    for i in k:
        n += 1
        #Grouping alignments by read name
        if rdata[i].split('\t')[4] == rdata[i+1].split('\t')[4]:
            if float(rdata[i].split('\t')[9]) <= 1: #Added line for future evalue filtering
                tmpread.append(rdata[i])
                c += 1
        else:
            if float(rdata[i].split('\t')[9]) <= 1: #Added line for future evalue filtering
                tmpread.append(rdata[i])
                c += 1
            false = '' #Reset false positive flag to null
            o += 1
            if c == 0:
                continue
            #Check if single alignment covers >=SVthres of read length
            if c == 1:
                if float(tmpread[0].split('\t')[6])/int(tmpread[0].split('\t')[8]) >= SVthres: 
                    tmpread = [] #Empty temp file
                    c = 0 #Reset alignment count
                    false = 'False-Pos' #Label as False positve alignment
                    #false_count += 1 #Don't count towards false positive
                    continue
                #Check if single alignment is large enough to qualify
                if int(tmpread[0].split('\t')[6]) < 200: 
                    tmpread = [] #Empty temp file
                    c = 0 #Reset alignment count
                    false = 'False-Pos' #Label as False positve alignment
                    #false_count += 1 #Don't count towards false positive
                    continue
            actualc = tmpread[0].split('\t')[14].split('=')[1]
            #Throw away multi-mapper reads (not in use, in order to capture LINE)
            if int(actualc) >= 1000:
                tmpread = [] #Empty temp file
                c = 0 #Reset alignment count
                false = 'False-Pos' #Label as False positve alignment
                #false_count += 1 #Don't count towards false positive
                continue
            r += 1
            chromocollect = tmpread[0].split('\t')[15].split('c')[0]
            #Limit number of maps
            if c > 4:
                sortdictlen = {}
                for line in tmpread:
                    #sortdictlen[int(line.split('\t')[6])] = line
                    sortdictlen[line] = int(line.split('\t')[6])
                #tmpread2 = [value for (key, value) in sorted(sortdictlen.items(), reverse=True)] #Sorting according to alignment length
                tmpread2 = [key for (key, value) in sorted(sortdictlen.items(), key=lambda x:x[1], reverse=True)]
                tmpread3 = tmpread2[0:4] #Selecting top 4 longest alignments
                #Sort back according to query start
                sortdict = {}
                for h in tmpread3:
                    sortdict[int(h.split('\t')[5])] = h
                tmpread = [value for (key, value) in sorted(sortdict.items())]
                c = 4
            #Set ranges
            j = range(c)
            p = range(c-1)
            #Define SV flags
            evalue_total, bitscore_ratio_total, piden_total, mismatch_ratio_total, gap_ratio_total, query, subject, chromorder, qurygaps, complex_SV, Iso_nov_ins, Iso_ins_size, Iso_ins_range, sv_range, intra_ins, intra_ins_range, dele, del_size, del_range, tdup, tdup_range, nov_ins, ins_size, ins_range, inv, inv_range, complex_nov_ins, complex_ins_size, inter_ins, inter_ins_range, intertx, intertx_break = ([] for i in range(32))
            #Satellite elements counter, Short alignments counter, Ref chromosome counter, Strandness counter, Second bp markers
            a = f = q = s = g = gtx = 0
            sat = ''
            #Choose tentative reference chromosome and strandness
            refchr = tmpread[0].split('\t')[0]
            refstrand = tmpread[0].split('\t')[7]
            #Retrieving read length
            readlength = int(tmpread[0].split('\t')[8])
            #Define short length alignment (Arbitrary definition as <5% of readlength)
            t = int(readlength)*0.05
            for i in j: #Collecting alignment information and filtering
                #Scale evalue 100000*, round to 5 dec places
                evalue = round(float(tmpread[i].split('\t')[9])*100000, 5)
                #Calculate bitscore ratio
                bitscore_ratio = round(float(tmpread[i].split('\t')[10])/int(tmpread[i].split('\t')[6]), 5)
                #Collect alignment percentage identity
                piden = tmpread[i].split('\t')[11]
                #Calculate mismatch ratio
                mismatch_ratio = round(float(tmpread[i].split('\t')[12])/int(tmpread[i].split('\t')[6]), 5)
                #Calculate gap ratio
                gap_ratio = round(float(tmpread[i].split('\t')[13])/int(tmpread[i].split('\t')[6]), 5)
                #Collect all alignment info
                evalue_total.append(str(evalue))
                bitscore_ratio_total.append(str(bitscore_ratio))
                piden_total.append(str(piden))
                mismatch_ratio_total.append(str(mismatch_ratio))
                gap_ratio_total.append(str(gap_ratio))
                #Collecting query and subject information
                query.append([int(tmpread[i].split('\t')[5]), int(tmpread[i].split('\t')[5]) + int(tmpread[i].split('\t')[6])])
                if tmpread[i].split('\t')[7] == '+':
                    subject.append([int(tmpread[i].split('\t')[1]), int(tmpread[i].split('\t')[1]) + int(tmpread[i].split('\t')[2])])
                elif tmpread[i].split('\t')[7] == '-':
                    subject.append([int(tmpread[i].split('\t')[1]) + int(tmpread[i].split('\t')[2]), int(tmpread[i].split('\t')[1])])
                #Filtering out strict telomeric and centromere regions
                chrm = tmpread[i].split('\t')[0].strip()
                chromorder.append(chrm)
                subjstart = int(tmpread[i].split('\t')[1])
                subjend = int(tmpread[i].split('\t')[1]) + int(tmpread[i].split('\t')[2])
                if hg38gap_path != '0':
                    try:
                        for v in gapdict[chrm]:
                            if v[0] <= (subjstart or subjend) <= v[1]:
                                a += 1
                    except:
                        pass
                #Count short alignments (Arbitrary definition as <5% of readlength)
                if int(tmpread[i].split('\t')[6]) < t:
                    f += 1
                #Count support for ref chromosome and strandness
                if tmpread[i].split('\t')[0] == refchr:
                    q += 1
                if tmpread[i].split('\t')[7] == refstrand:
                    s += 1
            #Calculating query middle gaps length
            for i in p:
                qurygaps.append(query[i+1][0] - query[i][1])
            #Check for Satellite element
            if a == c: #All alignments lies within satellite region
                false = 'False-Pos'
            elif 0 < a:
                false = 'False-Pos'
                #sat = 'Sate(' + str((float(a)/c)*100) + '%)'
            #Check for short alignment lengths
            short = float(f)/c
            #Defining major strandness
            if float(s)/c >= 0.5:
                majstrand = refstrand
            else:
                if refstrand == '+':
                    majstrand = '-'
                elif refstrand == '-':
                    majstrand = '+'
            #Counting number of inversions
            if majstrand == refstrand:
                if (c-s) > 1:
                    complex_SV.append('Multi_inv')
            else:
                if s > 1:
                    complex_SV.append('Multi_inv')
            #Calculating start/end gaps length and evaluation based on 100bp or more than ncutoff bases
            startgap = int(tmpread[0].split('\t')[5])
            endgap = int(tmpread[int(c)-1].split('\t')[8]) - int(tmpread[int(c)-1].split('\t')[6]) - int(tmpread[int(c)-1].split('\t')[5])
            termlength = 100
            if (startgap > termlength  and startgap > ncutoff):
                Iso_nov_ins.append('S-Nov-Ins')
                Iso_ins_size.append('S-Nov-Ins' + str(startgap))
                Iso_ins_range.append(str(tmpread[0].split('\t')[0].strip()) + ':' + str(subject[0][0]))
                sv_range.append('S-Nov_Ins:1-' + str(query[0][0]+1))
            if (endgap > termlength and endgap > ncutoff):
                Iso_nov_ins.append('E-Nov-Ins')
                Iso_ins_size.append('E-Nov-Ins' + str(endgap))
                Iso_ins_range.append(str(tmpread[-1].split('\t')[0].strip()) + ':' + str(subject[-1][1]))
                sv_range.append('E-Nov_Ins:' + str(query[-1][1]-1) + '-' + str(readlength))
            #Begin scanning all alignments on read
            if false == '':
                for i in p:
                    sbp = 0
                    qurygap = query[i+1][0] - query[i][1]
                    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ If same chromosome @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                    if tmpread[i].split('\t')[0] == tmpread[i+1].split('\t')[0]:
                        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Same strandness (1) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                        if tmpread[i].split('\t')[7] == '+' and tmpread[i+1].split('\t')[7] == '+':
                            if i != int(g-1): #Prevent duplicate record
                                #Gapstudy
                                subjgap = subject[i+1][0] - subject[i][1]
                                if subjgap > ncutoff: #If gap is positive and more than ncutoff
                                    if float(subjgap)/(qurygap+1) >= 2: #Subject gap needs to be at least twice the size of query gap
                                        #First breakpoint detected
                                        #Detecting second breakpoint
                                        for u in range(i+2, c):
                                            if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and any(subject[i][1] < y < subject[i+1][0] for y in subject[u]):
                                                #Second breakpoint detected, insertion event
                                                sbp = 1
                                                g = int(u)
                                                break
                                        if sbp == 1:
                                            #Second breakpoint present
                                            intra_ins.append('Intra-Ins(2)' + str(query[g][0]-query[i][1]))
                                            intra_ins.append('Intra-Ins(2)' + str(query[g][0]-query[i][1])) #Duplicate record for Intra-Ins(1)
                                            intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(subject[i+1][0]))
                                            intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[g-1][1]) + '-' + str(subject[g][0]))
                                            sv_range.append('Intra-Ins(1):' + str(query[i][1]) + '-' + str(query[i+1][0]) + ',' + 'Intra-Ins(2):' + str(query[g-1][1]) + '-' + str(query[g][0]))
                                        elif sbp == 0:                                           
                                            dele.append('Del')
                                            del_size.append('Del' + str(subjgap))
                                            del_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(subject[i+1][0]))
                                            sv_range.append('Del:' + str(query[i][1]) + '-' + str(query[i+1][0]))
                                elif -100000 < subjgap < -20: #If gap is negative and smaller than 100000 (Arbitrary tandemDup size gauge)
                                    tdup.append('TDupl')
                                    tdup_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(subject[i+1][0]))
                                    sv_range.append('TDupl:' + str(query[i][1]) + '-' + str(query[i+1][0]))
                                elif subjgap <= -100000: #If gap is negative and larger than or equal to 100000
                                    #Detecting second breakpoint
                                    for u in range(i+2, c):
                                        if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and all(subject[i][1] < y for y in subject[u]): #May still contain large deletion
                                            sbp = 1
                                            g = int(u)
                                            break
                                    if sbp == 1:
                                        #Second breakpoint present
                                        intra_ins.append('Intra-Ins(2)' + str(query[g][0]-query[i][1]))
                                        intra_ins.append('Intra-Ins(2)' + str(query[g][0]-query[i][1])) #Duplicate record for Intra-Ins(1)
                                        intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(subject[i+1][0]))  
                                        intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[g-1][1]) + '-' + str(subject[g][0]))
                                        sv_range.append('Intra-Ins(1):' + str(query[i][1]) + '-' + str(query[i+1][0]) + ',' + 'Intra-Ins(2):' + str(query[g-1][1]) + '-' + str(query[g][0]))
                                    elif sbp == 0:
                                        intra_ins.append('Intra-Ins')
                                        intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(subject[i+1][0]))  
                                        sv_range.append('Intra-Ins:' + str(query[i][1]) + '-' + str(query[i+1][0]))
                                if qurygap > ncutoff and subjgap >= 0: #If query gap more than cutoff and subject gap is positive
                                    if float(qurygap)/(subjgap+1) >= 2: #Test for novel insertion ratio
                                        nov_ins.append('Nov-Ins')
                                        ins_size.append('Nov-Ins' + str(qurygap))
                                        ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(int(subject[i][1])+1))
                                        sv_range.append('Nov_Ins:' + str(query[i][1]) + '-' + str(query[i+1][0]))
                                elif qurygap > ncutoff and subjgap < 0: #If query gap more than cutoff and subject gap is negative                      
                                    nov_ins.append('Nov-Ins')
                                    ins_size.append('Nov-Ins' + str(qurygap))
                                    ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(int(subject[i][1])+1))
                                    sv_range.append('Nov_Ins:' + str(query[i][1]) + '-' + str(query[i+1][0]))
                        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Same strandness (2) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                        if tmpread[i].split('\t')[7] == '-' and tmpread[i+1].split('\t')[7] == '-':
                            if i != int(g-1): #Prevent duplicate record
                                #Gapstudy
                                subjgap = subject[i][1] - subject[i+1][0]
                                if subjgap > ncutoff: #Is gap positive value and more than ncutoff? 
                                    if float(subjgap)/(qurygap+1) >= 2: #Test for significant gap using ratio
                                        #First breakpoint detected
                                        #Detecting second breakpoint
                                        for u in range(i+2, c):
                                            if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and any(subject[i+1][0] < y < subject[i][1] for y in subject[u]):
                                                #Second breakpoint detected, insertion event
                                                sbp = 1
                                                g = int(u)
                                                break
                                        if sbp == 1:
                                            #Second breakpoint present
                                            intra_ins.append('Intra-Ins(2)' + str(query[g][0]-query[i][1]))
                                            intra_ins.append('Intra-Ins(2)' + str(query[g][0]-query[i][1])) #Duplicate record for Intra-Ins(1)
                                            intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(subject[i+1][0]))  
                                            intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[g-1][1]) + '-' + str(subject[g][0]))
                                            sv_range.append('Intra-Ins(1):' + str(query[i][1]) + '-' + str(query[i+1][0]) + ',' +'Intra-Ins(2):' + str(query[g-1][1]) + '-' + str(query[g][0]))
                                        elif sbp == 0:                                    
                                            dele.append('Del')
                                            del_size.append('Del' + str(subjgap))
                                            del_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(subject[i+1][0]))
                                            sv_range.append('Del:' + str(query[i][1]) + '-' + str(query[i+1][0]))
                                elif -100000 < subjgap < -20: #If gap is negative and smaller than 100000 (Arbitrary tandemDup size gauge)
                                    tdup.append('TDupl')
                                    tdup_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(subject[i+1][0]))
                                    sv_range.append('TDupl:' + str(query[i][1]) + '-' + str(query[i+1][0]))
                                elif subjgap <= -100000: #If gap is negative and larger than or equal to 100000
                                    #Detecting second breakpoint
                                    for u in range(i+2, c):
                                        if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and all(subject[i+1][0] > y for y in subject[u]): #May still contain large deletion
                                            sbp = 1
                                            g = int(u)
                                            break
                                    if sbp == 1:
                                        #Second breakpoint present
                                        intra_ins.append('Intra-Ins(2)' + str(query[g][0]-query[i][1]))
                                        intra_ins.append('Intra-Ins(2)' + str(query[g][0]-query[i][1])) #Duplicate record for Intra-Ins(1)
                                        intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(subject[i+1][0]))  
                                        intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[g-1][1]) + '-' + str(subject[g][0]))
                                        sv_range.append('Intra-Ins(1):' + str(query[i][1]) + '-' + str(query[i+1][0]) + ',' + 'Intra-Ins(2):' + str(query[g-1][1]) + '-' + str(query[g][0]))
                                    elif sbp == 0: 
                                        intra_ins.append('Intra-Ins')
                                        intra_ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(subject[i+1][0]))
                                        sv_range.append('Intra-Ins:' + str(query[i][1]) + '-' + str(query[i+1][0]))
                                if qurygap > ncutoff and subjgap >= 0: #If query gap more than cutoff and subject gap is positive
                                    if float(qurygap)/(subjgap+1) >= 2: #Test for novel insertion ratio
                                        nov_ins.append('Nov-Ins')
                                        ins_size.append('Nov-Ins' + str(qurygap))
                                        ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(int(subject[i][1])+1))
                                        sv_range.append('Nov_Ins:' + str(query[i][1]) + '-' + str(query[i+1][0]))
                                elif qurygap > ncutoff and subjgap < 0: #If query gap more than cutoff and subject gap is negative 
                                    nov_ins.append('Nov-Ins')
                                    ins_size.append('Nov-Ins' + str(qurygap))
                                    ins_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(int(subject[i][1])+1))
                                    sv_range.append('Nov_Ins:' + str(query[i][1]) + '-' + str(query[i+1][0]))
                        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Different strandness (1) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                        elif tmpread[i].split('\t')[7] == '+' and tmpread[i+1].split('\t')[7] == '-': #Inversion
                            if i != int(g-1):
                                for u in range(i+2, c):
                                    if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and any(subject[i][1] < y for y in subject[u]):
                                        #Second breakpoint detected, paired inversion event
                                        sbp = 1
                                        g = int(u)
                                        break
                                if sbp == 1:
                                    #Second breakpoint present
                                    inv.append('Pair_Inv' + str(query[g][0]-query[i][1]))
                                    inv.append('Pair_Inv' + str(query[g][0]-query[i][1])) #Duplicate record for Inv(1)
                                    inv_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(subject[i+1][0]))  
                                    inv_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[g-1][1]) + '-' + str(subject[g][0]))
                                    sv_range.append('Inv(1):' + str(query[i][1]) + '-' + str(query[i+1][0]) + ',' +'Inv(2):' + str(query[g-1][1]) + '-' + str(query[g][0]))
                                elif sbp == 0:
                                    inv.append('Unpair_Inv')
                                    inv_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(subject[i+1][0]))
                                    sv_range.append('Inv:' + str(query[i][1]) + '-' + str(query[i+1][0]))
                            if qurygap > ncutoff: #If query gap more than ncutoff
                                complex_nov_ins.append('Inv-Nov-Ins')
                                complex_SV.append('Inv-Nov_Ins')
                                complex_ins_size.append('Inv-Nov_Ins' + str(qurygap))
                        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Different Strandness(2) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                        elif tmpread[i].split('\t')[7] == '-' and tmpread[i+1].split('\t')[7] == '+': #Inversion
                            if i != int(g-1):
                                for u in range(i+2, c):
                                    if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and any(y < subject[i][1] for y in subject[u]):
                                        #Second breakpoint detected, paired inversion event
                                        sbp = 1
                                        g = int(u)
                                        break
                                if sbp == 1:
                                    #Second breakpoint present
                                    inv.append('Pair_Inv' + str(query[g][0]-query[i][1]))
                                    inv.append('Pair_Inv' + str(query[g][0]-query[i][1])) #Duplicate record for Inv(1)
                                    inv_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(subject[i+1][0]))  
                                    inv_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[g-1][1]) + '-' + str(subject[g][0]))
                                    sv_range.append('Inv(1):' + str(query[i][1]) + '-' + str(query[i+1][0]) + ',' +'Inv(2):' + str(query[g-1][1]) + '-' + str(query[g][0]))
                                elif sbp == 0:
                                    inv.append('Unpair_Inv')
                                    inv_range.append(str(tmpread[i].split('\t')[0].strip()) + ':' + str(subject[i][1]) + '-' + str(subject[i+1][0]))
                                    sv_range.append('Inv:' + str(query[i][1]) + '-' + str(query[i+1][0]))
                            if qurygap > ncutoff: #If query gap more than ncutoff
                                complex_nov_ins.append('Inv-Nov-Ins')
                                complex_SV.append('Inv-Nov_Ins')
                                complex_ins_size.append('Inv-Nov_Ins' + str(qurygap))
                    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Different chromosome @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                    elif tmpread[i].split('\t')[0] != tmpread[i+1].split('\t')[0]:
                        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Strandness+ $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                        if tmpread[i].split('\t')[7] == '+':
                            if i != int(gtx-1):
                                for u in range(i+2, c):
                                    if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and any(subject[i][1] < y for y in subject[u]): #Detecting second breakpoint                                                               
                                        #Second breakpoint detected, insertion event
                                        sbp = 1
                                        gtx = int(u)
                                        break
                                if sbp == 1:
                                    #Second breakpoint present
                                    inter_ins.append('Inter-Ins(2)' + str(query[gtx][0]-query[i][1]))
                                    inter_ins.append('Inter-Ins(2)' + str(query[gtx][0]-query[i][1])) #Duplicate record for Inter-Ins(1)
                                    inter_ins_range.append(str(tmpread[i].split('\t')[0]) + ':' + str(subject[i][1]) + '~' + str(tmpread[i+1].split('\t')[0]) + ':' + str(subject[i+1][0]))
                                    inter_ins_range.append(str(tmpread[gtx-1].split('\t')[0]) + ':' + str(subject[gtx-1][1]) + '~' + str(tmpread[gtx].split('\t')[0]) + ':' + str(subject[gtx][0]))
                                    sv_range.append('Inter-Ins(1):' + str(query[i][1]) + '-' + str(query[i+1][0]) + ',' + 'Inter-Ins(2):' + str(query[gtx-1][1]) + '-' + str(query[gtx][0]))
                                elif sbp == 0:                                            
                                    intertx.append('InterTx')
                                    intertx_break.append(str(tmpread[i].split('\t')[0]) + ':' + str(subject[i][1]) + '~' + str(tmpread[i+1].split('\t')[0]) + ':' + str(subject[i+1][0]))
                                    sv_range.append('InterTx:' + str(query[i][1]) + '-' + str(query[i+1][0]))
                                if qurygap > ncutoff: #If query gap more than ncutoff
                                    complex_nov_ins.append('Inter-Nov-Ins')
                                    complex_SV.append('Inter-Nov_Ins')
                                    complex_ins_size.append('Inter-Nov_Ins' + str(qurygap))
                        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Strandness- $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                        if tmpread[i].split('\t')[7] == '-':
                            if i != int(gtx-1):
                                for u in range(i+2, c):
                                    if tmpread[i].split('\t')[0] == tmpread[u].split('\t')[0] and tmpread[i].split('\t')[7] == tmpread[u].split('\t')[7] and any(y < subject[i][1] for y in subject[u]): #Detecting second breakpoint                                                               
                                        #Second breakpoint detected, insertion event
                                        sbp = 1
                                        gtx = int(u)
                                        break
                                if sbp == 1:
                                    #Second breakpoint present
                                    inter_ins.append('Inter-Ins(2)' + str(query[gtx][0]-query[i][1]))
                                    inter_ins.append('Inter-Ins(2)' + str(query[gtx][0]-query[i][1])) #Duplicate record for Inter-Ins(1)
                                    inter_ins_range.append(str(tmpread[i].split('\t')[0]) + ':' + str(subject[i][1]) + '~' + str(tmpread[i+1].split('\t')[0]) + ':' + str(subject[i+1][0]))
                                    inter_ins_range.append(str(tmpread[gtx-1].split('\t')[0]) + ':' + str(subject[gtx-1][1]) + '~' + str(tmpread[gtx].split('\t')[0]) + ':' + str(subject[gtx][0]))
                                    sv_range.append('Inter-Ins(1):' + str(query[i][1]) + '-' + str(query[i+1][0]) + ',' + 'Inter-Ins(2):' + str(query[gtx-1][1]) + '-' + str(query[gtx][0]))
                                elif sbp == 0:                                            
                                    intertx.append('InterTx')
                                    intertx_break.append(str(tmpread[i].split('\t')[0]) + ':' + str(subject[i][1]) + '~' + str(tmpread[i+1].split('\t')[0]) + ':' + str(subject[i+1][0]))
                                    sv_range.append('InterTx:' + str(query[i][1]) + '-' + str(query[i+1][0]))
                                if qurygap > ncutoff: #If query gap more than ncutoff
                                    complex_nov_ins.append('Inter-Nov-Ins')
                                    complex_SV.append('Inter-Nov_Ins')
                                    complex_ins_size.append('Inter-Nov_Ins' + str(qurygap))
            if dele == [] and Iso_nov_ins == [] and nov_ins == [] and complex_nov_ins == [] and tdup == [] and inv == [] and intra_ins == [] and inter_ins == [] and intertx == [] and complex_SV == []:
                false = 'False-Pos'
            Total_SV_complex = ' '.join(dele) + ' ' + ' '.join(Iso_nov_ins) + ' ' + ' '.join(nov_ins) + ' ' + ' '.join(tdup) + ' ' + ' '.join(inv) + ' ' + ' '.join(intra_ins) + ' ' + ' '.join(inter_ins) + ' ' + ' '.join(intertx)
            if false == '':
                querypercent = []
                qurygappercent = []
                fullquerypercent = []
                for i in p: #Calculating gap percentages
                    qurygappercent.append(str((float(qurygaps[i])/readlength)*100))
                startgappercent = str((float(startgap)/readlength)*100)
                qurygappercent.append(str((float(endgap)/readlength)*100))
                for i in j:
                    querypercent.append(str((float(tmpread[i].split('\t')[6])/readlength)*100)) #Calculating query map percentage
                fullquerypercent.append('%.1f' % float(startgappercent) + '%')
                for i in j:
                    fullquerypercent.append('(' + '%.1f' % float(querypercent[i]) + '%)')
                    fullquerypercent.append('%.1f' % float(qurygappercent[i]) + '%')
                print str(tmpread[i].split('\t')[4].strip()) + '\t' + str(readlength) + '\t' + str(chromocollect) + '\t' + str(actualc) + ' maps' + '\t' + ','.join(fullquerypercent) + '\t' + ','.join(evalue_total) + '\t' + ','.join(bitscore_ratio_total) + '\t' + ','.join(del_size) + ' ' + ','.join(del_range) + '\t' + ','.join(Iso_ins_size) + ' ' + ','.join(Iso_ins_range) + '\t' + ','.join(ins_size) + ' ' + ','.join(ins_range) + '\t' + ','.join(complex_ins_size) + '\t' + ','.join(tdup) + ' ' + ','.join(tdup_range) + '\t' + ','.join(inv) + ' ' + ','.join(inv_range) + '\t' + ','.join(intra_ins) + ' ' + ','.join(intra_ins_range) + '\t' + ','.join(inter_ins) + ' ' + ','.join(inter_ins_range) + '\t' + ','.join(intertx) + ' ' + ','.join(intertx_break) + '\t' + majstrand + '\t' + str(short) + '\t' + sat + '\t' + ','.join(complex_SV) + '\t' + ", ".join(Total_SV_complex.split()) + '\t' + ','.join(sv_range) + '\t' + str(chromorder) + '\t' + ''.join(str(subject)) + '\t' + str(query) + '\t' + ','.join(piden_total) + '\t' + ','.join(mismatch_ratio_total) + '\t' + ','.join(gap_ratio_total)
            elif false == 'False-Pos':
                false_count = false_count + 1
            tmpread = []
            g = 0
            gtx = 0
            c = 0
    logging.info("Total Reads " + str(o))
    logging.info("Number of Passed Reads " + str(r))
    logging.info("Number of SV Reads " + str(r-false_count))
    logging.info("Number of False-Positive Reads " + str(false_count))
    logging.info("Job Finished")

if __name__ == '__main__':
    main()

#Ouput columns 1=readname, 2=readlength, 3=total_chromosome(s), 4=No._of_maps, 5=query_signature, 6=evalue, 7=bitscore, 8=deletion, 9=Isolated_insertion, 10=Insertion, 11=complex_insertion, 12=tandemDupl, 13=inversion, 14=intra_insertion, 15=inter_insertion, 16=inter_tx, 17=majstrand, 18=short_frag_perc, 19=satellite_perc, 20=complex_SV, 21=Total_SV_complex, 22=sv_range, 23=chromosome_order, 24=subject_string, 25=query_string, 26=%identity, 27=mismatch, 28=gap_ratio
