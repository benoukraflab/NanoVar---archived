"""
This script clusters long reads around each SV breakend.

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
if len(argv)<6 or len(argv)>=7:
    sys.exit("Usage: python nv_lr_overlap.py parse.tsv hsblast.tsv Buffer_value[INT] bedtoolspath nsplit[INT] > output.tsv")

#Assign variable to inputs
parse_file = argv[1]
hsblast_file = argv[2]
buf = int(argv[3])
bedpath = argv[4]
nsplit = int(argv[5])

#Function to convert non-digits to zero
def alpha(n):
    if str(n).isdigit():
        return int(n)
    else:
        return int(0)

#Function to convert . to zeroq
def dotter(dot):
    if dot == '.':
        return float(0)
    else:
        return float(dot)

#Function to set minimum value to 0
def mz(value):
    if value < 0:
        return 0
    else:
        return value

#Function to collect information of each breakpoint entry (e.g. chromosome, range, SV size, SV type)
def rangecollect(x):
    namerepeat = {}
    rangedict = {}
    svsizedict = {}
    DNNdict = {}
    infodict = {}
    bed1 = []
    bed2 = []
    for i in x:
        try:
            infodict[i.split('\t')[8]].append('\t'.join(i.split('\t')[0:10]))
        except:
            infodict[i.split('\t')[8]] = []
            infodict[i.split('\t')[8]].append('\t'.join(i.split('\t')[0:10]))
        try:
            namerepeat[i.split('\t')[8]] == 1
        except:
            namerepeat[i.split('\t')[8]] = 1
            rnameidx = i.split('\t')[8]
            svtype = str(svtypecorrector(i.split('\t')[3].split(' ')[0]))
            #idx = ''
            chm1 = ''
            l = ''
            r = ''
            chm2 = ''
            if len(i.split('\t')[6].split('~')) == 2: # Not Inter translocation
                if len(i.split('\t')[6].split('~')[1].split(':')[1].split('-')) == 1:
                    #idx = i.split('\t')[6].split('~')[0]
                    chm1 = i.split('\t')[6].split('~')[1].split(':')[0]
                    l = int(i.split('\t')[6].split('~')[1].split(':')[1])
                    #rangedict[rnameidx + '-l'] = [chm1, l - buf, l + buf, svtype]
                    svsizedict[rnameidx] = alpha(i.split('\t')[3].split(' ')[-1])
                    try:
                        DNNdict[rnameidx] = dotter(i.split('\t')[10])
                    except:
                        pass
                    bed1.append(chm1 + '\t' + str(l) + '\t' + str(l+1) + '\t' + rnameidx + '-l' + '\t' + svtype)
                    bed2.append(chm1 + '\t' + str(mz(l-buf)) + '\t' + str(l+buf) + '\t' + rnameidx + '-l' + '\t' + svtype)
                elif len(i.split('\t')[6].split('~')[1].split(':')[1].split('-')) == 2:
                    #idx = i.split('\t')[6].split('~')[0]
                    chm1 = i.split('\t')[6].split('~')[1].split(':')[0]
                    l = int(i.split('\t')[6].split('~')[1].split(':')[1].split('-')[0])
                    r = int(i.split('\t')[6].split('~')[1].split(':')[1].split('-')[1])
                    #rangedict[rnameidx + '-l'] = [chm1, l - buf, l + buf, svtype]
                    #rangedict[rnameidx + '-r'] = [chm1, r - buf, r + buf, svtype]
                    svsizedict[rnameidx] = alpha(i.split('\t')[3].split(' ')[-1])
                    try:
                        DNNdict[rnameidx] = dotter(i.split('\t')[10])
                    except:
                        pass
                    bed1.append(chm1 + '\t' + str(l) + '\t' + str(l+1) + '\t' + rnameidx + '-l' + '\t' + svtype)
                    bed1.append(chm1 + '\t' + str(r) + '\t' + str(r+1) + '\t' + rnameidx + '-r' + '\t' + svtype)
                    bed2.append(chm1 + '\t' + str(mz(l-buf)) + '\t' + str(l+buf) + '\t' + rnameidx + '-l' + '\t' + svtype)
                    bed2.append(chm1 + '\t' + str(mz(r-buf)) + '\t' + str(r+buf) + '\t' + rnameidx + '-r' + '\t' + svtype)
            elif len(i.split('\t')[6].split('~')) == 3: # Inter translocation
                #idx = i.split('\t')[6].split('~')[0]
                chm1 = i.split('\t')[6].split('~')[1].split(':')[0]
                chm2 = i.split('\t')[6].split('~')[2].split(':')[0]
                l = int(i.split('\t')[6].split('~')[1].split(':')[1])
                r = int(i.split('\t')[6].split('~')[2].split(':')[1])
                #rangedict[rnameidx + '-l'] = [chm1, l - buf, l + buf, svtype]
                #rangedict[rnameidx + '-r'] = [chm2, r - buf, r + buf, svtype]
                svsizedict[rnameidx] = alpha(i.split('\t')[3].split(' ')[-1])
                try:
                    DNNdict[rnameidx] = dotter(i.split('\t')[10])
                except:
                    pass
                bed1.append(chm1 + '\t' + str(l) + '\t' + str(l+1) + '\t' + rnameidx + '-l' + '\t' + svtype)
                bed1.append(chm2 + '\t' + str(r) + '\t' + str(r+1) + '\t' + rnameidx + '-r' + '\t' + svtype)
                bed2.append(chm1 + '\t' + str(mz(l-buf)) + '\t' + str(l+buf) + '\t' + rnameidx + '-l' + '\t' + svtype)
                bed2.append(chm2 + '\t' + str(mz(r-buf)) + '\t' + str(r+buf) + '\t' + rnameidx + '-r' + '\t' + svtype)
    return rangedict,svsizedict,DNNdict,infodict,bed1,bed2

#Function to standardize names of SV types
def svtypecorrector(svtype):
    if svtype == 'S-Nov_Ins_bp':
        return 'bp_Nov_Ins'
    elif svtype == 'E-Nov_Ins_bp':
        return 'bp_Nov_Ins'
    elif svtype == 'Inter-Ins(1)':
        return 'Inter-Ins'
    elif svtype == 'Inter-Ins(2)':
        return 'Inter-Ins'
    elif svtype == 'InterTx':
        return 'Inter'
    elif svtype == 'Inv(1)':
        return 'Inv2'
    elif svtype == 'Inv(2)':
        return 'Inv2'
    elif svtype == 'Intra-Ins(1)':
        return 'Intra-Ins2'
    elif svtype == 'Intra-Ins(2)':
        return 'Intra-Ins2'
    else:
        return svtype

def svgroup(svclass):
    if svclass == 'bp_Nov_Ins':
        return ['bp_Nov_Ins', 'Nov_Ins', 'Inter-Ins', 'Inter', 'Inv2', 'Inv', 'Intra-Ins2', 'Intra-Ins', 'Del', 'TDupl']
    elif svclass == 'Nov_Ins':
        return ['bp_Nov_Ins', 'Nov_Ins']
    elif svclass == 'Inter-Ins':
        return ['bp_Nov_Ins', 'Inter-Ins', 'Inter']
    elif svclass == 'Inter':
        return ['bp_Nov_Ins', 'Inter-Ins', 'Inter']
    elif svclass == 'Inv2':
        return ['bp_Nov_Ins', 'Inv2', 'Inv']
    elif svclass == 'Inv':
        return ['bp_Nov_Ins', 'Inv2', 'Inv']
    elif svclass == 'Intra-Ins2':
        return ['bp_Nov_Ins', 'Intra-Ins2', 'Intra-Ins']
    elif svclass == 'Intra-Ins':
        return ['bp_Nov_Ins', 'Intra-Ins2', 'Intra-Ins']
    elif svclass == 'Del':
        return ['bp_Nov_Ins', 'Del']
    elif svclass == 'TDupl':
        return ['bp_Nov_Ins', 'TDupl']

#Function to connect overlaping overlaping breakpoint entries
def intersection(intersect):
    connectdict = {}
    classdict = {}
    for i in bed1:
        connectdict[i.split('\t')[3]] = []
        classdict[i.split('\t')[3][0:-2]] = i.split('\t')[4]
    for i in intersect:
        if i.split('\t')[4] in svgroup(i.split('\t')[9]):
            connectdict[i.split('\t')[8]].append(i.split('\t')[3][0:-2])
    return connectdict,classdict

#Remove unique identifier from read name
def remove_uniq(x):
    connectdictname = {}
    for i in x:
        connectdictname[i[:-8]] = []
    return connectdictname

def intersection2(intersect2, x, nsplit):
    svnormalconnect = {}
    svnormalcov = {}
    for key in x:
        svnormalconnect[key] = []
    for i in intersect2:
        if not connectdictname.has_key(i.split('\t')[3]):
            if classdict[i.split('\t')[7]] == 'TDupl': #If Dup, make sure normal read covers both breakpoints of Dup
                for n in infodict[i.split('\t')[7]]:
                    if int(i.split('\t')[1]) <= int(infodict[i.split('\t')[7]][0].split('\t')[5]) <= int(i.split('\t')[2]) and int(i.split('\t')[1]) <= int(infodict[i.split('\t')[7]][1].split('\t')[5]) <= int(i.split('\t')[2]):
                        svnormalconnect[i.split('\t')[7]].append(i.split('\t')[3])
            else:
                svnormalconnect[i.split('\t')[7]].append(i.split('\t')[3])
    for key in svnormalconnect:
        svnormalcov[key] = float(len(set(svnormalconnect[key])))/nsplit
    return svnormalcov

#Ordering classes
def classranker(x):
    tier3sv = ['Inv2', 'Inter-Ins', 'Intra-Ins2']
    tier2sv = ['Inv', 'Del', 'TDupl', 'Nov_Ins', 'Intra-Ins', 'Inter']
    tier1sv = ['bp_Nov_Ins']
    tier1 = []
    tier2 = []
    tier3 = []
    for keys in x:
        if x[keys] in tier1sv:
            tier1.append(keys)
        elif x[keys] in tier2sv:
            tier2.append(keys)
        elif x[keys] in tier3sv:
            tier3.append(keys)
        else:
            print 'sv class error'
            sys.exit()
    return tier1,tier2,tier3

def connector(x):
    namerepeat = {}
    readdictA = {}
    readdictB = {}
    alltier = tier3 + tier2 + tier1
    for keys in alltier:
        if not namerepeat.has_key(keys):
            namerepeat[keys] = 1
            readdictA[keys] = []
            readdictB[keys] = []
            readdictB[keys].append(classdict[keys])
            for i in x[keys + '-l']:
                if not namerepeat.has_key(i):
                    if x.has_key(keys + '-r'):
                        if i in x[keys + '-r']:
                            namerepeat[i] = 1
                            readdictA[keys].append(i)
                            readdictB[keys].append(classdict[i])
                        else:
                            if i in tier1: #Only for single-bp-end reads
                                namerepeat[i] = 1
                                readdictA[keys].append(i)
                                readdictB[keys].append(classdict[i])
                    else:
                        namerepeat[i] = 1
                        readdictA[keys].append(i)
                        readdictB[keys].append(classdict[i])
            if x.has_key(keys + '-r'):
                for i in x[keys + '-r']:
                    if not namerepeat.has_key(i):
                        if i in tier1: #Only for single-bp-end reads
                            namerepeat[i] = 1
                            readdictA[keys].append(i)
                            readdictB[keys].append(classdict[i])
            for i in readdictA[keys]:
                for k in x[i + '-l']:
                    if not namerepeat.has_key(k):
                        if x.has_key(i + '-r'):
                            if k in x[i + '-r']:
                                namerepeat[k] = 1
                                readdictA[keys].append(k)
                                readdictB[keys].append(classdict[k])
                        else:
                            pass #Ignoring fellows from single-bp-end reads
    return readdictA,readdictB

#Function to rank and filter best sv type
def classupdate(readdictB):
    readdictC = {}
    d = {}
    tier3sv = ['Inv2', 'Inter-Ins', 'Intra-Ins2']
    tier2sv = ['Inv', 'Del', 'TDupl', 'Nov_Ins', 'Intra-Ins', 'Inter']
    tier1sv = ['bp_Nov_Ins']
    for keys in readdictB:
        for i in readdictB[keys]:
            if not d.has_key(i):
                d[i] = 1
            else:
                d[i] += 1
        if any(y in tier3sv for y in d):
            for s in tier2sv:
                try:
                    del d[s]
                except:
                    pass
            try:
                del d["bp_Nov_Ins"]
            except:
                pass
        elif any(y in tier2sv for y in d):
            try:
                del d["bp_Nov_Ins"]
            except:
                pass
        elif any(y in tier1sv for y in d):
            pass
        else:
            sys.exit("ERROR: UNKNOWN SV CLASS")
            print "Uknown class"
            print keys
            print i
            print d
            break
        mainsv = [key for (key, value) in sorted(d.items(), key=lambda x:x[1], reverse=True)][0]
        readdictC[keys] = mainsv
        d = {}
    return readdictC

#Function to find the lead breakpoint entry based on main sv type first come across
def bestread(x):
    newdict = {}
    for keys in x:
        main = ''
        mainsv = readdictC[keys]
        if classdict[keys] == mainsv:
            main = keys
            newdict[main] = []
            for k in x[main]:
                newdict[main].append(k)
        else:
            for k in x[keys]:
                if classdict[k] == mainsv:
                    main = k
                    break
            newdict[main] = []
            newdict[main].append(keys)
            for r in x[keys]:
                if r != main:
                    newdict[main].append(r)
    return newdict

#Function to generate bed file from hsblast.tsv
def normalbed(hsblast):
    totalbed = []
    with open(hsblast) as hs:
        for line in hs:
            coord1 = int(line.split('\t')[1]) + 400
            coord2 = int(line.split('\t')[1]) + int(line.split('\t')[2]) - 400
            if coord2 - coord1 > 0:
                totalbed.append(line.split('\t')[0] + '\t' + str(coord1) + '\t' + str(coord2) + '\t' + line.split('\t')[4])
            else:
                pass
    return totalbed

def svbed(x):
    totalsvbed = []
    for key in x:
        for n in infodict[key]:
            totalsvbed.append('\t'.join(n.split('\t')[4:6]) + '\t' + str(int(n.split('\t')[5]) + 1) + '\t' + key)
    return totalsvbed

#Function to count total number of unique reads in an overlaping breakpoint entry
def countcov(key, dict):
    readdict = {}
    readdict[key[:-6]] = 1
    for i in dict[key]:
        readdict[i[:-6]] = 1
    return len(readdict)

#Function to calculate the average short read coverage of an overlaping breakpoint entry
def avgDNN(key, value):
    dnn = 0
    num = len(value) + 1
    read = []
    dnn = dnn + DNNdict[key]
    for k in value:
        dnn = dnn + DNNdict[k]
    avg = float(dnn)/num
    return avg

#Function to parse lead overlaping breakpoint entry and long read coverage into output
def arrange(x):
    output = []
    for key in x:
        lcov = countcov(key, x) #len(x[key]) + 1
        if lcov == 1:
            for n in infodict[key]:
                output.append(n + '\t' + str(lcov) + '\t.\t' + str(svnormalcov[key]) + '\n')
        elif lcov > 1:
            for n in infodict[key]:
                output.append(n + '\t' + str(lcov) + '\t' + ','.join(x[key]) + '\t' + str(svnormalcov[key]) + '\n')
    return output

#Function to parse lead overlaping breakpoint entry, long read coverage and average DNN score into output
def arrangeDNN(x):
    output = []
    for key in x:
        DNN = avgDNN(key, x[key])
        lcov = countcov(key, x) #len(x[key]) + 1
        if lcov == 1:
            for n in infodict[key]:
                output.append(n + '\t' + str(lcov) + '\t.\t' + str(svnormalcov[key]) + '\t' + str(DNN) + '\n')
        elif lcov > 1:
            for n in infodict[key]:
                output.append(n + '\t' + str(lcov) + '\t' + ','.join(x[key]) + '\t' + str(svnormalcov[key]) + '\t' + str(DNN) + '\n')
    return output

def main():
    global svsizedict
    global DNNdict
    global infodict
    global bed1
    global connectdict
    global classdict
    global tier1
    global tier2
    global tier3
    global readdictC
    global connectdictname
    global svnormalcov
    parsefile = open(parse_file, 'r').read().splitlines()
    #Test for DNN scoring
    try:
        if float(parsefile[0].split('\t')[10]) > 0:
            DNN = 1
    except:
        DNN = 0
    rangedict,svsizedict,DNNdict,infodict,bed1,bed2 = rangecollect(parsefile)
    bedout1 = open('bed1', 'w')
    bedout1.write('\n'.join(bed1))
    bedout1.close()
    bedout2 = open('bed2', 'w')
    bedout2.write('\n'.join(bed2))
    bedout2.close()
    os.system(bedpath + ' sort -i bed1 > bed1.sort')
    os.system(bedpath + ' sort -i bed2 > bed2.sort')
    os.system(bedpath + " intersect -wa -wb -a bed1.sort -b bed2.sort | awk -F'\t' '{if ($4 != $9) print $0}' > intersect.txt")
    intersect = open('intersect.txt', 'r').read().splitlines()
    connectdict,classdict = intersection(intersect)
    tier1,tier2,tier3 = classranker(classdict)
    readdictA,readdictB = connector(connectdict)
    readdictC = classupdate(readdictB)
    newdict = bestread(readdictA)
    totalbed = normalbed(hsblast_file)
    bedout3 = open('bed3', 'w')
    bedout3.write('\n'.join(totalbed))
    bedout3.close()
    os.system(bedpath + ' sort -i bed3 > bed3.sort')
    totalsvbed = svbed(newdict)
    bedout4 = open('bed4', 'w')
    bedout4.write('\n'.join(totalsvbed))
    bedout4.close()
    os.system(bedpath + ' sort -i bed4 > bed4.sort')
    os.system(bedpath + " intersect -wa -wb -a bed3.sort -b bed4.sort | awk -F'\t' '{if ($4 != $8) print $0}' > intersect2.txt")
    connectdictname = remove_uniq(connectdict) #Remove unique identifier from read name
    intersect2 = open('intersect2.txt', 'r').read().splitlines()
    svnormalcov = intersection2(intersect2, newdict, nsplit)
    if DNN == 0:
        output = arrange(newdict)
    elif DNN == 1:
        output = arrangeDNN(newdict)
    output[-1] = output[-1].strip('\n')
    print ''.join(output)

if __name__ == '__main__':
    main()

#Note: Some bps seen in parse_file would be missing in overlap_file because they are overlapped by another bp within the same read
