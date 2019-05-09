"""
This script carries out ANN model inferencing using Keras library.

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
import numpy as np
import random

if len(argv)<7 or len(argv)>=8:
    sys.exit("Usage: python keras_ANN_inference.py overlap.tsv parse.tsv model.h5 maxovl splitcov [parse or overlap] > parse/overlap.NN.tsv")
overlap_file = argv[1]
parse_file = argv[2]
modelh5 = argv[3]
maxovl = float(argv[4]) #not in use
splitcov = float(argv[5]) #not in use
mode = argv[6]

from keras.models import Sequential
from keras.layers import Dense, Activation
from keras.layers import Dropout
from keras.models import load_model
from keras.optimizers import SGD
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

#For each breakpoint line, assign the long read coverage and normal cov ratio
def getoverlap(overlap_data):
    overlapdict = {}
    normalcovratio = {}
    classdict = {}
    for i in overlap_data:
        if int(i.split('\t')[10]) >= maxovl: #Not in use due to pre-filter of maxovl
            overlapdict[i.split('\t')[8]] = maxovl #Not in use
            normalcovratio[i.split('\t')[8]] = float(i.split('\t')[12])/splitcov #Not in use
            classdict[i.split('\t')[8]] = [i.split('\t')[3].split(' ')[0], i.split('\t')[3].split(' ')[1]]  #Not in use
            if i.split('\t')[11] != '.': #Not in use
                for d in i.split('\t')[11].split(','): #Not in use
                    overlapdict[d] = maxovl #Not in use
                    normalcovratio[d] = float(i.split('\t')[12])/splitcov #Not in use
                    classdict[d] = [i.split('\t')[3].split(' ')[0], i.split('\t')[3].split(' ')[1]] #Not in use
        else:
            overlapdict[i.split('\t')[8]] = i.split('\t')[10]
            normalcovratio[i.split('\t')[8]] = float(i.split('\t')[12])/(int(i.split('\t')[10]) + float(i.split('\t')[12]))
            classdict[i.split('\t')[8]] = [i.split('\t')[3].split(' ')[0], i.split('\t')[3].split(' ')[1]]
            if i.split('\t')[11] != '.':
                for d in i.split('\t')[11].split(','):
                    overlapdict[d] = i.split('\t')[10]
                    normalcovratio[d] = float(i.split('\t')[12])/(int(i.split('\t')[10]) + float(i.split('\t')[12]))
                    classdict[d] = [i.split('\t')[3].split(' ')[0], i.split('\t')[3].split(' ')[1]]
    return overlapdict,normalcovratio,classdict

#Filter parse data for >0 long read coverage
def filterparse(parse_data,overlapdict):
    newlist = []
    for i in parse_data:
        try:
            if int(overlapdict[i.split('\t')[8]]) > 0:
                newlist.append(i)
        except:
            pass
    return newlist

#Del and Ins detector
def delins(sv):
    if str(sv) == 'S-Nov_Ins_bp' or sv == 'E-Nov_Ins_bp' or sv == 'Nov_Ins' or sv == 'Del':
        return True

#Collect features for each breakpoint line and rescale the features between 0 and 1
def scalefeature(parse_data,overlapdict,normalcovratio,classdict):
    signdict = {}
    for i in parse_data:
        signdict[i.split('\t')[8]] = []
        for s in i.split('\t')[9].split(','):
            signdict[i.split('\t')[8]].append(s.strip('(%)'))
        signdict[i.split('\t')[8]].append(overlapdict[i.split('\t')[8]])
        signdict[i.split('\t')[8]].append(normalcovratio[i.split('\t')[8]])
        if delins(classdict[i.split('\t')[8]][0]):
            if float(classdict[i.split('\t')[8]][1]) <= 50:
                signdict[i.split('\t')[8]].append(0)
            elif float(classdict[i.split('\t')[8]][1]) <= 100:
                signdict[i.split('\t')[8]].append(0.2)
            elif float(classdict[i.split('\t')[8]][1]) <= 200:
                signdict[i.split('\t')[8]].append(0.4)
            else:
                signdict[i.split('\t')[8]].append(1)
        else:
            signdict[i.split('\t')[8]].append(1)
    l = len(signdict[parse_data[0].split('\t')[8]])
    minmaxlist = []
    for i in range(l):
        tmp = []
        for key in signdict:
            tmp.append(float(signdict[key][i]))
        minmaxlist.append([min(tmp), max(tmp)])
    normsigndict = {}
    for key in signdict:
        normsigndict[key] = []
        for i in range(l):
            normsigndict[key].append(str((float(signdict[key][i]) - minmaxlist[i][0])/(minmaxlist[i][1] - minmaxlist[i][0] + 0.0001)))
    return normsigndict

#Function to average out NN score
def finalprob(leadname, fellow):
    leadprob = probdict[leadname]
    prob = float(leadprob)
    if fellow != '.':
        fellows = fellow.split(',')
        n = len(fellows) + 1
        for i in fellows:
            prob += float(probdict[i])
    elif fellow == '.':
        n = 1
    fprob = float(prob)/n
    return fprob

odata = open(overlap_file, 'r').read().splitlines()
pdata = open(parse_file, 'r').read().splitlines()
odict,normalcovratiodict,classdict = getoverlap(odata)
newpdata = filterparse(pdata,odict)
nsigndict = scalefeature(newpdata,odict,normalcovratiodict,classdict)


tmp = []
tmpkey = []
for key in nsigndict:
    tmp.append(nsigndict[key])
    tmpkey.append(key)

readarray = np.array(tmp, dtype=np.float64)
model = load_model(modelh5)
predictions = model.predict(readarray)
predlist = [float(x[0]) for x in predictions]

probdict = {}
for i in range(len(tmpkey)):
    probdict[tmpkey[i]] = predlist[i]

if mode == 'parse':
    for i in newpdata:
        print i + '\t' + str(np.tanh(0.4*float(odict[i.split('\t')[8]]))*probdict[i.split('\t')[8]]) #tanh((((x-1)/9)+1)) or tanh((10(x-1)/31)+0.8) scale lcov between 0.8 and 0.8 and 10.8, and then tanh transform -0.4507
elif mode == 'overlap':
    for i in odata:
        finalp = finalprob(i.split('\t')[8], i.split('\t')[11])
        print i + '\t' + str(np.tanh(0.4*float(odict[i.split('\t')[8]]))*finalp) #tanh(((x-1)/9)+1) scale lcov between 1 and 3, and then tanh transform, try this tanh(20((x-1)/31)+0.5)
