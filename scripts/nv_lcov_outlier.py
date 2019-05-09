"""
This script calculates the upper limit for long-read SV depth coverage.

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
import numpy as np
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import spline

#Ensure correct number of inputs
if len(argv)<6 or len(argv)>=7:
    sys.exit("Usage: python nv_lcov_outlier.py hsblast.tsv total_genome_size genome.sizes nsplit bedtoolspath > output")

#Assign variable to inputs
hsblast_file = argv[1]
tgsize = int(argv[2])
gsizes = argv[3]
nsplit = int(argv[4])
bedpath = argv[5]

#generate number of genomic points
def ngenerate(gsize):
    if gsize > 10000:
        return 10000
    else:
        return int(gsize*0.9)

#Function to generate bed file from hsblast.tsv
def normalbed(hsblast, nsplits):
    totalbed = []
    if nsplits == 1:
        with open(hsblast) as hs:
            for line in hs:
                coord1 = int(line.split('\t')[1])
                coord2 = int(line.split('\t')[1]) + int(line.split('\t')[2])
                totalbed.append(line.split('\t')[0] + '\t' + str(coord1) + '\t' + str(coord2) + '\t' + line.split('\t')[4])
    else:
        random.seed(3)
        tmp = []
        temp = []
        with open(hsblast) as hs:
            for line in hs:
                tmp.append('\t'.join(line.split('\t')[0:5]))
        tmp.append('null\tnull\tnull\tnull\tnull')
        for i in range(len(tmp)-1):
            if tmp[i].split('\t')[4] == tmp[i+1].split('\t')[4]:
                temp.append(tmp[i])
            else:
                temp.append(tmp[i])
                if random.randint(1,nsplits) == 1:
                    for k in temp:
                        coord1 = int(k.split('\t')[1])
                        coord2 = int(k.split('\t')[1]) + int(k.split('\t')[2])
                        totalbed.append(k.split('\t')[0] + '\t' + str(coord1) + '\t' + str(coord2) + '\t' + k.split('\t')[4])
                temp = []
    return totalbed

#Function to generate bed file from hsblast.tsv
def normalbed_old(hsblast):
    totalbed = []
    with open(hsblast) as hs:
        for line in hs:
            coord1 = int(line.split('\t')[1])
            coord2 = int(line.split('\t')[1]) + int(line.split('\t')[2])
            totalbed.append(line.split('\t')[0] + '\t' + str(coord1) + '\t' + str(coord2) + '\t' + line.split('\t')[4])
    return totalbed

#Median absolute deviation
def mad(x):
    b = 1.4826
    return b*np.mean(abs(x-np.median(x)))

#Plot curve
def curve_old(data, n, nsplit):
    c = range(21)
    p = []
    for i in c:
        p.append(data.count(i))
    p.append(n - sum(p))
    y = [(float(z)/n) for z in p]
    theoretical = [0.0915, 0.0441, 0.1032, 0.1498, 0.1739, 0.1626, 0.1132, 0.0808, 0.0412, 0.0247, 0.0097, 0.0028, 0.0015, 0.0006, 0.0002, 0.0002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    c.append('21')
    xnew = np.linspace(0,21,100)
    smooth = spline(c,y,xnew)
    tsmooth = spline(c,theoretical,xnew)
    params = {'axes.labelsize': 14,'axes.titlesize':17, 'legend.fontsize': 10, 'xtick.labelsize': 12, 'ytick.labelsize': 12, 'font.family': 'Arial, Helvetica, sans-serif'}
    matplotlib.rcParams.update(params)
    fig = plt.figure(figsize=(8, 6))
    fig.patch.set_facecolor('#f6f7f9')
    ax = fig.add_subplot(111)
    ax.plot(xnew, smooth, color='#403f7d', linewidth=2.0)
    ax.plot(xnew, tsmooth, color='#7c7d3f', linewidth=2.0, alpha=1)
    ax.set_facecolor('#ebebff')
    plt.xticks(range(22), ('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','>20'))
    ax.grid(color='w', linestyle='-', linewidth=1)
    vals = ax.get_yticks()
    ax.set_yticklabels(['{:,.1%}'.format(x) for x in vals])
    plt.ylabel('Percentage')
    plt.xlabel('Depth of coverage')
    #plt.title("Read depth across " + str(n) + " points in genome")
    ax.legend(['Input FASTA (spilt_factor=' + str(nsplit) + ')', 'Theoretical 4x hg38 depth'], loc='upper right', ncol=1, fancybox=True)
    plt.savefig('../../nanovar_results/figures/depth_of_coverage.png',bbox_inches='tight', dpi=100, facecolor=fig.get_facecolor(), edgecolor='none')

def curve(data, n, upper_limit):
    c = range(int(upper_limit))
    p = []
    for i in c:
        p.append(data.count(i))
    p.append(n - sum(p))
    y = [(float(z)/n) for z in p]
    #theoretical = [0.0915, 0.0441, 0.1032, 0.1498, 0.1739, 0.1626, 0.1132, 0.0808, 0.0412, 0.0247, 0.0097, 0.0028, 0.0015, 0.0006, 0.0002, 0.0002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    c = range(int(upper_limit)+1) #c.append('21')
    xnew = np.linspace(0,int(upper_limit),100)
    smooth = spline(c,y,xnew)
    #tsmooth = spline(c,theoretical,xnew)
    params = {'axes.labelsize': 14,'axes.titlesize':17, 'legend.fontsize': 10, 'xtick.labelsize': 12, 'ytick.labelsize': 12, 'font.family': 'Arial, Helvetica, sans-serif'}
    matplotlib.rcParams.update(params)
    fig = plt.figure(figsize=(8, 6))
    fig.patch.set_facecolor('#f6f7f9')
    ax = fig.add_subplot(111)
    ax.plot(xnew, smooth, color='#403f7d', linewidth=2.0)
    #ax.plot(xnew, tsmooth, color='#7c7d3f', linewidth=2.0, alpha=1)
    ax.set_facecolor('#ebebff')
    plt.text(int(upper_limit), y[-1], '>=' + str(int(upper_limit)))
    #plt.xticks(np.arange(int(upper_limit),int(upper_limit)+1), ('>' + str(upper_limit)))
    #fig.canvas.draw()
    #labels = [item.get_text() for item in ax.get_xticklabels()]
    #labels[5] = '>' #+ str(int(upper_limit))
    ax.grid(color='w', linestyle='-', linewidth=1)
    vals = ax.get_yticks()
    ax.set_yticklabels(['{:,.1%}'.format(x) for x in vals])
    plt.ylabel('Percentage')
    plt.xlabel('Depth of coverage')
    #plt.title("Read depth across " + str(n) + " points in genome")
    #ax.legend(['Input FASTA (spilt_factor=' + str(nsplit) + ')', 'Theoretical 4x hg38 depth'], loc='upper right', ncol=1, fancybox=True)
    ax.legend(['Input sequencing data'], loc='upper right', ncol=1, fancybox=True)
    plt.savefig('../../nanovar_results/figures/depth_of_coverage.png',bbox_inches='tight', dpi=100, facecolor=fig.get_facecolor(), edgecolor='none')

def main():
    totalbed = normalbed(hsblast_file, nsplit)
    bedout5 = open('bed5', 'w')
    bedout5.write('\n'.join(totalbed))
    bedout5.close()
    os.system(bedpath + ' sort -i bed5 > bed5.sort')
    n = ngenerate(tgsize)
    os.system(bedpath + " random -l 100 -n " + str(n) + " -seed 3 -g " + str(gsizes) + " | " + bedpath + " sort -i - > bed6.sort")
    os.system(bedpath + " intersect -wa -wb -a bed6.sort -b bed5.sort | cut -f 4 | sort | uniq -c | awk -F' ' '{print $1}' > readcovcounts.txt")
    data = open('readcovcounts.txt', 'r').read().splitlines()
    data = [float(i) for i in data]
    zerolist = [0.0]*(n-len(data))
    data2 = data + zerolist
    #curve(data2, n, nsplit)
    med = np.median(data2)
    medad = mad(data2)
    curve(data2, n, round((medad*6) + med, 0))
    #outliers = []
    #for i in data2:
        #if (i-med)/medad > 4: #set empirical cutoff at 4 * mad
            #outliers.append(i)
    #return min(outliers)
    return round((medad*4) + med, 1)

if __name__ == '__main__':
    print main()
