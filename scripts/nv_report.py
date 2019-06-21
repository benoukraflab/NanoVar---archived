"""
This script creates the output HTML report.

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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import seaborn as sns
import os
import math
import datetime
cwd = os.getcwd()

now = datetime.datetime.now()
#This script needs to be run in "nanovar_results" directory
if len(argv) != 8:
    sys.exit("Usage: python nv_report.py sample.total.vcf sample.filtered.vcf num threshold fasta reference_genome fasta_path")

vcf_file = argv[1]
vcf_file2 = argv[2]
num = int(argv[3])
threshold = float(argv[4])
fasta = str(argv[5])
ref = str(argv[6])
fastapath = str(argv[7])

fastaname = os.path.basename(fastapath)

#Setting global figure parameters
params = {'axes.labelsize': 14,'axes.titlesize':17, 'legend.fontsize': 10, 'xtick.labelsize': 12, 'ytick.labelsize': 12, 'font.family': 'Arial, Helvetica, sans-serif'}
matplotlib.rcParams.update(params)

vcf_data = open(vcf_file, 'r').read().splitlines()
vcf = sorted(vcf_data[num:], key=lambda line: float(line.split('\t')[5]), reverse=True)
scorelist = []
ratiolist = []
lcovlist = []

for i in vcf:
    scorelist.append(float(i.split('\t')[5]))
    ratiolist.append(float(i.split('\t')[7].split(';')[6].split('=')[1]))
    lcovlist.append(float(i.split('\t')[7].split(';')[4].split('=')[1]))

fig = plt.figure(figsize=(8, 6))
fig.patch.set_facecolor('#f6f7f9')
ax = fig.add_subplot(111)
ax.scatter(ratiolist,scorelist, c='#7d3f5d', alpha=0.1)
ax.axhline(y=threshold, linewidth=1, color='firebrick')
ax.annotate("Threshold=" + str(threshold), xy=(0.05,threshold+0.2))
ax.set_facecolor('#ebebff')#e6e6ff
plt.ylabel('Confidence score')
plt.xlabel('Breakend read ratio')
plt.ylim(ymin=-0.3)
#plt.title("SV score vs SV read ratio")
plt.savefig('figures/scatter1.png',bbox_inches='tight', dpi=100, facecolor=fig.get_facecolor(), edgecolor='none')

fig = plt.figure(figsize=(8, 6))
fig.patch.set_facecolor('#f6f7f9')
ax = fig.add_subplot(111)
mcov = max(lcovlist)
ax.scatter(lcovlist,scorelist, c='#3f5d7d', alpha=0.1)
ax.axhline(y=threshold, linewidth=1, color='firebrick')
ax.annotate("Threshold=" + str(threshold), xy=(mcov - 2.8,threshold+0.2))
ax.set_facecolor('#ebebff')
plt.ylabel('Confidence score')
plt.xlabel('Number of breakend-supporting reads')
plt.ylim(ymin=-0.3)
#plt.title("SV score vs SV read depth")
plt.savefig('figures/scatter2.png',bbox_inches='tight', dpi=100, facecolor=fig.get_facecolor(), edgecolor='none')

scorelist = []
ratiolist = []
lcovlist = []
svdict = {}
svdict['DEL'] = 0
svdict['INV'] = 0
svdict['DUP'] = 0
svdict['INS'] = 0
svdict['BND'] = 0
dellen = []
inslen = []
invlen = []
duplen = []
bndlen = []
delnolen = 0
insnolen = 0
invnolen = 0
dupnolen = 0
bndnolen = 0
data = []
totalsv = 0
n = 0
for i in vcf:
    if float(i.split('\t')[5]) >= threshold:
        if n <= 999 and float(i.split('\t')[7].split(';')[6].split('=')[1]) <= 1:
            n += 1
            try:
                data.append([str(n), str(i.split('\t')[4]), str(i.split('\t')[0]), str(i.split('\t')[1]), str(i.split('\t')[7].split(';')[1].split('=')[1]), str(i.split('\t')[7].split(';')[2].split('=')[1]), str(i.split('\t')[7].split(';')[-2].split('=')[1]), str(i.split('\t')[5]), str(i.split('\t')[7].split(';')[4].split('=')[1]), str(i.split('\t')[7].split(';')[5].split('=')[1]), str(i.split('\t')[7].split(';')[6].split('=')[1]), str(i.split('\t')[7].split(';')[3].split('=')[1]), str(i.split('\t')[2][:14])])
            except:
                data.append([str(n), str(i.split('\t')[4]), str(i.split('\t')[0]), str(i.split('\t')[1]), str(i.split('\t')[7].split(';')[1].split('=')[1]), str(i.split('\t')[7].split(';')[2].split('=')[1]), '-', str(i.split('\t')[5]), str(i.split('\t')[7].split(';')[4].split('=')[1]), str(i.split('\t')[7].split(';')[5].split('=')[1]), str(i.split('\t')[7].split(';')[6].split('=')[1]), str(i.split('\t')[7].split(';')[3].split('=')[1]), str(i.split('\t')[2][:14])])
        svdict[str(i.split('\t')[4])] += 1
        if str(i.split('\t')[4]) == 'DEL':
            try:
                dellen.append(int(i.split('\t')[7].split(';')[-2].split('=')[1]))
            except:
                delnolen += 1
        elif str(i.split('\t')[4]) == 'INS':
            try:
                inslen.append(int(i.split('\t')[7].split(';')[-2].split('=')[1]))
            except:
                insnolen += 1
        elif str(i.split('\t')[4]) == 'INV':
            try:
                invlen.append(int(i.split('\t')[7].split(';')[-2].split('=')[1]))
            except:
                invnolen += 1
        elif str(i.split('\t')[4]) == 'DUP':
            try:
                duplen.append(int(i.split('\t')[7].split(';')[-2].split('=')[1]))
            except:
                dupnolen += 1
        elif str(i.split('\t')[4]) == 'BND':
            try:
                bndlen.append(int(i.split('\t')[7].split(';')[-2].split('=')[1]))
            except:
                bndnolen += 1
        else:
            sys.exit('SVTYPE error at line ' + str(i))
        totalsv += 1

totalsvlen = []
totalsvlen.append(dellen)
totalsvlen.append(inslen)
totalsvlen.append(invlen)
totalsvlen.append(bndlen)
totalsvlen.append(duplen)
'''
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(ratiolist,scorelist, c='r', alpha=0.1)
plt.savefig('fig3.png',bbox_inches='tight', dpi=100)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(lcovlist,scorelist, c='b', alpha=0.1)
plt.savefig('fig4.png',bbox_inches='tight', dpi=100)
'''
#Boxplot
fig = plt.figure(figsize=(8, 6))
fig.patch.set_facecolor('#f6f7f9')
ax = fig.add_subplot(111)
label = ["Deletion", "Insertion", "Inversion", "Breakend", "TandemDup"]
colors = ['#7d3f5d', '#3f5d7d', '#7d5f3f', '#403f7d', '#3f7c7d']#colors = ['#ff9999', '#9999ff', '#ff99ff', '#ffcc99', '#99ff99']
medianprop = dict(linewidth=1.5,color='black')
bp = ax.boxplot(totalsvlen,whis=[5,95], medianprops=medianprop, showfliers=False, patch_artist=True)
ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

table = []
table.append([len(dellen), len(inslen), len(invlen), len(bndlen), len(duplen)])
table.append([delnolen, insnolen, invnolen, bndnolen, dupnolen])
ax.table(cellText=table, rowLabels=['Known size', 'Unknown size'], colLabels=label, loc='bottom')
ax.set_facecolor('#ebebff')
ax.xaxis.grid(False)
ax.yaxis.grid(color='white', linewidth=1)
plt.ylabel('Length (base pair)')
#plt.title("Distribution of known SV lengths across SV types")
plt.savefig('figures/sv_lengths.png',bbox_inches='tight', dpi=100, facecolor=fig.get_facecolor(), edgecolor='none')

#Donut SV types
fig, ax = plt.subplots(figsize=(10, 5), subplot_kw=dict(aspect="equal"))
fig.patch.set_facecolor('#f6f7f9')
label = [str(svdict['DEL']) + " Deletions",
         str(svdict['INS']) + " Insertions", 
         str(svdict['INV']) + " Inversions", 
         str(svdict['DUP']) + " TandemDups", 
         str(svdict['BND']) + " Breakends"]
x = [svdict['DEL'], svdict['INS'], svdict['INV'], svdict['DUP'], svdict['BND']]

def func(pct, allvals):
    absolute = int(pct/100.*np.sum(allvals))
    return "{:.1f}%".format(pct, absolute)

explode = (0.02,0.02,0.02,0.02, 0.02)
wedges, text, autotexts = ax.pie(x, wedgeprops=dict(width=0.5), startangle=-40, autopct=lambda pct: func(pct,x), textprops=dict(color="white"), counterclock=False, pctdistance=0.75, colors=['#7d3f5d', '#3f5d7d', '#7d5f3f', '#3f7c7d', '#403f7d'], explode = explode) #'#ff9999', '#66b3ff', '#ffb3e6', '#99ff99', '#ffcc99'
bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
kw = dict(xycoords='data', textcoords='data', arrowprops=dict(arrowstyle="-"),
          bbox=bbox_props, zorder=0, va="center")

for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax.annotate(label[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
                 horizontalalignment=horizontalalignment, **kw)

plt.setp(autotexts, size=11)#, weight="bold")
ax.text(0, 0, 'SV types', ha='center')
#plt.title("Distribution of SV types", y=1.10)
#ax.axis('equal')  
#fig.tight_layout()
plt.savefig('figures/sv_type_donut.png',bbox_inches='tight', dpi=100, facecolor=fig.get_facecolor(), edgecolor='none')

#Read length distribution
def measureQlen(fasta):
    qlen = []
    with open(fasta) as f:
        for line in f:
            qlen.append(len(next(f).strip()))
    return qlen

qlen = measureQlen(fasta)
m = math.ceil(max(np.log10(np.array(qlen))))
e = math.floor(min(np.log10(np.array(qlen))))
bins=10**np.linspace(e, m, m*20)
fig = plt.figure(figsize=(8, 6))
fig.patch.set_facecolor('#f6f7f9')
ax = fig.add_subplot(111)
ax.hist(qlen, bins=bins, color="#403f7d", edgecolor='black', linewidth=0.5)
ax.set_xscale('log')
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.ticklabel_format(useOffset=False, style='plain')
ax.set_facecolor('#ebebff')
plt.xlabel("Read length (bases)")
plt.ylabel("Number of reads") 
#plt.title("Distribution of read lengths")
plt.savefig('figures/read_length_dist.png',bbox_inches='tight', dpi=100, facecolor=fig.get_facecolor(), edgecolor='none')

'''
fig = plt.figure()
ax = fig.add_subplot(111)
label = ["Deletion", "Insertion", "Inversion", "Breakend", "TandemDup"]
table = []
table.append([svdict['DEL']])
table.append([svdict['INS']])
table.append([svdict['INV']])
table.append([svdict['BND']])
table.append([svdict['DUP']])
t = ax.table(cellText=table, rowLabels=label, colLabels=['Number of SVs'], loc='center')
ax.axis('tight')
ax.axis('off')
t.scale(0.3, 1)
plt.savefig('fig5.png',bbox_inches='tight', dpi=100)
'''

#Create html
vcf_name = open(vcf_file, 'r').name
html = open(str(vcf_name[:-9]) + 'report.html', 'w')
begin = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <meta http-equiv="x-ua-compatible" content="ie=edge">
    <title>NanoVar Report</title>
    <!-- Font Awesome -->
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css">
    <!-- Bootstrap core CSS -->
    <link href="css/bootstrap.min.css" rel="stylesheet">
    <!-- Material Design Bootstrap -->
    <link href="css/mdb.min.css" rel="stylesheet">
    <!-- DataTable CSS -->
    <link href="css/jquery.dataTables.min.css" rel="stylesheet">
    <!-- DataTable buttons CSS  -->
    <link href="css/buttons.dataTables.min.css" rel="stylesheet">
    <style>
        body {
            background-color: #f6f7f9;
        }
        h1 { 
            margin-left: 4px;
            font-weight: bold;
            font-size: 40px;
            font-family: Arial, Helvetica, sans-serif;
        }
        h2, h3, h4, h5, p { 
            margin-left: 4px;
            font-family: Arial, Helvetica, sans-serif;
        }
        table {
            font-family: Arial, Helvetica, sans-serif;
        }
        #image {
            text-align:center;
        }
        #margin {
            margin-left: 35px;
            margin-right: 35px;
            margin-top: 25px;
        }
        .floatleft {
            float: left;
            width: 45%;
            margin-right: 5%;
            padding-top: 80px
        }
        .floatright {
            float: right;
            width: 45%;
            margin-left: 5%;
        }
        .container {
            overflow;
            max-width: 1400px;
        }
        .aligncenter {
            float: center;
        }
        .alignleft {
            float: left;
        }
        .alignright {
            float: right;
        }
    </style>
</head>
<body>
<div id="margin">
    <div>
        <h1 style="text-align:center;">NanoVar Report</h1>
        <h5 style="text-align:center;">Version: NanoVar 1.0.1</h5>
        <h5 style="text-align:center;">""" + str(now.strftime("%a %d %B %Y")) + """</h5>
        <h5 style="text-align:center;">""" + fastaname + """</h5>
        <br>
        <br>
        <div style="background-color:#3f5d7d;color:#f6f7f9">
            <h4 style="text-align:center;">Run details</h4>
        </div>
    </div>
    <div style="clear: both;"></div>
    <br>
    <table width="50%" class="table-sm table-bordered" style="table-layout:fixed" align="center">
        <col width="150">
        <thead>
        <tr>
            <th style="background-color:#2e6096; color:#f6f7f9; font-weight: bold; text-align:center;">Variable</th>
            <th style="background-color:#2e6096; color:#f6f7f9; font-weight: bold; text-align:center;">Value</th>
        </tr>
        </thead>
        <tbody>
        <tr>
            <td>Input FASTA:</td>
            <td style="word-wrap:break-word;">""" + fastapath + """</td>
        </tr>
        <tr>
            <td>Reference genome:</td>
            <td style="word-wrap:break-word;">""" + ref + """</td>
        </tr>
        <tr>
            <td>Output VCF:</td>
            <td style="word-wrap:break-word;">""" + cwd + '/' + vcf_file2 + """</td>
        </tr>
        <tr>
            <td>SV score threshold:</td>
            <td>""" + str(threshold) +  """</td>
        </tr>
        </tbody>
    </table>
    <br>
    <br>
    <br>
    <br>
    <div style="background-color:#366b6c;color:#f6f7f9;">
        <h3 style="text-align:center;">Results</h3>
    </div>
    <p style="clear: both;">
    <br>
    <h4 style="text-align:center;"><u>1. Table of output SVs</u></h4>
    <h5 style="text-align:center;">Showing """ + str(len(data)) + """ out of """ + str(totalsv) + """ total SVs</h5>
    <table id="NanoVar_report_table" class="table table-striped table-bordered table-sm" cellspacing="0" width="100%">
        <thead>
            <tr>
                <th class="th-sm">#
                    <i class="fa fa-sort float-right" aria-hidden="true"></i>
                </th>
                <th class="th-sm">SV type
                    <i class="fa fa-sort float-right" aria-hidden="true"></i>
                </th>
                <th class="th-sm">Chrom1
                    <i class="fa fa-sort float-right" aria-hidden="true"></i>
                </th>
                <th class="th-sm">Pos1
                    <i class="fa fa-sort float-right" aria-hidden="true"></i>
                </th>
                <th class="th-sm">Chrom2
                    <i class="fa fa-sort float-right" aria-hidden="true"></i>
                </th>
                <th class="th-sm">Pos2
                    <i class="fa fa-sort float-right" aria-hidden="true"></i>
                </th>
                <th class="th-sm">Length
                    <i class="fa fa-sort float-right" aria-hidden="true"></i>
                </th>
                <th class="th-sm">Score
                    <i class="fa fa-sort float-right" aria-hidden="true"></i>
                </th>
                <th class="th-sm">No. of breakend-supporting reads
                    <i class="fa fa-sort float-right" aria-hidden="true"></i>
                </th>
                <th class="th-sm">No. of breakend-opposing reads
                    <i class="fa fa-sort float-right" aria-hidden="true"></i>
                </th>
                <th title="No. of breakend-supporting reads/total reads at breakend" class="th-sm">Breakend read ratio
                    <i class="fa fa-sort float-right" aria-hidden="true"></i>
                </th>
                <th class="th-sm">No. of breakend-supporting short reads
                    <i class="fa fa-sort float-right" aria-hidden="true"></i>
                </th>
                <th title="ID_code ~ readname" class="th-sm">SV ID
                    <i class="fa fa-sort float-right" aria-hidden="true"></i>
                </th>
            </tr>
        </thead>
        <tbody>
"""
html.write(begin)

#Create table body
row = ""
for i in data:
    row += """          <tr>
"""
    for j in i:
        row += "                <td>" + j + """</td>
"""
    row += """          </tr>
"""
row += """      </tbody>
    </table>
    <br>
    <br>
    <div id="image">
        <figure>
            <h4 style="text-align:center;"><u>2. Distribution of SV types</u></h4>
            <img src="figures/sv_type_donut.png" alt="2. Distribution of SV types" title=""" + '"' + cwd + '/figures/sv_type_donut.png' + '"' + """>
        </figure>
        <br>
        <br>
        <figure>
            <h4 style="text-align:center;"><u>3. Size distribution of SVs across SV types</u></h4>
            <img src="figures/sv_lengths.png" alt="3. Size distribution of SVs across SV types" title=""" + '"' + cwd + '/figures/sv_lengths.png' + '"' + """>
        </figure>
        <br>
        <br>
        <figure>
            <h4 style="text-align:center;"><u>4. Scatter plot between SV confidence score and read ratio</u></h4>
            <img src="figures/scatter1.png" alt="4. Scatter plot between SV confidence score and read ratio" title=""" + '"' + cwd + '/figures/scatter1.png' + '"' + """>
        </figure>
        <br>
        <br>
        <figure>
            <h4 style="text-align:center;"><u>5. Scatter plot between SV confidence score and read depth</u></h4>
            <img src="figures/scatter2.png" alt="5. Scatter plot between SV confidence score and read depth" title=""" + '"' + cwd + '/figures/scatter2.png' + '"' + """>
        </figure>
        <br>
        <br>
        <div style="background-color:#403f7d;color:#f6f7f9;">
            <h3 style="text-align:center;">QC Statistics</h3>
        </div>
        <br>
        <figure>
            <h4 style="text-align:center;"><u>6. Distribution of read lengths of input FASTA</u></h4>
            <img src="figures/read_length_dist.png" alt="6. Distribution of read lengths of input FASTA" title=""" + '"' + cwd + '/figures/read_length_dist.png' + '"' + """>
        </figure>
        <br>
        <br>
        <figure>
            <h4 style="text-align:center;"><u>7. Depth of coverage across genome</u></h4>
            <img src="figures/depth_of_coverage.png" alt="7. Depth of coverage across genome" title=""" + '"' + cwd + '/figures/depth_of_coverage.png' + '"' + """>
        </figure>
        <br>
    </div>
    <br>
    <br>
    <!-- SCRIPTS -->
    <!-- JQuery -->
    <script type="text/javascript" src="js/jquery-3.3.1.min.js"></script>
    <!-- datatable javascript  -->
    <script type="text/javascript" src="js/jquery.dataTables.min.js"></script>
    <!-- datatable buttons javascript  -->
    <script type="text/javascript" src="js/dataTables.buttons.min.js"></script>
    <!-- buttons flash javascript  -->
    <script type="text/javascript" src="js/buttons.flash.min.js"></script>
    <!-- buttons jszip javascript  -->
    <script type="text/javascript" src="js/jszip.min.js"></script>
    <!-- buttons html5 buttons javascript  -->
    <script type="text/javascript" src="js/buttons.html5.min.js"></script>
    <!-- Bootstrap tooltips -->
    <script type="text/javascript" src="js/popper.min.js"></script>
    <!-- Bootstrap core JavaScript -->
    <script type="text/javascript" src="js/bootstrap.min.js"></script>
    <!-- MDB core JavaScript -->
    <script type="text/javascript" src="js/mdb.min.js"></script>
    <script type="text/javascript">
        $('#NanoVar_report_table').DataTable({
            "scrollX": true,
            "scrollY": 200,
            dom: 'Bfrtip',
            buttons: [
                'copy', 'csv', 'excel'
            ]
        });
        $('.dataTables_length').addClass('bs-select');
    </script>
</div>
</body>
</html>
"""
html.write(row)
html.close()


'''
#<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css">
#scroll {
        position: relative;
        width: 1000px;
        height: 150px;
        overflow: auto;
        }
        
        #scroll table {
        width: 100%;
        }
        
        #scroll table thead th {
        position: absolute;   
        border:1px solid red;
        }
'''