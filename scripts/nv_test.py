"""
This script is used for make check.

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

if len(argv)<5 or len(argv)>=6:
    sys.exit("Usage: python nv_test.py md5sumfile [b1 or b0] maindir [linux or mac]")

md5_file = argv[1]
bowtie = argv[2]
maindir = argv[3]
machtype = argv[4]
md5_names = open(md5_file, 'r').read().splitlines()
md5 = []
for i in md5_names:
    md5.append(i.split(' ')[0])

def checker(test, ref):
    if test == ref:
        return 1
    else:
        return 0

def main():
    if machtype == "linux":
        if bowtie == "b1": #Enabled bowtie2
            if checker(md5[0], "f846237582adf9340d40f9c23e2563e4"):
                print "<Reference> ---------: Checked"
            else:
                print "<Reference> ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[1], "353914b702cdb564ebfa41fa8819bc77"):
                print "<longfa>    ---------: Checked"
            else:
                print "<longfa>    ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[2], "3168ccfb91a64a6d0e3174b05bab9adf"):
                print "<shortfa1>  ---------: Checked"
            else:
                print "<shortfa1>  ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[3], "399fa36013ebec490df7923963ed6892"):
                print "<shortfa2>  ---------: Checked"
            else:
                print "<shortfa2>  ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[4], "c76d425ecd5af4ef0c6059da6e4bb847"):
                print "<fai>       ---------: Passed"
            else:
                print "<fai>       ---------: Failed"
                sys.exit("Samtools faidx failed, please ensure samtools is functional")
            if checker(md5[5], "b95f9a36c118811952a149562a7f2d61"):
                print "<nhr>       ---------: Passed"
            else:
                print "<nhr>       ---------: Failed"
                sys.exit("makeblastdb failed, please download Blast v2.7.1+ source, build it and copy 'makeblastdb' binaries into " + maindir + "/blast/")
            if checker(md5[7], "45d73fff1fd267ff5606831f99a0952d"):
                print "<nsq>       ---------: Passed"
            else:
                print "<nsq>       ---------: Failed"
                sys.exit("makeblastdb failed, please download Blast v2.7.1+ source, build it and copy 'makeblastdb' binaries into " + maindir + "/blast/")
            if checker(md5[8], "dd20f9ca6221d6ccb3ac9f686e7564b1"):
                print "<bwt>       ---------: Passed"
            else:
                print "<bwt>       ---------: Failed"
                sys.exit("hsblast index failed, please re-configure and re-compile NanoVar")
            if checker(md5[10], "3c4bddf8ddc525bc0fee6ea2ad1ac714"):
                print "<sa>        ---------: Passed"
            else:
                print "<sa>        ---------: Failed"
                sys.exit("hsblast index failed, please re-configure and re-compile NanoVar")
            if checker(md5[11], "8506e7d6ecdb13a205dfb5ddc0ee6870"):
                print "<sequence>  ---------: Passed"
            else:
                print "<sequence>  ---------: Failed"
                sys.exit("hsblast index failed, please re-configure and re-compile NanoVar")
            if checker(md5[12], "2c5c01b1dcfe571a074120d5ae151051"):
                print "<counts>    ---------: Passed"
            else:
                print "<counts>    ---------: Failed"
                sys.exit("windowmasker failed, please download Blast v2.7.1+ source, build it and copy 'windowmasker' binaries into " + maindir + "/blast/")
            if checker(md5[13], "0b30579e91aeb001cdb57701cbe31579"):
                print "<obinary>   ---------: Passed"
            else:
                print "<obinary>   ---------: Failed"
                sys.exit("windowmasker failed, please download Blast v2.7.1+ source, build it and copy 'windowmasker' binaries into " + maindir + "/blast/")
            if checker(md5[14], "8c68ba585e3584c0c8e8b2ad13964faf"):
                print "<align>     ---------: Passed"
            else:
                print "<align>     ---------: Failed"
                sys.exit("hsblast align failed, please re-configure and re-compile NanoVar")
            if checker(md5[15], "c070563dbae248542179c5080146af5d"):
                print "<parse>     ---------: Passed"
            else:
                print "<parse>     ---------: Failed"
                sys.exit("bedtools sort failed, please ensure bedtools version >= 2.26.0 and functional")
            if checker(md5[16], "0227aa5acaec17ddabf2aa1f454c48f3"):
                print "<overlap>   ---------: Passed"
            else:
                print "<overlap>   ---------: Failed"
                sys.exit("bedtools intersect or awk failed, please ensure bedtools version >= 2.26.0 and awk are functional")
            if checker(md5[17], "632abc262071733089f262752f374852"):
                print "<ANN>       ---------: Passed"
            else:
                print "<ANN>       ---------: Failed"
                sys.exit("python2.7 failed, please ensure python2.7 is in PATH and functional")
            if checker(md5[18], "2ce838d744897896ccc3294f86f9ad8e"):
                print "<ANN0>      ---------: Passed"
            else:
                print "<ANN0>      ---------: Failed"
                sys.exit("awk failed, please ensure awk is functional")
            if checker(md5[19], "106e03e6661155dd2c650f62156da114"):
                print "<bt1>       ---------: Passed"
            else:
                print "<bt1>       ---------: Failed"
                sys.exit("bowtie2-build failed, please ensure bowtie2-build is functional")
            if checker(md5[20], "eae1a5d74420fb4373fa9cd1c22fa330"):
                print "<bt2>       ---------: Passed"
            else:
                print "<bt2>       ---------: Failed"
                sys.exit("bowtie2-build failed, please ensure bowtie2-build is functional")
            if checker(md5[21], "28226cd64614e0e25cbd3bc9fa02602e"):
                print "<bt3>       ---------: Passed"
            else:
                print "<bt3>       ---------: Failed"
                sys.exit("bowtie2-build failed, please ensure bowtie2-build is functional")
            if checker(md5[22], "c1bd49fe25edf2cf908d8b7a3fa1fa4b"):
                print "<bt4>       ---------: Passed"
            else:
                print "<bt4>       ---------: Failed"
                sys.exit("bowtie2-build failed, please ensure bowtie2-build is functional")
            if checker(md5[23], "f1d1e3d4306eae03933109cb3bb7db6f"):
                print "<revbt1>    ---------: Passed"
            else:
                print "<revbt1>    ---------: Failed"
                sys.exit("bowtie2-build failed, please ensure bowtie2-build is functional")
            if checker(md5[24], "91da4bff0e657e42c000b7bd412c347b"):
                print "<revbt2>    ---------: Passed"
            else:
                print "<revbt2>    ---------: Failed"
                sys.exit("bowtie2-build failed, please ensure bowtie2-build is functional")
            if checker(md5[25], "1757a67825471a4a86a37fd1b13fbbf3"):
                print "<bam-sam>   ---------: Passed"
            else:
                print "<bam-sam>   ---------: Failed"
                sys.exit("bowtie2 align or samtools view/sort failed, please ensure bowtie2 and samtools are functional")
            if checker(md5[26], "ca5191a9bf508edd401ad1fcd103c229"):
                print "<fq1>       ---------: Passed"
            else:
                print "<fq1>       ---------: Failed"
                sys.exit("samtools fastq failed, please ensure samtools is functional")
            if checker(md5[27], "296303b787a8d9c08989fec3c32073eb"):
                print "<fq2>       ---------: Passed"
            else:
                print "<fq2>       ---------: Failed"
                sys.exit("samtools fastq failed, please ensure samtools is functional")
            if checker(md5[28], "e49d43c4147ec01e996594a0fc7c863c"):
                print "<cov>       ---------: Passed"
            else:
                print "<cov>       ---------: Failed"
                sys.exit("bedtools map or awk or sort failed, please ensure bedtools and awk and sort are functional")
            if checker(md5[29], "24438027e0f8cf21643ab1019d6b96a0"):
                print "<t-vcf>     ---------: Passed"
            else:
                print "<t-vcf>     ---------: Failed"
                sys.exit("bedtools sort failed, please ensure bedtools is functional")
            if checker(md5[30], "207b8eaf23dca13920bc4b2e04fee283"):
                print "<filt-vcf>  ---------: Passed"
            else:
                print "<filt-vcf>  ---------: Failed"
                sys.exit("bedtools sort or awk failed, please ensure bedtools and awk are functional")
            if checker(md5[31], "05d1a814056f2795113d3858ea2e094c"):
                print "<svoverlap> ---------: Passed"
            else:
                print "<svoverlap> ---------: Failed"
                sys.exit("awk or sed failed, please ensure awk and sed are functional")
            os.system('echo "AllPassed" > ' + maindir + '/test/AllPassed')
        elif bowtie == "b0": #Disabled bowtie2
            if checker(md5[0], "f846237582adf9340d40f9c23e2563e4"):
                print "<Reference> ---------: Checked"
            else:
                print "<Reference> ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[1], "353914b702cdb564ebfa41fa8819bc77"):
                print "<longfa>    ---------: Checked"
            else:
                print "<longfa>    ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[2], "3168ccfb91a64a6d0e3174b05bab9adf"):
                print "<shortfa1>  ---------: Checked"
            else:
                print "<shortfa1>  ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[3], "399fa36013ebec490df7923963ed6892"):
                print "<shortfa2>  ---------: Checked"
            else:
                print "<shortfa2>  ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[4], "c76d425ecd5af4ef0c6059da6e4bb847"):
                print "<fai>       ---------: Passed"
            else:
                print "<fai>       ---------: Failed"
                sys.exit("Samtools faidx failed, please ensure samtools is functional")
            if checker(md5[5], "b95f9a36c118811952a149562a7f2d61"):
                print "<nhr>       ---------: Passed"
            else:
                print "<nhr>       ---------: Failed"
                sys.exit("makeblastdb failed, please download Blast v2.7.1+ source, build it and copy 'makeblastdb' binaries into " + maindir + "/blast/")
            if checker(md5[7], "45d73fff1fd267ff5606831f99a0952d"):
                print "<nsq>       ---------: Passed"
            else:
                print "<nsq>       ---------: Failed"
                sys.exit("makeblastdb failed, please download Blast v2.7.1+ source, build it and copy 'makeblastdb' binaries into " + maindir + "/blast/")
            if checker(md5[8], "dd20f9ca6221d6ccb3ac9f686e7564b1"):
                print "<bwt>       ---------: Passed"
            else:
                print "<bwt>       ---------: Failed"
                sys.exit("hsblast index failed, please re-configure and re-compile NanoVar")
            if checker(md5[10], "3c4bddf8ddc525bc0fee6ea2ad1ac714"):
                print "<sa>        ---------: Passed"
            else:
                print "<sa>        ---------: Failed"
                sys.exit("hsblast index failed, please re-configure and re-compile NanoVar")
            if checker(md5[11], "8506e7d6ecdb13a205dfb5ddc0ee6870"):
                print "<sequence>  ---------: Passed"
            else:
                print "<sequence>  ---------: Failed"
                sys.exit("hsblast index failed, please re-configure and re-compile NanoVar")
            if checker(md5[12], "2c5c01b1dcfe571a074120d5ae151051"):
                print "<counts>    ---------: Passed"
            else:
                print "<counts>    ---------: Failed"
                sys.exit("windowmasker failed, please download Blast v2.7.1+ source, build it and copy 'windowmasker' binaries into " + maindir + "/blast/")
            if checker(md5[13], "0b30579e91aeb001cdb57701cbe31579"):
                print "<obinary>   ---------: Passed"
            else:
                print "<obinary>   ---------: Failed"
                sys.exit("windowmasker failed, please download Blast v2.7.1+ source, build it and copy 'windowmasker' binaries into " + maindir + "/blast/")
            if checker(md5[14], "8c68ba585e3584c0c8e8b2ad13964faf"):
                print "<align>     ---------: Passed"
            else:
                print "<align>     ---------: Failed"
                sys.exit("hsblast align failed, please re-configure and re-compile NanoVar")
            if checker(md5[15], "c070563dbae248542179c5080146af5d"):
                print "<parse>     ---------: Passed"
            else:
                print "<parse>     ---------: Failed"
                sys.exit("bedtools sort failed, please ensure bedtools version >= 2.26.0 and functional")
            if checker(md5[16], "0227aa5acaec17ddabf2aa1f454c48f3"):
                print "<overlap>   ---------: Passed"
            else:
                print "<overlap>   ---------: Failed"
                sys.exit("bedtools intersect or awk failed, please ensure bedtools version >= 2.26.0 and awk are functional")
            if checker(md5[17], "632abc262071733089f262752f374852"):
                print "<ANN>       ---------: Passed"
            else:
                print "<ANN>       ---------: Failed"
                sys.exit("python2.7 failed, please ensure python2.7 is in PATH and functional")
            if checker(md5[18], "2ce838d744897896ccc3294f86f9ad8e"):
                print "<ANN0>      ---------: Passed"
            else:
                print "<ANN0>      ---------: Failed"
                sys.exit("awk failed, please ensure awk is functional")
            if checker(md5[19], "8d4196ea6d2553ebc4b7b2c7cf9b4bb3"):
                print "<cov>       ---------: Passed"
            else:
                print "<cov>       ---------: Failed"
                sys.exit("awk or sort failed, please ensure awk and sort are functional")
            if checker(md5[20], "ef6fdc3833dd26ca460b652fb42f3481"):
                print "<t-vcf>     ---------: Passed"
            else:
                print "<t-vcf>     ---------: Failed"
                sys.exit("bedtools sort failed, please ensure bedtools is functional")
            if checker(md5[21], "8766c25990aee00fe07e7b57aa5b1daa"):
                print "<filt-vcf>  ---------: Passed"
            else:
                print "<filt-vcf>  ---------: Failed"
                sys.exit("bedtools sort or awk failed, please ensure bedtools and awk are functional")
            if checker(md5[22], "05d1a814056f2795113d3858ea2e094c"):
                print "<svoverlap> ---------: Passed"
            else:
                print "<svoverlap> ---------: Failed"
                sys.exit("awk or sed failed, please ensure awk and sed are functional")
            os.system('echo "AllPassed" > ' + maindir + '/test/AllPassed')
    elif machtype == "mac":
        if bowtie == "b1": #Enabled bowtie2
            if checker(md5[0], "f846237582adf9340d40f9c23e2563e4"):
                print "<Reference> ---------: Checked"
            else:
                print "<Reference> ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[1], "353914b702cdb564ebfa41fa8819bc77"):
                print "<longfa>    ---------: Checked"
            else:
                print "<longfa>    ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[2], "3168ccfb91a64a6d0e3174b05bab9adf"):
                print "<shortfa1>  ---------: Checked"
            else:
                print "<shortfa1>  ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[3], "399fa36013ebec490df7923963ed6892"):
                print "<shortfa2>  ---------: Checked"
            else:
                print "<shortfa2>  ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[4], "c76d425ecd5af4ef0c6059da6e4bb847"):
                print "<fai>       ---------: Passed"
            else:
                print "<fai>       ---------: Failed"
                sys.exit("Samtools faidx failed, please ensure samtools is functional")
            if checker(md5[5], "b95f9a36c118811952a149562a7f2d61"):
                print "<nhr>       ---------: Passed"
            else:
                print "<nhr>       ---------: Failed"
                sys.exit("makeblastdb failed, please download Blast v2.7.1+ source, build it and copy 'makeblastdb' binaries into " + maindir + "/blast/")
            if checker(md5[7], "45d73fff1fd267ff5606831f99a0952d"):
                print "<nsq>       ---------: Passed"
            else:
                print "<nsq>       ---------: Failed"
                sys.exit("makeblastdb failed, please download Blast v2.7.1+ source, build it and copy 'makeblastdb' binaries into " + maindir + "/blast/")
            if checker(md5[8], "dd20f9ca6221d6ccb3ac9f686e7564b1"):
                print "<bwt>       ---------: Passed"
            else:
                print "<bwt>       ---------: Failed"
                sys.exit("hsblast index failed, please re-configure and re-compile NanoVar")
            if checker(md5[10], "3c4bddf8ddc525bc0fee6ea2ad1ac714"):
                print "<sa>        ---------: Passed"
            else:
                print "<sa>        ---------: Failed"
                sys.exit("hsblast index failed, please re-configure and re-compile NanoVar")
            if checker(md5[11], "8506e7d6ecdb13a205dfb5ddc0ee6870"):
                print "<sequence>  ---------: Passed"
            else:
                print "<sequence>  ---------: Failed"
                sys.exit("hsblast index failed, please re-configure and re-compile NanoVar")
            if checker(md5[12], "2c5c01b1dcfe571a074120d5ae151051"):
                print "<counts>    ---------: Passed"
            else:
                print "<counts>    ---------: Failed"
                sys.exit("windowmasker failed, please download Blast v2.7.1+ source, build it and copy 'windowmasker' binaries into " + maindir + "/blast/")
            if checker(md5[13], "0b30579e91aeb001cdb57701cbe31579"):
                print "<obinary>   ---------: Passed"
            else:
                print "<obinary>   ---------: Failed"
                sys.exit("windowmasker failed, please download Blast v2.7.1+ source, build it and copy 'windowmasker' binaries into " + maindir + "/blast/")
            if checker(md5[14], ""):
                print "<align>     ---------: Passed"
            else:
                print "<align>     ---------: Failed"
                sys.exit("hsblast align failed, please re-configure and re-compile NanoVar")
            if checker(md5[15], ""):
                print "<parse>     ---------: Passed"
            else:
                print "<parse>     ---------: Failed"
                sys.exit("bedtools sort failed, please ensure bedtools version >= 2.26.0 and functional")
            if checker(md5[16], ""):
                print "<overlap>   ---------: Passed"
            else:
                print "<overlap>   ---------: Failed"
                sys.exit("bedtools intersect or awk failed, please ensure bedtools version >= 2.26.0 and awk are functional")
            if checker(md5[17], ""):
                print "<ANN>       ---------: Passed"
            else:
                print "<ANN>       ---------: Failed"
                sys.exit("python2.7 failed, please ensure python2.7 is in PATH and functional")
            if checker(md5[18], ""):
                print "<ANN0>      ---------: Passed"
            else:
                print "<ANN0>      ---------: Failed"
                sys.exit("awk failed, please ensure awk is functional")
            if checker(md5[19], ""):
                print "<bt1>       ---------: Passed"
            else:
                print "<bt1>       ---------: Failed"
                sys.exit("bowtie2-build failed, please ensure bowtie2-build is functional")
            if checker(md5[20], ""):
                print "<bt2>       ---------: Passed"
            else:
                print "<bt2>       ---------: Failed"
                sys.exit("bowtie2-build failed, please ensure bowtie2-build is functional")
            if checker(md5[21], ""):
                print "<bt3>       ---------: Passed"
            else:
                print "<bt3>       ---------: Failed"
                sys.exit("bowtie2-build failed, please ensure bowtie2-build is functional")
            if checker(md5[22], ""):
                print "<bt4>       ---------: Passed"
            else:
                print "<bt4>       ---------: Failed"
                sys.exit("bowtie2-build failed, please ensure bowtie2-build is functional")
            if checker(md5[23], ""):
                print "<revbt1>    ---------: Passed"
            else:
                print "<revbt1>    ---------: Failed"
                sys.exit("bowtie2-build failed, please ensure bowtie2-build is functional")
            if checker(md5[24], ""):
                print "<revbt2>    ---------: Passed"
            else:
                print "<revbt2>    ---------: Failed"
                sys.exit("bowtie2-build failed, please ensure bowtie2-build is functional")
            if checker(md5[25], ""):
                print "<bam-sam>   ---------: Passed"
            else:
                print "<bam-sam>   ---------: Failed"
                sys.exit("bowtie2 align or samtools view/sort failed, please ensure bowtie2 and samtools are functional")
            if checker(md5[26], ""):
                print "<fq1>       ---------: Passed"
            else:
                print "<fq1>       ---------: Failed"
                sys.exit("samtools fastq failed, please ensure samtools is functional")
            if checker(md5[27], ""):
                print "<fq2>       ---------: Passed"
            else:
                print "<fq2>       ---------: Failed"
                sys.exit("samtools fastq failed, please ensure samtools is functional")
            if checker(md5[28], ""):
                print "<cov>       ---------: Passed"
            else:
                print "<cov>       ---------: Failed"
                sys.exit("bedtools map or awk or sort failed, please ensure bedtools and awk and sort are functional")
            if checker(md5[29], ""):
                print "<t-vcf>     ---------: Passed"
            else:
                print "<t-vcf>     ---------: Failed"
                sys.exit("bedtools sort failed, please ensure bedtools is functional")
            if checker(md5[30], ""):
                print "<filt-vcf>  ---------: Passed"
            else:
                print "<filt-vcf>  ---------: Failed"
                sys.exit("bedtools sort or awk failed, please ensure bedtools and awk are functional")
            if checker(md5[31], ""):
                print "<svoverlap> ---------: Passed"
            else:
                print "<svoverlap> ---------: Failed"
                sys.exit("awk or sed failed, please ensure awk and sed are functional")
            os.system('echo "AllPassed" > ' + maindir + '/test/AllPassed')
        elif bowtie == "b0": #Disabled bowtie2
            if checker(md5[0], "f846237582adf9340d40f9c23e2563e4"):
                print "<Reference> ---------: Checked"
            else:
                print "<Reference> ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[1], "353914b702cdb564ebfa41fa8819bc77"):
                print "<longfa>    ---------: Checked"
            else:
                print "<longfa>    ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[2], "3168ccfb91a64a6d0e3174b05bab9adf"):
                print "<shortfa1>  ---------: Checked"
            else:
                print "<shortfa1>  ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[3], "399fa36013ebec490df7923963ed6892"):
                print "<shortfa2>  ---------: Checked"
            else:
                print "<shortfa2>  ---------: Failed"
                sys.exit("Test data file corrupted, please re-download NanoVar")
            if checker(md5[4], "c76d425ecd5af4ef0c6059da6e4bb847"):
                print "<fai>       ---------: Passed"
            else:
                print "<fai>       ---------: Failed"
                sys.exit("Samtools faidx failed, please ensure samtools is functional")
            if checker(md5[5], "b95f9a36c118811952a149562a7f2d61"):
                print "<nhr>       ---------: Passed"
            else:
                print "<nhr>       ---------: Failed"
                sys.exit("makeblastdb failed, please download Blast v2.7.1+ source, build it and copy 'makeblastdb' binaries into " + maindir + "/blast/")
            if checker(md5[7], "45d73fff1fd267ff5606831f99a0952d"):
                print "<nsq>       ---------: Passed"
            else:
                print "<nsq>       ---------: Failed"
                sys.exit("makeblastdb failed, please download Blast v2.7.1+ source, build it and copy 'makeblastdb' binaries into " + maindir + "/blast/")
            if checker(md5[8], "dd20f9ca6221d6ccb3ac9f686e7564b1"):
                print "<bwt>       ---------: Passed"
            else:
                print "<bwt>       ---------: Failed"
                sys.exit("hsblast index failed, please re-configure and re-compile NanoVar")
            if checker(md5[10], "3c4bddf8ddc525bc0fee6ea2ad1ac714"):
                print "<sa>        ---------: Passed"
            else:
                print "<sa>        ---------: Failed"
                sys.exit("hsblast index failed, please re-configure and re-compile NanoVar")
            if checker(md5[11], "8506e7d6ecdb13a205dfb5ddc0ee6870"):
                print "<sequence>  ---------: Passed"
            else:
                print "<sequence>  ---------: Failed"
                sys.exit("hsblast index failed, please re-configure and re-compile NanoVar")
            if checker(md5[12], "2c5c01b1dcfe571a074120d5ae151051"):
                print "<counts>    ---------: Passed"
            else:
                print "<counts>    ---------: Failed"
                sys.exit("windowmasker failed, please download Blast v2.7.1+ source, build it and copy 'windowmasker' binaries into " + maindir + "/blast/")
            if checker(md5[13], "0b30579e91aeb001cdb57701cbe31579"):
                print "<obinary>   ---------: Passed"
            else:
                print "<obinary>   ---------: Failed"
                sys.exit("windowmasker failed, please download Blast v2.7.1+ source, build it and copy 'windowmasker' binaries into " + maindir + "/blast/")
            if checker(md5[14], ""):
                print "<align>     ---------: Passed"
            else:
                print "<align>     ---------: Failed"
                sys.exit("hsblast align failed, please re-configure and re-compile NanoVar")
            if checker(md5[15], ""):
                print "<parse>     ---------: Passed"
            else:
                print "<parse>     ---------: Failed"
                sys.exit("bedtools sort failed, please ensure bedtools version >= 2.26.0 and functional")
            if checker(md5[16], ""):
                print "<overlap>   ---------: Passed"
            else:
                print "<overlap>   ---------: Failed"
                sys.exit("bedtools intersect or awk failed, please ensure bedtools version >= 2.26.0 and awk are functional")
            if checker(md5[17], ""):
                print "<ANN>       ---------: Passed"
            else:
                print "<ANN>       ---------: Failed"
                sys.exit("python2.7 failed, please ensure python2.7 is in PATH and functional")
            if checker(md5[18], ""):
                print "<ANN0>      ---------: Passed"
            else:
                print "<ANN0>      ---------: Failed"
                sys.exit("awk failed, please ensure awk is functional")
            if checker(md5[19], ""):
                print "<cov>       ---------: Passed"
            else:
                print "<cov>       ---------: Failed"
                sys.exit("awk or sort failed, please ensure awk and sort are functional")
            if checker(md5[20], ""):
                print "<t-vcf>     ---------: Passed"
            else:
                print "<t-vcf>     ---------: Failed"
                sys.exit("bedtools sort failed, please ensure bedtools is functional")
            if checker(md5[21], ""):
                print "<filt-vcf>  ---------: Passed"
            else:
                print "<filt-vcf>  ---------: Failed"
                sys.exit("bedtools sort or awk failed, please ensure bedtools and awk are functional")
            if checker(md5[22], ""):
                print "<svoverlap> ---------: Passed"
            else:
                print "<svoverlap> ---------: Failed"
                sys.exit("awk or sed failed, please ensure awk and sed are functional")
            os.system('echo "AllPassed" > ' + maindir + '/test/AllPassed')

if __name__ == "__main__":
    main()
