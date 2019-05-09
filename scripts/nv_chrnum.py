"""
This script calculates the number of FASTA sequences in the reference file

Copyright (C) 2019 Tham Cheng Yong, Roberto Tirado Magallanes, Touati Benoukraf.

This file is part of NanoVar.
"""

from sys import argv

gensize = open(argv[1], 'r').read().split('\n')
lg = len(gensize) - 1
print lg