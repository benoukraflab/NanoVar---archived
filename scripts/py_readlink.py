"""
This script prints the full path of files.

Copyright (C) 2019 Tham Cheng Yong, Roberto Tirado Magallanes, Touati Benoukraf.

This file is part of NanoVar.
"""

from sys import argv
import sys
import os

if len(argv)<2 or len(argv)>=3:
    sys.exit("Usage: python py_readlink.py filepath")

print os.path.realpath(sys.argv[1])