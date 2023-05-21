# -*- coding: utf-8 -*-
"""
Created on 2023\05\12

@author: liugaojing2012

Email:liugaojing2012@foxmail.com

Any questions email me.

"""

import argparse
import pysam
import os
import re
import operator
import numpy
import multiprocessing
import glob

parser = argparse.ArgumentParser(description="This script is used for transforming format from narrowpeak to saf")
parser.add_argument('--In',help='Input file',required=True)
parser.add_argument('--Out',help='Output file',required=True)

argv=vars(parser.parse_args())
In  = argv['In']
Out = argv['Out']

fout = open(Out,"w+")

PeakNumber = 0
for eachline in open(In,"r"):
    chrName = eachline.strip().split("\t")[0]
    start   = str(int(eachline.strip().split("\t")[1])+1)
    end     = eachline.strip().split("\t")[2]
    fout.write("Peak"+str(PeakNumber)+"\t"+chrName+"\t"+start+"\t"+end+"\t"+"."+"\n")
    PeakNumber += 1
fout.close()
