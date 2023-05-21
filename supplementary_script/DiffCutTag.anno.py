# -*- coding: utf-8 -*-
"""
Created on 2023\05\13

Author: liugaojing2012

Email:liugaojing2012@foxmal.com
"""

import argparse
import os
import re
import multiprocessing
import glob

parser = argparse.ArgumentParser(description="Chansfer peak to bed and annotation")
parser.add_argument('--InputFile',help='Inputfiles',required=True)
parser.add_argument('--OutPutPath',help='Out put path',required=True)
parser.add_argument('--PeakPosition',help='Peak position file',required=True)

argv=vars(parser.parse_args())

InputFile    = argv['InputFile']
OutPutPath   = argv['OutPutPath']
PeakPosition = argv['PeakPosition']

dictPeakPosition = {} #{Peak0:[1,28940-1,29328]}
for eachline in open(PeakPosition,"r"):
    peakNumber = eachline.strip().split("\t")[0]
    chrName    = eachline.strip().split("\t")[1]
    start      = str(int(eachline.strip().split("\t")[2])-1)
    end        = eachline.strip().split("\t")[3]
    dictPeakPosition[peakNumber] = [chrName,start,end]

fout = open(OutPutPath,"w+")
for line in open(InputFile,"r"):
    if line.strip().split("\t")[0] == "baseMean":
        continue
    else:
        peakNumber = line.strip().split("\t")[0]
        chrName = dictPeakPosition[peakNumber][0]
        start   = dictPeakPosition[peakNumber][1]
        end     = dictPeakPosition[peakNumber][2]
        fout.write(chrName+"\t"+start+"\t"+end+"\n")

fout.close()











