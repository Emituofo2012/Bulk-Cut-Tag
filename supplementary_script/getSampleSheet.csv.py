# -*- coding: utf-8 -*-
"""
Created on 2023\01\01

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
import time

parser = argparse.ArgumentParser(description="This script is used for creating samplesheet for diffbind")
parser.add_argument('--condition1',help='condition1',required=True)
parser.add_argument('--condition2',help='condition2',required=True)
parser.add_argument('--samplecfg',help='sample congigure file',required=True)
parser.add_argument('--OutPath',help='Out put path',required=True)
parser.add_argument('--InputPath',help='In put path',required=True)

argv=vars(parser.parse_args())
condition1 = argv['condition1']
condition2 = argv['condition2']
samplecfg  = argv['samplecfg']
OutPath    = argv['OutPath']
InputPath  = argv['InputPath']

rootdir = os.getcwd()

dictTypeSample = {} #{"condit1":[sample1,sample2,sample3]}
listAllSample = []  #[sample1,sample2,sample3]
listAllType = []    #[condition1,condition2,condition3...]
for eachline in open(samplecfg,"r"):
    sampleName = eachline.strip().split("\t")[0]
    sampleType = eachline.strip().split("\t")[2]
    listAllSample.append(sampleName)
    listAllType.append(sampleType)
    if sampleType in dictTypeSample.keys():
        dictTypeSample[sampleType].append(sampleName)
    else:
        dictTypeSample[sampleType] = [sampleName]

fout = open(OutPath+"/SampleSheet.csv","w+")
fout.write("SampleID,Tissue,Factor,Condition,Treatment,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller"+"\n")
for eachsample in dictTypeSample[condition1]:
    number = dictTypeSample[condition1].index(eachsample)+1
    list_write = [eachsample,"NA","NA","NA",condition1,str(number),InputPath+"/Bowtie2/"+eachsample+".filetMycopy.filterLowQuality.sorted.bam","NA","NA",InputPath+"/PeakCalling/"+eachsample+"_peaks.narrowPeak","narrow"]
    fout.write(",".join(list_write)+"\n")

for eachsample2 in dictTypeSample[condition2]:
    number = dictTypeSample[condition2].index(eachsample2)+1
    list_write = [eachsample2,"NA","NA","NA",condition2,str(number),InputPath+"/Bowtie2/"+eachsample2+".filetMycopy.filterLowQuality.sorted.bam","NA","NA",InputPath+"/PeakCalling/"+eachsample2+"_peaks.narrowPeak","narrow"]
    fout.write(",".join(list_write)+"\n")
fout.close()
