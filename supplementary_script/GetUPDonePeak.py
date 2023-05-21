# -*- coding: utf-8 -*-
"""
Created on 2022\05\13

Author: liugaojing2012

Email:liugaojing2012@foxmal.com
"""

import argparse
import os
import re
import multiprocessing
import glob

parser = argparse.ArgumentParser(description="Get Down or up cuttag peak")
parser.add_argument('--InputFile',help='Inputfiles',required=True)
parser.add_argument('--OutPutPath',help='Out put path',required=True)
parser.add_argument('-p','--pvalueTh',help='Threshold p value, default=0.05',default=0.05)
parser.add_argument('-q','--qvalueTh',help='Threshold q value, p value adjust, default=0.05',default=0.05)
parser.add_argument('--foldChangeTh',help='Threshold of log2(fold change)',default=1)
parser.add_argument('--PeakPosition',help='Peak Position file',required=True)
parser.add_argument('--EnsemblTOsymbl',help='Ensembl to symbl',required=True)
parser.add_argument('--DeseqFile',help='Deseq2 file that contain pvalue and q value',required=True)

argv=vars(parser.parse_args())

InputFile    = argv['InputFile']
OutPutPath   = argv['OutPutPath']
pvalueTh     = argv['pvalueTh']
qvalueTh     = argv['qvalueTh']
PeakPosition = argv['PeakPosition']
EnsemblTOsymbl = argv['EnsemblTOsymbl']
DeseqFile    = argv['DeseqFile']
foldChangeTh = argv['foldChangeTh']

foutUP = open(OutPutPath+"/UP.sig.different.CutTag.region.txt","w+")
foutDOWN = open(OutPutPath+"/DOWN.sig.different.CutTag.region.txt","w+")
foutNS = open(OutPutPath+"/NO.sig.different.CutTag.region.txt","w+")

dictPeakPosition = {} #{1_28940_29328:Peak0} #Not bed position
for eachline in open(PeakPosition,"r"):
    peakNumber = eachline.strip().split("\t")[0]
    chrName    = eachline.strip().split("\t")[1]
    start      = eachline.strip().split("\t")[2]
    end        = eachline.strip().split("\t")[3]
    key = chrName+"_"+start+"_"+end
    if key in dictPeakPosition.keys():
        print("Script:GetUPDonePeak.py error, the peak dup")
    else:
        dictPeakPosition[key] = peakNumber

ensembleIDGeneName = {} #{ENSG00000143514:TP53BP2}
for eachline in open(EnsemblTOsymbl,"r"):
    ensembleID = eachline.strip().split("\t")[0]
    geneName   = eachline.strip().split("\t")[1]
    ensembleIDGeneName[ensembleID] = geneName

dictPeakTh = {} #{"Peak17947":[log2FoldChange,pvalue,padj]}
for eachdegline in open(DeseqFile,"r"):
    if eachdegline.strip().split("\t")[0] == "baseMean":
        continue
    else:
        peaknumber = eachdegline.strip().split("\t")[0]
        log2FoldChange = eachdegline.strip().split("\t")[2]
        pvalue = eachdegline.strip().split("\t")[5]
        qvalue = eachdegline.strip().split("\t")[6]
        dictPeakTh[peaknumber] = [log2FoldChange,pvalue,qvalue]
#print(dictPeakTh)
for eachline in open(InputFile,"r"):
    seqnames = eachline.strip().split("\t")[0]
    start    = eachline.strip().split("\t")[1]
    end      = eachline.strip().split("\t")[2]
    ensembl  = eachline.strip().split("\t")[11]
    keyfile = seqnames+"_"+start+"_"+end
    if seqnames == "seqnames":
        header = eachline.strip().split("\t")+["peakNumber","log2FoldChange","pvalue","padj","genesymbl"]
        foutUP.write("\t".join(header)+"\n")
        foutDOWN.write("\t".join(header)+"\n")
        foutNS.write("\t".join(header)+"\n")
    else:
        if keyfile in dictPeakPosition.keys():
            peakNumber = dictPeakPosition[keyfile]
#            print(peakNumber)
#            print(dictPeakTh[peakNumber])
            if dictPeakTh[peakNumber][0] == "NA" or dictPeakTh[peakNumber][1] == "NA" or dictPeakTh[peakNumber][2] == "NA":
                genesymbl = ensembleIDGeneName[ensembl]
                resultOut = eachline.strip().split("\t")+[peakNumber,dictPeakTh[peakNumber][0],dictPeakTh[peakNumber][1],dictPeakTh[peakNumber][2],genesymbl]
                foutNS.write("\t".join(resultOut)+"\n")
            else:
                log2FoldChange = float(dictPeakTh[peakNumber][0])
                pvalue = float(dictPeakTh[peakNumber][1])
                padj   = float(dictPeakTh[peakNumber][2])
                genesymbl = ensembleIDGeneName[ensembl]
                resultOut = eachline.strip().split("\t")+[peakNumber,str(log2FoldChange),str(pvalue),str(padj),genesymbl]
                if log2FoldChange > float(foldChangeTh) and  pvalue < float(pvalueTh) and padj < float(qvalueTh):
                    foutUP.write("\t".join(resultOut)+"\n")
                elif log2FoldChange < -float(foldChangeTh) and  pvalue < float(pvalueTh) and padj < float(qvalueTh):
                    foutDOWN.write("\t".join(resultOut)+"\n")
                else:
                    foutNS.write("\t".join(resultOut)+"\n")
        else:
            print("Script:GetUPDonePeak.py error, the key error")

foutUP.close()
foutDOWN.close()
foutNS.close()
