# -*- coding: utf-8 -*-
"""
Created on 2022\05\14

Author: liugaojing2012

Email:liugaojing2012@foxmal.com
"""

import argparse
import os
import re
import multiprocessing
import glob

parser = argparse.ArgumentParser(description="Preparing computeMatrix script")
parser.add_argument('--allsamplestr',help='all sample name,sep by ,',required=True)
parser.add_argument('--InPutPath',help='In put path',required=True)
parser.add_argument('--OutPutPath',help='Out put path',required=True)
parser.add_argument('--bedfile',help='bed file',required=True)
parser.add_argument('--scriptPath',help='script path',required=True)

argv=vars(parser.parse_args())

allsamplestr = argv['allsamplestr']
InPutPath    = argv['InPutPath']
OutPutPath   = argv['OutPutPath']
bedfile      = argv['bedfile']
scriptPath   = argv['scriptPath']

listSample = allsamplestr.strip().split(",")

liatBWfile,listBamfile = [],[]
for eachsample in listSample:
     BWfile  = InPutPath+"/"+eachsample+".filetMycopy.filterLowQuality.sorted.bam.bw"
     Bamfile = InPutPath+"/"+eachsample+".filetMycopy.filterLowQuality.sorted.bam"
     liatBWfile.append(BWfile)
     listBamfile.append(Bamfile)
cmd1 = "computeMatrix scale-regions -S "+" ".join(liatBWfile)+" -R "+bedfile+" --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o "+OutPutPath+"/matrix_gene.mat.computeMatrix.gz -p 30"
cmd2 = "plotHeatmap -m "+OutPutPath+"/matrix_gene.mat.computeMatrix.gz -out "+OutPutPath+"/computeMatrix.cuttag.pdf --sortUsing sum --plotFileFormat pdf"
cmd3 = "/home/liugaojing/anaconda3/bin/multiBamSummary bins -p 30 --bamfiles "+" ".join(listBamfile)+" --labels "+" ".join(listSample)+" -o "+OutPutPath+"/Correlationallsample.npz --outRawCounts "+OutPutPath+"/allsample.correlation.tab"
cmd4 = "plotCorrelation -in "+OutPutPath+"/Correlationallsample.npz"+" --corMethod pearson --removeOutliers --skipZeros --plotTitle "+'"pearson Correlation of Read Counts"'+" --whatToPlot heatmap --colorMap RdYlBu_r --plotNumbers -o "+OutPutPath+"/heatmap_SpearmanCorr_readCounts.pdf --outFileCorMatrix "+OutPutPath+"/SpearmanCorr_readCounts.allsample.tab"

fout = open(scriptPath+"/Visual.sh","w+")
fout.write(cmd1+"\n")
fout.write(cmd2+"\n")
fout.write(cmd3+"\n")
fout.write(cmd4+"\n")
fout.close()
#os.system("computeMatrix scale-regions -S "+" ".join(BWfile)+" -R "+bedfile+" --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o "+OutPutPath+"/matrix_gene.mat.gz -p 75")
#os.system("plotHeatmap "+OutPutPath+"/matrix_gene.mat.gz -out "+OutPutPath+"/Heatmap.cuttag.png --sortUsing sum")
#os.system("/home/liugaojing/anaconda3/bin/multiBamSummary bins --bamfiles "+" ".join(listBamfile)+" --labels "+" ".join(listSample)+" -o "+OutPutPath+"/allsample.npz --outRawCounts "+OutPutPath+"/allsample.tab")
#os.system("plotCorrelation -in "+OutPutPath+"/allsample.npz"+" --corMethod pearson --skipZeros --plotTitle "+'"pearson Correlation of Read Counts"'+" --whatToPlot heatmap --colorMap RdYlBu_r --plotNumbers -o "+OutPutPath+"/heatmap_SpearmanCorr_readCounts.png --outFileCorMatrix SpearmanCorr_readCounts.allsample.tab")

