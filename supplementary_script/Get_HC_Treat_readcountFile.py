# -*- coding: utf-8 -*-
"""
Created on 2021\11\08

@author: liugaojing2012

Email:liugaojing2012@foxmail.com

"""
import argparse
import os
import re
import multiprocessing
import glob

rootdir=os.getcwd()
parser = argparse.ArgumentParser(description="Get HC and treatment read count file ")
#parser.add_argument('--conditionFile',help='conditionFile',required=True)
parser.add_argument('--strHCsample',help='All HC sample,string',required=True)
parser.add_argument('--strTreatmentsample',help='All treatment sample,string',required=True)
parser.add_argument('--targetPath',help='Target path',required=True)
parser.add_argument('--outPath',help='output path',required=True)
argv=vars(parser.parse_args())
#conditionFile  = argv['conditionFile']
strHCsample    = argv['strHCsample']
strTreatmentsample = argv['strTreatmentsample']
targetPath     = argv['targetPath']
outPath        = argv['outPath']

#list_HC = [line.strip().split("\t")[0] for line in open(conditionFile,"r") if line.strip().split("\t")[1] == "HC"]
#list_Treat = [line.strip().split("\t")[0] for line in open(conditionFile,"r") if line.strip().split("\t")[1] == "Treat"]

list_HC = strHCsample.strip().split(",")
list_Treat = strTreatmentsample.strip().split(",")

dict_HC,dict_Treat = {},{}  #{sampleName:{Gene1:readcount,Gene2:readcount......}}
for eachsample in list_HC:
    for line in open(targetPath+"/"+eachsample+".readcount.txt","r"):
        if line.strip().strip("\n").startswith("#") or line.strip().strip("\n").split("\t")[0] == "Geneid":
            continue
        else:
            list_line = line.strip().split('\t')
            GeneID    = list_line[0]
            ReadCount = list_line[-1]
            if eachsample in dict_HC.keys():
                dict_HC[eachsample][GeneID] = ReadCount
            else:
                dict_HC[eachsample] = {GeneID:ReadCount}
for eachsample in list_Treat:
    for line in open(targetPath+"/"+eachsample+".readcount.txt","r"):
        if line.strip().strip("\n").startswith("#") or line.strip().strip("\n").split("\t")[0] == "Geneid":
            continue
        else:
            list_line = line.strip().split('\t')
            GeneID    = list_line[0]
            ReadCount = list_line[-1]
            if eachsample in dict_Treat.keys():
                dict_Treat[eachsample][GeneID] = ReadCount
            else:
                dict_Treat[eachsample] = {GeneID:ReadCount}
list_all_gene = dict_Treat[eachsample].keys()
foutHC    = open(outPath+"/"+"HC.readcount.txt","w+")
foutHC.write("Gene"+"\t"+"\t".join(list_HC)+"\n")
for eachGene in list_all_gene:
    list_write = []
    for eachsample in list_HC:
        list_write.append(dict_HC[eachsample][eachGene])
    foutHC.write(eachGene+"\t"+"\t".join(list_write)+"\n")

foutTreat = open(outPath+"/"+"Treat.readcount.txt","w+")
foutTreat.write("Gene"+"\t"+"\t".join(list_Treat)+"\n")
for eachGene in list_all_gene:
    list_write = []
    for eachsample in list_Treat:
        list_write.append(dict_Treat[eachsample][eachGene])
    foutTreat.write(eachGene+"\t"+"\t".join(list_write)+"\n")
foutTreat.close()
foutHC.close()
