# -*- coding: utf-8 -*-
"""
Created on 2023\04\17

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

parser = argparse.ArgumentParser(description="This script is used for analysis of Atac-seq")
parser.add_argument('-cfg','--samplecfg',help='sample congigure file',required=True)
parser.add_argument('-qc','--fastqc',help='fastqc',default='/home/liugaojing/soft/fastqc/FastQC/fastqc')
parser.add_argument('-p','--pvalue',help='DEG threshold p value, default=0.05',default=0.05)
parser.add_argument('-q','--qvalue',help='DEG threshold q value, p value adjust, default=0.05',default=0.05)
parser.add_argument('-fdc','--foldchange',help='DEG threshold expression fold change(log2(folachange)), default=1',default=1)
parser.add_argument('--dataType',help='cleanData or rawData, default=cleanData',default="cleanData")
parser.add_argument('--species',help='Human or Rat or Zokor,default=Human',default="Human")
parser.add_argument('--conditionCompare','-CC',help='Condition files contain paired condition which used to get different atac region',required=True)
parser.add_argument('--IgGsample',help='The text file contain sample and IgG blan sample',required=True)

argv=vars(parser.parse_args())
samplecfg        = argv['samplecfg']
fastqc           = argv['fastqc']
pvalue           = argv['pvalue']
qvalue           = argv['qvalue']
foldchange       = argv['foldchange']
dataType         = argv['dataType']
species          = argv['species']
conditionCompare = argv['conditionCompare']
IgGsample        = argv['IgGsample']

rootdir = os.getcwd()

class BulkCutTag():
    def __init__(self):
        self.rootdir = os.getcwd()
        self.log = open("LogFile.txt","w+")

        self.dictTypeSample = {} #{"condit1":[sample1,sample2,sample3]}
        self.listAllSample = []  #[sample1,sample2,sample3]
        self.listAllType = []    #[condition1,condition2,condition3...]
        self.dictIgG     = {}    #{"sample1":"IgGsample"}   {"SRR12246717":"SRR11923224"}
        for eachline in open(samplecfg,"r"):
            sampleName = eachline.strip().split("\t")[0]
            sampleType = eachline.strip().split("\t")[2]
            self.listAllSample.append(sampleName)
            self.listAllType.append(sampleType)
            if sampleType in self.dictTypeSample.keys():
                self.dictTypeSample[sampleType].append(sampleName)
            else:
                self.dictTypeSample[sampleType] = [sampleName]
        print(self.dictTypeSample)

        for eachIgGline in open(IgGsample,"r"):
            if eachIgGline.strip().split("\t")[0] == "Sample":
                continue
            else:
                self.dictIgG[eachIgGline.strip().split("\t")[0]] = eachIgGline.strip().split("\t")[1]

        self.list_compair = []
        for eachconditionline in open(conditionCompare,"r"):
            if eachconditionline.strip().split("\t")[0] == "Control":
                continue
            else:
                condition1 = eachconditionline.strip().split("\t")[0]
                condition2 = eachconditionline.strip().split("\t")[1]
                self.list_compair.append(condition1+"_"+condition2)
        if species == "HumanHG38":
            self.gtf = "/home/liugaojing/DataBase/Genome/Human/ensemble/hg38_sm_primary_assembly/release-105/Homo_sapiens.GRCh38.105.gtf"
            self.genome = "/home/liugaojing/DataBase/Genome/Human/ensemble/hg38_sm_primary_assembly/release-105/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
            self.bowtie2Genome = "/home/liugaojing/DataBase/Genome/Human/ensemble/hg38_sm_primary_assembly/release-105/bowtie2Index/hg38"
            self.genomeSize = "/home/liugaojing/DataBase/Genome/Human/ensemble/hg38_sm_primary_assembly/release-105/Human38.chrom.sizes"
            #samtools faidx self.genome get genome.fa.fai  
            #cut -f 1,2 genome.fa.fai > Human38.chrom.sizes
            self.GenomeEffectiveGenome = "2.9e9"
            self.TxDB = "/home/liugaojing/DataBase/Genome/Human/ensemble/hg38_sm_primary_assembly/release-105/TxDb/TxDbHuman38.FromGtf.sqlite"
            self.geneTransFile = "/home/liugaojing/DataBase/Genome/Human/ensemble/hg38_sm_primary_assembly/release-105/ensembleID_geneName/EnsembleID_geneName.txt"
            self.hintName = "human38"
            self.genePosition = "/home/liugaojing/DataBase/Genome/Human/ensemble/hg38_sm_primary_assembly/release-105/GenePosition/GenePosition.bed"
        if species == "Mmus":
            self.gtf = "/home/liugaojing/DataBase/Genome/Mus_musculus/release105/Mus_musculus.GRCm39.105.gtf"
            self.genome = "/home/liugaojing/DataBase/Genome/Mus_musculus/release105/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa"
            self.bowtie2Genome = "/home/liugaojing/DataBase/Genome/Mus_musculus/release105/bowtie2Index/MouseGenome"
            self.genomeSize = "/home/liugaojing/DataBase/Genome/Mus_musculus/release105/Mouse.chrom.sizes"
            self.GenomeEffectiveGenome = "2.6e9"
            self.TxDB = "/home/liugaojing/DataBase/Genome/Mus_musculus/release105/TxDb/TxDbMmus.FromGtf.sqlite"
            self.geneTransFile = "/home/liugaojing/DataBase/Genome/Mus_musculus/release105/ensembleID_geneName/EnsembleID_geneName.txt"
            self.hintName = "Mmus"
            self.genePosition = "/home/liugaojing/DataBase/Genome/Mus_musculus/release105/GenePosition/GenePosition.bed"

    def lnData(self,sampleName,Datapath,filewrite):
        sampleFiles = sorted(glob.glob(Datapath+"/"+sampleName+"*.fq.gz"))
        print (sampleFiles)
        if len(sampleFiles) != 2:
            self.log.write("lnData:This sample read1 and read2 cant paired"+"\n")
        else:
            cmd1 = "ln -s "+sampleFiles[0]+" "+self.rootdir+"/data/"+sampleName+".R1.fq.gz"
            cmd2 = "ln -s "+sampleFiles[1]+" "+self.rootdir+"/data/"+sampleName+".R2.fq.gz"
            filewrite.write(cmd1+"\n")
            filewrite.write(cmd2+"\n")

    def lnCleanData(self,sampleName,Datapath,outPath,filewrite):
        sampleFiles = sorted(glob.glob(Datapath+"/"+sampleName+"*.fq.gz"))
        print (sampleFiles)
        if len(sampleFiles) != 2:
            self.log.write("lnData:This sample read1 and read2 cant paired"+"\n")
        else:
            cmd1 = "ln -s "+sampleFiles[0]+" "+outPath+"/"+sampleName+"_fastpTrimmed.R1.fq.gz"
            cmd2 = "ln -s "+sampleFiles[1]+" "+outPath+"/"+sampleName+"_fastpTrimmed.R2.fq.gz"
            filewrite.write(cmd1+"\n")
            filewrite.write(cmd2+"\n")

    def FastQC(self,sampleName,Datapath,filewrite,fileOutPath):
        sampleFiles = [Datapath+"/"+sampleName+".R1.fq.gz",
                       Datapath+"/"+sampleName+".R2.fq.gz"]
        if len(sampleFiles) != 2:
            self.log.write("FastQC:This sample read1 and read2 cant paired"+"\n")
        else:
            cmd = fastqc+" --format fastq --noextract --outdir "+fileOutPath+" "+sampleFiles[0]+" "+sampleFiles[1]
            filewrite.write(cmd+"\n")

    def Multiqc(self,TargetPath,OutPutPath,filewrite):
        cmd = "/home/liugaojing/anaconda3/bin/multiqc "+TargetPath+"/"+" -o "+OutPutPath
        filewrite.write(cmd+"\n")

    def fastp(self,sampleName,Datapath,filewrite,fileOutPath):
        #CUTTAG office suggestion:There is no need to trim reads from out standard 25x25 PE sequencing, as adapter sequences will not be included in reads of inserts >25 bp. 
        #However, for users performing longer sequencing, reads will need to be trimmed by Cutadap
        sampleFiles = [Datapath+"/"+sampleName+".R1.fq.gz",
                       Datapath+"/"+sampleName+".R2.fq.gz"]
        if len(sampleFiles) != 2:
            self.log.write("Fastp:This sample read1 and read2 cant paired"+"\n")
        else:
            cmd = "/home/liugaojing/soft/fastp/fastp"+" --in1 "+sampleFiles[0]+" --in2 "+sampleFiles[1]+" --out1 "+fileOutPath+"/"+sampleName+"_fastpTrimmed.R1.fq.gz --out2 "+fileOutPath+"/"+sampleName+"_fastpTrimmed.R2.fq.gz --json "+fileOutPath+"/"+sampleName+".result.json --html "+fileOutPath+"/"+sampleName+".result.html --length_required 20 --qualified_quality_phred 20"
            filewrite.write(cmd+"\n")

    def trim_fastqc(self,filepath,name1,name2,filewrite,fileOutPath):
        cmd = fastqc+" --format fastq --noextract --outdir "+fileOutPath+" "+filepath+"/"+name1+" "+filepath+"/"+name2
        filewrite.write(cmd+"\n")

    def Multiqc(self,TargetPath,OutPutPath,filewrite):
        cmd = "/home/liugaojing/anaconda3/bin/multiqc "+TargetPath+"/"+" -o "+OutPutPath
        filewrite.write(cmd+"\n")

    def filteMycoplasma(self,sampleName,filewrite,outPath,targetPath):
        sampleFiles = [targetPath+"/"+sampleName+"_fastpTrimmed.R1.fq.gz",
                       targetPath+"/"+sampleName+"_fastpTrimmed.R2.fq.gz"]
        cmd = "/home/liugaojing/soft/bowtie2/bowtie2-2.4.1-linux-x86_64/bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 10 -p 10 -x /home/liugaojing/DataBase/Genome/Mycoplasma_hyorhinis/bowtie2Index/Myco"+" -1 "+sampleFiles[0]+" -2 "+sampleFiles[1]+" -S "+outPath+"/"+sampleName+".sam"+" --un-conc-gz "+outPath+"/"+sampleName+".filteMycoplasma.fq.gz"
        filewrite.write(cmd+"\n")

    def alignment(self,sampleName,filewrite,outPath,targetPath):
        #CUT&Tag office suggestion:There is no need to trim reads from out standard 25x25 PE sequencing, as adapter sequences will not be included in reads of inserts >25 bp. However, for users performing longer sequencing, reads will need to be trimmed by Cutadapt and mapped by --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 to ignore any remaining adapter sequence at the 3’ ends of reads during mapping.

        sampleFiles = [targetPath+"/"+sampleName+".filteMycoplasma.fq.1.gz",
                       targetPath+"/"+sampleName+".filteMycoplasma.fq.2.gz"]
#        cmd = "/home/liugaojing/soft/bowtie2/bowtie2-2.4.1-linux-x86_64/bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 10 -X 700 -p 10 -x "+self.bowtie2Genome+" -1 "+sampleFiles[0]+" -2 "+sampleFiles[1]+" -S "+outPath+"/"+sampleName+".sam"
        cmd = "/home/liugaojing/soft/bowtie2/bowtie2-2.4.1-linux-x86_64/bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 10 -x "+self.bowtie2Genome+" -1 "+sampleFiles[0]+" -2 "+sampleFiles[1]+" -S "+outPath+"/"+sampleName+".sam"
        filewrite.write(cmd+"\n")

    def deDup(self):
        suggestion = "CUT&Tag integrates adapters into DNA in the vicinity of the antibody-tethered pA-Tn5, and the exact sites of integration are affected by the accessibility of surrounding DNA. For this reason fragments that share exact starting and ending positions are expected to be common, and such ‘duplicates’ may not be due to duplication during PCR. In practice, we have found that the apparent duplication rate is low for high quality CUT&Tag datasets, and even the apparent ‘duplicate’ fragments are likely to be true fragments. Thus, we do not recommend removing the duplicates."
        result = "Thus, we do not recommend removing the duplicates"

    def samtools(self,sampleName,filewrite,outPath,targetPath):
        #eliminate all the alignment results that are below the minQualityScore defined by user.
        cmd1 = "/home/liugaojing/soft/samtools/samtools-1.14/bin/bin/samtools view -h -q 2 "+targetPath+"/"+sampleName+".sam" +" > "+ outPath+"/"+sampleName+".filetMycopy.filterLowQuality.sam"
        #Filter and keep the mapped read pairs
        cmd2 = "/home/liugaojing/soft/samtools/samtools-1.14/bin/bin/samtools view -bS -F 0x04 "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.sam > "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.bam"
        #sam to bam
        cmd3 = "/home/liugaojing/soft/samtools/samtools-1.14/bin/bin/samtools sort "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.bam"+" -o "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.sorted.bam"
        cmd31 = "/home/liugaojing/soft/samtools/samtools-1.14/bin/bin/samtools index "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.sorted.bam"
        cmd32 = "/home/liugaojing/anaconda3/bin/bamCoverage -b "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.sorted.bam -o "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.sorted.bam.bw"
        cmd33 = "/home/liugaojing/soft/samtools/samtools-1.14/bin/bin/samtools sort -n "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.bam"+" -o "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.sortedByName.bam"
        #Convert into bed file format
        cmd4 = "/home/liugaojing/soft/bedtools/bedtools2-2.25.0/bin/bedtools bamtobed -i "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.sortedByName.bam -bedpe > "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.bed"
        #Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
        cmd5 = "awk '$1==$4 && $6-$2 < 1000 {print $0}' "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.bed > "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.clean.bed"
        #Only extract the fragment related columns
        cmd6 = "cut -f 1,2,6 "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.clean.bed | sort -k1,1 -k2,2n -k3,3n >" + outPath+"/"+sampleName+".filetMycopy.filterLowQuality.fragment.bed"
        #bed to bedgraph
        cmd7 = "/home/liugaojing/soft/bedtools/bedtools2-2.25.0/bin/bedtools genomecov -bg -i "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.fragment.bed -g "+self.genomeSize+" > "+outPath+"/"+sampleName+".filetMycopy.filterLowQuality.fragment.bedgraph"

        filewrite.write(cmd1+"\n")
        filewrite.write(cmd2+"\n")
        filewrite.write(cmd3+"\n")
        filewrite.write(cmd31+"\n")
        filewrite.write(cmd32+"\n")
        filewrite.write(cmd33+"\n")
        filewrite.write(cmd4+"\n")
        filewrite.write(cmd5+"\n")
        filewrite.write(cmd6+"\n")
        filewrite.write(cmd7+"\n")

    def MACS2(self,sampleName,IgGsampleName,filewrite,outPath,targetPath):
        cmd1 = "/home/liugaojing/anaconda3/envs/py37/bin/macs2 callpeak "+"-g "+str(self.GenomeEffectiveGenome)+" -f BAMPE --keep-dup all --bdg --call-summits -t "+targetPath+"/"+sampleName+".merge.sorted.bam -c "+targetPath+"/"+IgGsampleName+".merge.sorted.bam"+" -n "+sampleName+" --outdir "+outPath+" -q 0.1"
        cmd2 = "sort -k8,8nr "+outPath+"/"+sampleName+"_peaks.narrowPeak > "+outPath+"/"+sampleName+"_peaks.sortByPvalue.narrowPeak"
        filewrite.write(cmd1+"\n")
        filewrite.write(cmd2+"\n")

    def ConditionPeakMerge(self,conditionList,compairName,filewrite,outPath,targetPath,targetPathBam):
        listFile = []
        listFileBam = []
        for eachcondition in conditionList:
            fileName = targetPath+"/"+eachcondition+"_peaks.narrowPeak"
            listFile.append(fileName)
            
            for eachsample in self.dictTypeSample[eachcondition]:
                fileBam = targetPathBam+"/"+eachsample+".filetMycopy.filterLowQuality.sorted.bam"
                listFileBam.append(fileBam)
        cmd1 = "cat "+" ".join(listFile)+" > "+outPath+"/"+compairName+".mergeAll.peaks.narrowPeak"
        cmd2 = "sort -k1,1 -k2,2n "+outPath+"/"+compairName+".mergeAll.peaks.narrowPeak"+" > "+outPath+"/"+compairName+".mergeAll.peaks.sortedByPosition.narrowPeak"
        cmd3 = "bedtools merge -i "+outPath+"/"+compairName+".mergeAll.peaks.sortedByPosition.narrowPeak > "+outPath+"/"+compairName+".mergeAll.peaks.merge.narrowPeak"
        cmd4 = "python /home/liugaojing/Pipeline/CUTTAG/supplementary_script/NarrowPeak_to_SAF.py --In "+outPath+"/"+compairName+".mergeAll.peaks.merge.narrowPeak"+" --Out "+outPath+"/"+compairName+".mergeAll.peaks.merge.saf"

        filewrite.write(cmd1+"\n")
        filewrite.write(cmd2+"\n")
        filewrite.write(cmd3+"\n")
        filewrite.write(cmd4+"\n")
        for eachbamfile in listFileBam:
            sampleName = eachbamfile.strip().split("/")[-1].strip().split(".")[0]
            cmd5 = "/home/liugaojing/soft/subread/subread-2.0.2-Linux-x86_64/bin/featureCounts -T 2 -F SAF -a "+outPath+"/"+compairName+".mergeAll.peaks.merge.saf -p -B --countReadPairs"+" -o "+outPath+"/"+sampleName+".readcount.txt"+" -f "+eachbamfile
            filewrite.write(cmd5+"\n")

    def DESeq2(self,HCcondition,Treatcondition,filewrite,targetPath,outPath):
        cmd1 = "python /home/liugaojing/Pipeline/CUTTAG/supplementary_script/Get_HC_Treat_readcountFile.py --strHCsample "+",".join(self.dictTypeSample[HCcondition])+" --strTreatmentsample "+",".join(self.dictTypeSample[Treatcondition])+" --targetPath "+targetPath+" --outPath "+outPath
        cmd2 = "Rscript /home/liugaojing/Pipeline/RNA-seq/script/Deseq2/Deseq2.R --HCFile "+outPath+"/HC.readcount.txt --TreatFile "+outPath+"/Treat.readcount.txt --HCNumber "+str(len(self.dictTypeSample[HCcondition]))+" --TreatmentNumber "+str(len(self.dictTypeSample[Treatcondition]))+" --outPath "+outPath
        filewrite.write(cmd1+"\n")
        filewrite.write(cmd2+"\n")


    def IDR(self,sample1,sample2,condition,filewrite,outPath,targetPath):
        cmd1 = "idr --samples "+targetPath+"/"+sample1+"_peaks.sortByPvalue.narrowPeak "+targetPath+"/"+sample2+"_peaks.sortByPvalue.narrowPeak --input-file-type narrowPeak --rank p.value --output-file "+outPath+"/"+condition+".idr.peak.txt"+" --plot --log-output-file "+outPath+"/"+condition+".idr.log"
        filewrite.write(cmd1+"\n")

#    def DiffCUTTAGAnno(self,InputPath,OutPutPath,filewrite):
#        cmd1 = "python /home/liugaojing/Pipeline/Atac-seq/bulk-supple-script/DiffAtacRegion/AnnoPeak.1.py --InputPath "+InputPath+" --OutPutPath "+OutPutPath+" --GeneIDTrans "+self.geneTransFile
#        filewrite.write(cmd1+"\n")

    def HomerDiffCUTTAG(self,inputPath,outPutPath,compairName,filewrite):
        cmd1 =  'awk \'{print "Peak_"NR"\\t"$1"\\t"$2"\\t"$3"\\t""."}\' '+inputPath+"/"+compairName+"/UP.sig.different.Atac.region.bed > "+outPutPath+"/"+compairName+"/UP.sig.different.Atac.region.bed.saf"
        cmd2 =  'awk \'{print "Peak_"NR"\\t"$1"\\t"$2"\\t"$3"\\t""."}\' '+inputPath+"/"+compairName+"/DOWN.sig.different.Atac.region.bed > "+outPutPath+"/"+compairName+"/DOWN.sig.different.Atac.region.bed.saf"
        cmd3 = "/home/liugaojing/soft/homer/bin/findMotifsGenome.pl "+outPutPath+"/"+compairName+"/UP.sig.different.Atac.region.bed.saf "+self.genome+" "+outPutPath+"/"+compairName+"/HomerUPResult -p 4"
        cmd4 = "/home/liugaojing/soft/homer/bin/findMotifsGenome.pl "+outPutPath+"/"+compairName+"/DOWN.sig.different.Atac.region.bed.saf "+self.genome+" "+outPutPath+"/"+compairName+"/HomerDOWNResult -p 4" 
        filewrite.write(cmd1+"\n")
        filewrite.write(cmd2+"\n")
        filewrite.write(cmd3+"\n")
        filewrite.write(cmd4+"\n")    

    def DiffCutTagAnno(self,inputFile,outPutFile,PeakPosition,filewrite,outPutPath):
        cmd1 = "python /home/liugaojing/Pipeline/CUTTAG/supplementary_script/DiffCutTag.anno.py --InputFile "+inputFile+" --OutPutPath "+outPutFile+" --PeakPosition "+PeakPosition
        cmd2 = "Rscript /home/liugaojing/Pipeline/Atac-seq/bulk-supple-script/DiffAtacRegion/chipseekerAnno.r --PeakFile "+outPutFile+" --TxDB "+self.TxDB+" --OutPath "+outPutPath
        filewrite.write(cmd1+"\n")
        filewrite.write(cmd2+"\n")

    def getUPDownPeak(self,inputFile,OutPutPath,pvalue,qvalue,FoldChange,PeakPosition,Deseq2File,filewrite):
        cmd1 = "python /home/liugaojing/Pipeline/CUTTAG/supplementary_script/GetUPDonePeak.py --InputFile "+inputFile+" --OutPutPath "+OutPutPath+" --pvalueTh "+str(pvalue)+" --qvalueTh "+str(qvalue)+" --foldChangeTh "+str(FoldChange)+" --PeakPosition "+PeakPosition+" --EnsemblTOsymbl "+self.geneTransFile+" --DeseqFile "+Deseq2File
        filewrite.write(cmd1+"\n")

    def Visualization(self,allsamplestr,InPutPath,OutPutPath,scriptPath):
        os.system("python /home/liugaojing/Pipeline/CUTTAG/supplementary_script/PrepareComputeMatrix.py --allsamplestr "+allsamplestr+" --InPutPath "+InPutPath+" --OutPutPath "+OutPutPath+" --bedfile "+self.genePosition+" --scriptPath "+scriptPath)
        
    def MergeBam(self,OutPath,conditionName,filewrite):
        listIn = []
#        listBlankBam = []
        for eachsample in self.dictTypeSample[conditionName]:
            bamFile = self.rootdir+"/Bowtie2/"+eachsample+".filetMycopy.filterLowQuality.sorted.bam"
            listIn.append(bamFile)
#            BlankbamFile = self.rootdir+"/Bowtie2/"+self.dictIgG[eachsample]+".filetMycopy.filterLowQuality.sorted.bam"
#            listBlankBam.append()
        cmd1 = "/home/liugaojing/soft/samtools/samtools-1.14/bin/bin/samtools merge -@ 20 -o "+OutPath+"/"+conditionName+".merge.bam "+" ".join(listIn)
        cmd2 = "/home/liugaojing/soft/samtools/samtools-1.14/bin/bin/samtools sort -@ 20 "+OutPath+"/"+conditionName+".merge.bam -o "+OutPath+"/"+conditionName+".merge.sorted.bam "
#        cmd2 = "/home/liugaojing/soft/samtools/samtools-1.14/bin/bin/samtools merge -@ 4 -o "+OutPath+"/"+conditionName+".Blank.merge.bam "+" ".join(listBlankBam)
        filewrite.write(cmd1+"\n")
        filewrite.write(cmd2+"\n")

    def CreatFile(self):
        if not os.path.isdir(self.rootdir+"/data"):
            os.system("mkdir "+self.rootdir+"/data")

        if not os.path.isdir(self.rootdir+"/FirstFastqc"):
            os.system("mkdir "+self.rootdir+"/FirstFastqc")

        if not os.path.isdir(self.rootdir+"/FirstMultiQC"):
            os.system("mkdir "+self.rootdir+"/FirstMultiQC")

        if not os.path.isdir(self.rootdir+"/Fastp"):
            os.system("mkdir "+self.rootdir+"/Fastp")

        if not os.path.isdir(self.rootdir+"/SecondFastqc"):
            os.system("mkdir "+self.rootdir+"/SecondFastqc")

        if not os.path.isdir(self.rootdir+"/SecondMultiQC"):
            os.system("mkdir "+self.rootdir+"/SecondMultiQC")

        if not os.path.isdir(self.rootdir+"/Bowtie2"):
            os.system("mkdir "+self.rootdir+"/Bowtie2")

        if not os.path.isdir(self.rootdir+"/FilteMycoplasma"):
            os.system("mkdir "+self.rootdir+"/FilteMycoplasma")

        if not os.path.isdir(self.rootdir+"/Samtools"):
            os.system("mkdir "+self.rootdir+"/Samtools")

        if not os.path.isdir(self.rootdir+"/PeakCalling"):
            os.system("mkdir "+self.rootdir+"/PeakCalling")

        if not os.path.isdir(self.rootdir+"/DiffPeak"):
            os.system("mkdir "+self.rootdir+"/DiffPeak")

        if not os.path.isdir(self.rootdir+"/IDR"):
            os.system("mkdir "+self.rootdir+"/IDR")

        if not os.path.isdir(self.rootdir+"/Visualization"):
            os.system("mkdir "+self.rootdir+"/Visualization")

        if not os.path.isdir(self.rootdir+"/FirstMultiQC"):
            os.system("mkdir "+self.rootdir+"/FirstMultiQC")

        if not os.path.isdir(self.rootdir+"/BamMerge"):
            os.system("mkdir "+self.rootdir+"/BamMerge")

detailObject = BulkCutTag()

detailObject.CreatFile()

dictTypeSample = {} #{"50K":[SRR891268,SRR891269,SRR891270,SRR891271]}
listCondition = []
listAllsampleAll = []
for line in open(samplecfg,"r"):
    sampleName = line.strip().split("\t")[0]
    samplePath = line.strip().split("\t")[1]
    condition = line.strip().split("\t")[-1]
    listCondition.append(condition)
    listAllsampleAll.append(sampleName)

    if condition in dictTypeSample.keys():
        dictTypeSample[condition].append(sampleName)
    else:
        dictTypeSample[condition] = [sampleName]

    fileOpen = open("First."+sampleName+".analysis.sh","w+")
    if dataType == "rawData":
        detailObject.lnData(sampleName,samplePath,fileOpen)
        detailObject.FastQC(sampleName,rootdir+"/data",fileOpen,rootdir+"/FirstFastqc")
        detailObject.fastp(sampleName,rootdir+"/data",fileOpen,rootdir+"/Fastp")
        detailObject.trim_fastqc(rootdir+"/Fastp",sampleName+"_fastpTrimmed.R1.fq.gz",sampleName+"_fastpTrimmed.R2.fq.gz",fileOpen,rootdir+"/SecondFastqc")
        detailObject.alignment(sampleName,fileOpen,rootdir+"/Bowtie2",rootdir+"/Fastp")
        detailObject.filterReads(sampleName,fileOpen,rootdir+"/Bowtie2",rootdir+"/Bowtie2")
        detailObject.RMdup(sampleName,fileOpen,rootdir+"/Bowtie2",rootdir+"/Bowtie2")
        detailObject.Tn5shift(sampleName,fileOpen,rootdir+"/TN5OffsetAdjustment",rootdir+"/Bowtie2")
        detailObject.BamtoBed(sampleName,fileOpen,rootdir+"/TN5OffsetAdjustment",rootdir+"/TN5OffsetAdjustment")
        detailObject.MACS2(sampleName,fileOpen,rootdir+"/PeakCalling",rootdir+"/TN5OffsetAdjustment")
        detailObject.sortNarrowPeak(sampleName,fileOpen,rootdir+"/PeakCalling",rootdir+"/PeakCalling")
        detailObject.bamToBW(sampleName,rootdir+"/TN5OffsetAdjustment",rootdir+"/TN5OffsetAdjustment",fileOpen)
        detailObject.sortNarrowByPosition(sampleName,fileOpen,rootdir+"/PeakCalling",rootdir+"/PeakCalling")
        detailObject.peakAnnoSinglePeak(sampleName,rootdir+"/PeakCalling",rootdir+"/PeakCalling",fileOpen)
        detailObject.homerSingleSample(rootdir+"/PeakCalling",rootdir+"/PeakCalling/"+sampleName+".HomerResult",sampleName,fileOpen)
        detailObject.HintAtacsingleSample(sampleName,rootdir+"/TFfootprinting/"+sampleName,rootdir+"/Bowtie2",rootdir+"/PeakCalling",fileOpen)
        
    elif dataType == "cleanData":
        detailObject.lnData(sampleName,samplePath,fileOpen)
        detailObject.FastQC(sampleName,rootdir+"/data",fileOpen,rootdir+"/FirstFastqc")
        detailObject.lnCleanData(sampleName,samplePath,rootdir+"/Fastp",fileOpen)
##        detailObject.trim_fastqc(rootdir+"/Fastp",sampleName+"_fastpTrimmed.R1.fq.gz",sampleName+"_fastpTrimmed.R2.fq.gz",fileOpen,rootdir+"/SecondFastqc")
        detailObject.filteMycoplasma(sampleName,fileOpen,rootdir+"/FilteMycoplasma",rootdir+"/Fastp")
        detailObject.alignment(sampleName,fileOpen,rootdir+"/Bowtie2",rootdir+"/FilteMycoplasma")
        detailObject.samtools(sampleName,fileOpen,rootdir+"/Bowtie2",rootdir+"/Bowtie2")
    else:
        print("Error sample Type")
    fileOpen.close()

#fileIDR = open("Third.IDR.analysis.sh","w+")
#for eachcondition in dictTypeSample.keys():
#    listsample_condition = dictTypeSample[eachcondition]
#    detailObject.IDR(listsample_condition[0],listsample_condition[1],eachcondition,fileIDR,rootdir+"/IDR",rootdir+"/PeakCalling")
#fileIDR.close()
if dataType == "rawData":
    fileshell = open("Second.MultiQc.analysis.sh","w+")
    detailObject.Multiqc(rootdir+"/FirstFastqc",rootdir+"/FirstMultiQC",fileshell)
    detailObject.Multiqc(rootdir+"/SecondFastqc",rootdir+"/SecondMultiQC",fileshell)
    fileshell.close()
if dataType == "cleanData":
    fileshell = open("Second.MultiQc.analysis.sh","w+")
    detailObject.Multiqc(rootdir+"/FirstFastqc",rootdir+"/FirstMultiQC",fileshell)
#    detailObject.Multiqc(rootdir+"/SecondFastqc",rootdir+"/SecondMultiQC",fileshell)
    fileshell.close()

fThird = open("Second.MergeSameConditionBam.analysis.sh","w+")
for eachcondition in list(set(listCondition)):
#    if not os.path.isdir(rootdir+"/BamMerge/"+eachcondition):
#        os.system("mkdir "+rootdir+"/BamMerge/"+eachcondition)
    detailObject.MergeBam(rootdir+"/BamMerge",eachcondition,fThird)
fThird.close()

fThird = open("Third.Peakcalling.analysis.sh","w+")
fThird.write("source /home/liugaojing/anaconda3/etc/profile.d/conda.sh"+"\n")
fThird.write("conda activate py37"+"\n")
for eachIgGline in open(IgGsample,"r"):
    sampleName = eachIgGline.strip().split("\t")[0]
    if sampleName == "Sample":
        continue
    else:
        condition1Name = eachIgGline.strip().split("\t")[0]
        IgGName = eachIgGline.strip().split("\t")[1]
        detailObject.MACS2(condition1Name,IgGName,fThird,rootdir+"/PeakCalling",rootdir+"/BamMerge")
fThird.write("conda activate base")
fThird.close()

conditionPairList = []
for eachconditionline in open(conditionCompare,"r"):
    if eachconditionline.strip().split("\t")[0] == "Control":
        continue
    else:
        condition1 = eachconditionline.strip().split("\t")[0]
        condition2 = eachconditionline.strip().split("\t")[1]
        conditionPairList.append(condition1+"_"+condition2)
for eachPair in conditionPairList:
    condition1 = eachPair.strip().split("_")[0]
    condition2 = eachPair.strip().split("_")[1]
#    listsample_condition1 = dictTypeSample[condition1]
#    listsample_condition2 = dictTypeSample[condition2]

    compairName = str(condition1)+"_"+str(condition2)
    if not os.path.isdir(rootdir+"/DiffPeak/"+compairName):
        os.system("mkdir "+rootdir+"/DiffPeak/"+compairName)
#    if not os.path.isdir(rootdir+"/DiffPeak/"+compairName+"/HomerDOWNResult"):
#        os.system("mkdir "+rootdir+"/DiffPeak/"+compairName+"/HomerDOWNResult")

    fileOpenpair = open("Fourth."+compairName+".diff.CutTag.analysis.sh","w+")
#    listAllsampleCondition = listsample_condition1 + listsample_condition2
    listAllsampleCondition = [condition1,condition2]
    detailObject.ConditionPeakMerge(listAllsampleCondition,compairName,fileOpenpair,rootdir+"/DiffPeak/"+compairName,rootdir+"/PeakCalling",rootdir+"/Bowtie2")
    detailObject.DESeq2(condition1,condition2,fileOpenpair,rootdir+"/DiffPeak/"+compairName,rootdir+"/DiffPeak/"+compairName)
    detailObject.DiffCutTagAnno(rootdir+"/DiffPeak/"+compairName+"/Deseq2_HC_vs_Treat_diff.xls",rootdir+"/DiffPeak/"+compairName+"/All.DESeq2.Peask.bed",rootdir+"/DiffPeak/"+compairName+"/"+compairName+".mergeAll.peaks.merge.saf",fileOpenpair,rootdir+"/DiffPeak/"+compairName)
    detailObject.getUPDownPeak(rootdir+"/DiffPeak/"+compairName+"/PeakAnno.result.tsv",rootdir+"/DiffPeak/"+compairName,0.05,1,0.585,rootdir+"/DiffPeak/"+compairName+"/"+compairName+".mergeAll.peaks.merge.saf",rootdir+"/DiffPeak/"+compairName+"/Deseq2_HC_vs_Treat_diff.xls",fileOpenpair)
    
    fileOpenpair.close()

detailObject.Visualization(",".join(listAllsampleAll),rootdir+"/Bowtie2",rootdir+"/Visualization",rootdir)
