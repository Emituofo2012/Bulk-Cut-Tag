import multiprocessing
import time
import os
import glob
import subprocess

allFileFirst = glob.glob("First*analysis.sh")
allFileSecond = glob.glob("Second*analysis.sh")
allFileThird = glob.glob("Third*analysis.sh")
allFileFourth = glob.glob("Fourth*analysis.sh")
allFileVisual = ["Visual.sh"]

def fuctionRun(shellName):
    sampleName = shellName.strip().split(".")[0]+"."+shellName.strip().split(".")[1]
    statas = subprocess.getstatusoutput("bash "+shellName+" >> "+sampleName+".log.txt 2>&1 ")
    if statas[0] == 0:
        pass
    else:
        print("ERROR:"+shellName)

def main():
    pool = multiprocessing.Pool(processes = 7)
    for eachfile in allFileFirst:
        pool.apply_async(fuctionRun, args=(eachfile,))
    pool.close()
    pool.join()

    pool = multiprocessing.Pool(processes = 5)
    for eachfile in allFileThird:
        pool.apply_async(fuctionRun, args=(eachfile,))
    pool.close()
    pool.join()

    pool = multiprocessing.Pool(processes = 5)
    for eachfile in allFileFourth:
        pool.apply_async(fuctionRun, args=(eachfile,))
    pool.close()
    pool.join()

    pool = multiprocessing.Pool(processes = 2)
    for eachfile in allFileSecond:
        pool.apply_async(fuctionRun, args=(eachfile,))
    pool.close()
    pool.join()

    pool = multiprocessing.Pool(processes = 1)
    for eachfile in allFileVisual:
        pool.apply_async(fuctionRun, args=(eachfile,))
    pool.close()
    pool.join()

if __name__ == "__main__":
    main()
