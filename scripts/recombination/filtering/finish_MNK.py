#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 2/1/2018
# compareDatabases.py

import sys
import os
import datetime
import numpy
from numpy import random
import gzip
import math
import re

##########################
##### MAIN FUNCTIONS #####
##########################

#Run 3seq first

def addPVals():
    keyToP = {}
    with open('filtering/data/mnk.log') as f:
        for line in f:
            splitLine = (line.strip()).split()
            if len(splitLine) > 1:
                if splitLine[0].startswith('Enter'):
                    myKey = '_'.join(splitLine[-3:])
                else:
                    keyToP[myKey] = float(splitLine[-1])

    myOutString = ''
    alreadyUsed = {}
    with open('filtering/data/allRelevantNodesMNK.txt') as f:
        for line in f:
            splitLine = (line.strip()).split()
            myKey = '_'.join(splitLine[-3:])
            if not myKey in keyToP:
                print(myKey)
            else:
                myOutString += joiner(splitLine)+'\t'+str(keyToP[myKey])+'\n'
    open('filtering/data/allRelevantNodesMNKPval.txt','w').write(myOutString)


def combinePValueFiles():
    recombToParents = {}
    with open('filtering/data/combinedCatOnlyBestWithPVals.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                if not int(splitLine[0]) in recombToParents:
                    recombToParents[int(splitLine[0])] = {}
                recombToParents[int(splitLine[0])][str(splitLine[3])+'_'+str(splitLine[6])] = True

    recombTo3seqPval = {}
    recombToBestParents = {}
    recombToAB = {}
    recombToPrinted = {}
    with open('filtering/data/allRelevantNodesMNKPval.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                if int(splitLine[0]) in recombToParents:
                    if (str(splitLine[1])+'_'+str(splitLine[2]) in recombToParents[int(splitLine[0])]):
                        if int(splitLine[0]) not in recombTo3seqPval or float(splitLine[11]) < recombTo3seqPval[int(splitLine[0])]:
                            recombTo3seqPval[int(splitLine[0])] = float(splitLine[11])
                            recombToBestParents[int(splitLine[0])] = str(splitLine[1])+'_'+str(splitLine[2])
                            recombToAB[int(splitLine[0])] = splitLine[7]
                            recombToPrinted[int(splitLine[0])] = False

    myOutString = ''
    myOutString2 = ''
    with open('filtering/data/combinedCatOnlyBestWithPVals.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                if int(splitLine[0]) in recombTo3seqPval:
                    if recombToBestParents[int(splitLine[0])] == str(splitLine[3])+'_'+str(splitLine[6]):
                        myOutString += '\t'.join(splitLine[:-1])+'\t'+str(recombTo3seqPval[int(splitLine[0])])+'\t'+recombToAB[int(splitLine[0])]+'\t'+str(splitLine[-1])+'\n'
                        if splitLine[12].startswith('0/'):
                            splitLine[12] = '0.0'
                        if splitLine[13].startswith('0/'):
                            splitLine[13] = '0.0'
                        myOutString2 += '\t'.join(splitLine[:-1])+'\t'+str(recombTo3seqPval[int(splitLine[0])])+'\t'+recombToAB[int(splitLine[0])]+'\t'+str(splitLine[-1])+'\n'
                        recombToPrinted[int(splitLine[0])] = True
    open('filtering/data/combinedCatOnlyBestWithAll3PValsTiesBroken.txt','w').write(myOutString)
    open('filtering/data/combinedCatOnlyBestWithAll3PValsRealTiesBroken.txt','w').write(myOutString2)

    for k in recombToPrinted:
        if recombToPrinted[k] == False:
            print(k, recombToBestParents[k])

def addInfSites():
    finalReportTrios = {}
    with open('filtering/data/final_report.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            finalReportTrios[str(splitLine[0])+'_'+str(splitLine[1])+'_'+str(splitLine[2])] = True

    trioToInfSites = {}
    with open('filtering/data/allRelevantNodesInfSites.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            trioToInfSites[str(splitLine[0])+'_'+str(splitLine[1])+'_'+str(splitLine[2])] = splitLine[7]

    trioToInfSeq = {}
    with open('filtering/data/allRelevantNodesInfSeq.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            trioToInfSeq[str(splitLine[0])+'_'+str(splitLine[1])+'_'+str(splitLine[2])] = splitLine[7]

    myOutString = ''
    with open('filtering/data/combinedCatOnlyBestWithPVals.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6]) in trioToInfSites and str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6]) in finalReportTrios:
                splitLine.append(trioToInfSites[str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])])
                splitLine.append(trioToInfSeq[str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])])
                myOutString += joiner(splitLine)+'\n'
    open('filtering/data/combinedCatOnlyBestWithPValsFinalReportWithInfSites.txt','w').write(myOutString)



##########################
#### HELPER FUNCTIONS ####
##########################

def getK(seq, a, b):
    myPath = []
    currentPlace = 0
    for k in seq:
        if k == a:
            currentPlace += 1
        else:
            currentPlace -= 1
        myPath.append(currentPlace)
    maxDesc = 0
    for i in range(1,len(myPath)):
        if max(myPath[:i])-myPath[i] > maxDesc:
            maxDesc = max(myPath[:i])-myPath[i]
    return(maxDesc)


def getPos(myInds, intLineNumToPos):
    myReturn = []
    for k in myInds:
        myReturn.append(intLineNumToPos[k])
    return(myReturn)

def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return('\t'.join(newList))

def joinerC(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return(','.join(newList))

#########################
##### FUNCTION CALL #####
#########################

def main():
    addPVals()
    combinePValueFiles()
    addInfSites()


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit
