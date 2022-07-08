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

def checkClusters():
    trioTo3P = {}
    with open('filtering/data/allRelevantNodesMNKPval.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            trioTo3P[str(splitLine[0])+'_'+str(splitLine[1])+'_'+str(splitLine[2])] = splitLine[8:]

    goodTrios = {}
    with open('filtering/data/combinedCatOnlyBestWithPValsFinalReportWithInfSites.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            goodTrios[str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])] = False

    myOutString = ''
    bp1 = {}
    bp2 = {}
    myOKs = {28881:True,28882:True,28883:True,28280:True,28281:True,28282:True}
    counter = 0
    with open('filtering/data/combinedCatOnlyBestWithPValsFinalReportWithInfSites.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            myTrio = str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])
            myStart = int(splitLine[1].split(',')[0][1:])
            myEnd = int(splitLine[2].split(',')[1][:-1])
            mySeq = splitLine[16]
            mySites = toInt(splitLine[15].split(','))
            myA = []
            myB = []
            for i in range(0,len(mySites)):
                if mySeq[i] == 'A':
                    myA.append(mySites[i])
                elif mySeq[i] == 'B':
                    myB.append(mySites[i])
            if max(myA)-min(myA) > 20 and max(myB)-min(myB) > 20:
                myOutString += line.strip()+'\t'+joiner(trioTo3P[myTrio])+'\n'
                if myStart == 0 or myEnd == 29903:
                    bp1[splitLine[0]] = True
                else:
                    bp2[splitLine[0]] = True

            if str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6]) in goodTrios:
                counter += 1
                goodTrios[str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])] = True

    for k in goodTrios:
        if goodTrios[k] == False:
            print(k.split('_')[0])

    print(len(bp1),len(bp2))
    print(counter)
    open('filtering/data/combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClusters.txt','w').write(myOutString)


##########################
#### HELPER FUNCTIONS ####
##########################

def toInt(myList):
    myReturn = []
    for k in myList:
        myReturn.append(int(k))
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
    checkClusters()


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit
