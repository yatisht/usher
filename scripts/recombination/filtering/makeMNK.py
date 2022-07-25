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


def makeMNK():
    myOutString = ''
    with open('filtering/data/allRelevantNodesInfSeq.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            # seq = BAABAAAABBABBBBAAAABBB
            seq = splitLine[-1]
            if seq.startswith('A'):
                myOutString += joiner(splitLine)+'\t'+joiner([seq.count('A'),seq.count('B'),getK(seq,'A','B')])+'\n'
            else:
                myOutString += joiner(splitLine)+'\t'+joiner([seq.count('B'),seq.count('A'),getK(seq,'B','A')])+'\n'
    open('filtering/data/allRelevantNodesMNK.txt','w').write(myOutString)

def removeDups():
    alreadyDone = {}
    myOutString = ''
    with open('filtering/data/allRelevantNodesMNK.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not str(splitLine[-3])+'_'+str(splitLine[-2])+'_'+str(splitLine[-1]) in alreadyDone:
                myOutString += str(splitLine[-3])+' '+str(splitLine[-2])+' '+str(splitLine[-1])+'\n'
            alreadyDone[str(splitLine[-3])+'_'+str(splitLine[-2])+'_'+str(splitLine[-1])] = True
    open('filtering/data/mnk_no_dups.txt','w').write(myOutString)

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
    makeMNK()
    removeDups()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit

