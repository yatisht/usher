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

"""
goal of this is to give him the information in a digestible form.
he needs:
- list of samples to look at
- list of recombinant nodes and their parents
- list of nearest sample sequences for each recombinant/parent node
- list of informative sites for each sample


^ recombinant and parents are in combinedCatOnlyBestGt1 -- give him that
^ add informative sites to closest samples file

"""

##########################
##### MAIN FUNCTIONS #####
##########################

def makeSampleInfo():
    nodeToSites = {}
    nodeToClosestSamples = {}
    with open('filtering/data/allRelevantNodesToDescendants.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            nodeToSites[int(splitLine[0][1:-1])] = {}
            nodeToClosestSamples[int(splitLine[0][1:-1])] = splitLine[1]

    nodeToRelatives = {}
    with open('filtering/data/combinedCatOnlyBestWithPVals.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                for i in [0,3,6]:
                    if not int(splitLine[i]) in nodeToRelatives:
                        nodeToRelatives[int(splitLine[i])] = {}
                    if i == 0:
                        nodeToRelatives[int(splitLine[i])][int(splitLine[3])] = True
                        nodeToRelatives[int(splitLine[i])][int(splitLine[6])] = True
                    elif i == 3:
                        nodeToRelatives[int(splitLine[i])][int(splitLine[0])] = True
                        nodeToRelatives[int(splitLine[i])][int(splitLine[6])] = True
                    elif i == 6:
                        nodeToRelatives[int(splitLine[i])][int(splitLine[3])] = True
                        nodeToRelatives[int(splitLine[i])][int(splitLine[6])] = True

    with open('filtering/data/allRelevantNodesInfSites.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            for k in splitLine[-1].split(','):
                if int(splitLine[0]) in nodeToSites:
                    nodeToSites[int(splitLine[0])][int(k)] = True
                if int(splitLine[1]) in nodeToSites:
                    nodeToSites[int(splitLine[1])][int(k)] = True
                if int(splitLine[2]) in nodeToSites:
                    nodeToSites[int(splitLine[2])][int(k)] = True

    myOutString = 'node\tdescendants\tinformative_sites\n'
    for n in sorted(nodeToSites.keys()):
        myOutString += str(n)+'\t'+nodeToClosestSamples[n]+'\t'+joinerC(nodeToSites[n].keys())+'\n'
    open('filtering/data/sampleInfo.txt','w').write(myOutString)

# def getSubset():
#     mySamples = {}
#     with open('last_batch.txt') as f:
#         for line in f:
#             mySamples[line.strip()] = True

#     myOutString = 'node\tdescendants\tinformative_sites\n'
#     with open('sampleInfoEdit.txt') as f:
#         for line in f:
#             toPrint = False
#             splitLine = (line.strip()).split('\t')
#             if not splitLine[0] == 'node':
#                 for k in splitLine[1].split(','):
#                     if k in mySamples:
#                         toPrint = True
#             if toPrint == True:
#                 myOutString += line.strip()+'\n'
#     open('sampleInfoEditSubset.txt','w').write(myOutString)

#     myTrios = {}
#     with open('full_report.txt') as f:
#         for line in f:
#             splitLine = (line.strip()).split('\t')
#             if splitLine[0] != 'recomb_node_id':
#                 if 'Missing' in line.strip():
#                     myTrios[str(splitLine[0])+'_'+str(splitLine[1])+'_'+str(splitLine[2])] = True

#     myOutString = ''
#     with open('combinedCatOnlyBestWithPVals.txt') as f:
#         for line in f:
#             splitLine = (line.strip()).split('\t')
#             if splitLine[0]+'_'+str(splitLine[3])+'_'+str(splitLine[6]) in myTrios:
#                 myOutString += line.strip()+'\n'
#     open('combinedCatOnlyBestWithPValsSubset.txt','w').write(myOutString)


##########################
#### HELPER FUNCTIONS ####
##########################

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
    makeSampleInfo()
    #getSubset()


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit

