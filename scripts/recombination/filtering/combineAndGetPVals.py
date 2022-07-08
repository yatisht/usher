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

def catOnlyBest():
    nodeToLines = {}
    nodeToMinStart = {}
    with open('filtering/data/recombination.tsv') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):

                ### REPLACE TEXT
                if 'GENOME_SIZE' in splitLine[2]:
                    splitLine[2] = splitLine[2].replace('GENOME_SIZE', '29903')
                for i in [0,3,6]:
                    if 'node_' in splitLine[i]:
                        splitLine[i] = splitLine[i].replace('node_', '')

                if not str(splitLine[0]) in nodeToLines:
                    nodeToLines[str(splitLine[0])] = []
                    nodeToMinStart[str(splitLine[0])] = int(splitLine[-2])
                nodeToLines[str(splitLine[0])].append(splitLine)
                if int(splitLine[-2]) < nodeToMinStart[str(splitLine[0])]:
                    nodeToMinStart[str(splitLine[0])] = int(splitLine[-2])

    myOutString = ''
    for k in nodeToLines:
        for l in nodeToLines[k]:
            if int(l[-2]) > nodeToMinStart[k]:
                print(l)
                l[-2] = nodeToMinStart[k]
            myOutString += joiner(l)+'\n'
    open('filtering/data/catRecombinationReplacedMinStartingPars.tsv','w').write(myOutString)

    nodeToKeepLines = {}
    nodeToBestScore = {}
    with open('filtering/data/catRecombinationReplacedMinStartingPars.tsv') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                myImprovement = int(splitLine[-2])-int(splitLine[-1])
                myNode = str(splitLine[0])
                if myNode not in nodeToBestScore or nodeToBestScore[myNode] < myImprovement:
                    nodeToBestScore[myNode] = myImprovement
                    nodeToKeepLines[myNode] = ''
                if myImprovement == nodeToBestScore[myNode]:
                    nodeToKeepLines[myNode] += joiner(splitLine)+'\n'
    myOutString = ''
    for n in sorted(nodeToKeepLines.keys()):
        myOutString += nodeToKeepLines[n]
    open('filtering/data/catRecombOnlyBestScoresBeforeCombining.txt','w').write(myOutString)

    myFlag = ''
    for ITERATION in range(0, 10):
        if not myFlag == 'printed':
            recombToLines = {}
            lineCounter = 0
            if ITERATION == 0:
                ### Read in initial recombination output
                for line in myOutString.split('\n'):
                    splitLine = line.split('\t')
                    if len(splitLine) > 1:
                        lineCounter += 1
                        #splitLine[6] = splitLine[5][:-1] # correct for missed tab, retain only relevant fields
                        if not str(splitLine[0]) in recombToLines:
                            recombToLines[str(splitLine[0])] = []
                        recombToLines[str(splitLine[0])].append(splitLine)#splitLine[:4]+splitLine[6:])


            else: # Otherwise, use output from previous round that we didn't write yet
                for line in myOutString.split('\n'):
                    splitLine = line.split('\t')
                    if len(splitLine) > 1:
                        lineCounter += 1
                        if not str(splitLine[0]) in recombToLines:
                            recombToLines[str(splitLine[0])] = []
                        recombToLines[str(splitLine[0])].append(splitLine)
            currentLen = lineCounter
        
            keyToCombinedLines = {}
            myOutString = ''
            for r in recombToLines.keys():
                keyToIndices = {}
                indicesToKeys = {}
                for i in range(0,len(recombToLines[r])):
                    justCombined = False # reset for new i
                    for j in range(i,len(recombToLines[r])):
                        myKey = str(r)+'_'+str(i)+'_'+str(j)
                        l1 = (recombToLines[r])[i]
                        keyToCombinedLines[str(r)+'_'+str(i)] = l1
                        keyToIndices[str(r)+'_'+str(i)] = [i]
                        l2 = (recombToLines[r])[j]
                        keyToCombinedLines[str(r)+'_'+str(j)] = l2
                        keyToIndices[str(r)+'_'+str(j)] = [j]

                        if justCombined == True: # if we just combined i with the previous j, use that last recombinant
                            l1 = keyToCombinedLines[prevKey]
                            myKey = prevKey+'_'+str(j)
                        
                        if l1[3:] == l2[3:]: # if same parents
                            ## combine if: one interval is identical and other overlaps in any way
                            bp1a = toInt(l1[1][1:-1].split(','))
                            bp1b = toInt(l1[2][1:-1].split(','))
                            bp2a = toInt(l2[1][1:-1].split(','))
                            bp2b = toInt(l2[2][1:-1].split(','))
                            if (bp2a[0] >= bp1a[0] and bp2a[0] <= bp1a[1]+1) and (bp2b[0] >= bp1b[0] and bp2b[0] <= bp1b[1]+1):
                                newBP1 = '('+str(min(bp1a[0],bp2a[0]))+','+str(max(bp1a[1],bp2a[1]))+')'
                                newBP2 = '('+str(min(bp1b[0],bp2b[0]))+','+str(max(bp1b[1],bp2b[1]))+')'
                                keyToCombinedLines[myKey] = l1[:1]+[newBP1,newBP2]+l1[3:]
                                keyToIndices[myKey] = []
                                for ind in myKey.split('_')[1:]:
                                    keyToIndices[myKey].append(int(ind))
                                    if not ind in indicesToKeys:
                                        indicesToKeys[ind] = []
                                    indicesToKeys[ind].append(myKey)
                                justCombined = True
                                prevKey = myKey
                            else:
                                justCombined = False

                alreadyPrinted = {}
                myKeepKeys = {}
                for k in sorted(keyToIndices, key=lambda k: len(keyToIndices[k]), reverse=True):
                    keep = True
                    for ind in keyToIndices[k]:
                        if ind in alreadyPrinted:
                            keep = False
                    if keep == True:
                        myKeepKeys[k] = [str(keyToCombinedLines[k][3]), str(keyToCombinedLines[k][6]), int(keyToCombinedLines[k][1][1:-1].split(',')[0]), int(keyToCombinedLines[k][2][1:-1].split(',')[0])]
                        for ind in keyToIndices[k]:
                            alreadyPrinted[ind] = True
                for k in sorted(myKeepKeys, key=lambda k: (myKeepKeys[k][0], myKeepKeys[k][1], myKeepKeys[k][2], myKeepKeys[k][3])):
                    #print(k, keyToCombinedLines[k])
                    myOutString += joiner(keyToCombinedLines[k])+'\n'
            print(myOutString.count('\n'), currentLen)
            if myOutString.count('\n') == currentLen:
                sys.stderr.write('Converged on final output. Printing combined file.\n')
                open('filtering/data/combinedCatOnlyBest.txt','w').write(myOutString)
                myFlag = 'printed'

    if myFlag == '':
        sys.stderr.write('Did not converge after 10 iterations of combining. Printing combined file.\n')
        open('filtering/data/combinedCatOnlyBest.txt','w').write(myOutString)
    nodeToDesc = {}
    with open('filtering/data/descendants.tsv') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                if 'node_' in splitLine[0]:
                    splitLine[0] = splitLine[0].replace('node_', '')
                nodeToDesc[int(splitLine[0])] = splitLine[1]

    russNull = {}
    russOrigParsToTotal = {}
    with open('filtering/russ_null.txt') as f:
        for line in f:
            splitLine = (line.strip()).split()
            if len(splitLine) == 1 and splitLine[0].isdigit():
                myOrigPars = int(splitLine[0])
                russNull[myOrigPars] = {}
                russOrigParsToTotal[myOrigPars] = 0
            elif len(splitLine) == 2:
                (russNull[myOrigPars])[int(splitLine[0])] = int(splitLine[1])
                russOrigParsToTotal[myOrigPars] += int(splitLine[1])

    robNull = {}
    robOrigParsToTotal = {}
    with open('filtering/rob_null.txt') as f:
        for line in f:
            splitLine = (line.strip()).split()
            print(splitLine)
            if len(splitLine) == 1 and splitLine[0].isdigit():
                myOrigPars = int(splitLine[0])
                robNull[myOrigPars] = {}
                robOrigParsToTotal[myOrigPars] = 0
            elif len(splitLine) == 2:
                (robNull[myOrigPars])[int(splitLine[0])] = int(splitLine[1])
                robOrigParsToTotal[myOrigPars] += int(splitLine[1])

    print(robOrigParsToTotal)
    print(robNull)

    myOutString = '#recomb_node_id\tbreakpoint-1_interval\tbreakpoint-2_interval\tdonor_node_id\tdonor_is_sibling\tdonor_parsimony\t'
    myOutString += 'acceptor_node_id\tacceptor_is_sibling\tacceptor_parsimony\toriginal_parsimony\tmin_starting_parsimony\trecomb_parsimony'
    myOutString += '\trob_pval\truss_pval\tdescendants\n'
    with open('filtering/data/combinedCatOnlyBest.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if int(splitLine[-2]) > 0 and (int(splitLine[-2])-int(splitLine[-1])) >= 3:

                if int(splitLine[-2]) not in robNull:
                    myRobNull = 'NA'
                else:
                    myTotal = robOrigParsToTotal[int(splitLine[-2])]
                    myImprovement = int(splitLine[-2])-int(splitLine[-1])
                    for k in sorted(robNull[(int(splitLine[-2]))].keys()):
                        if k < myImprovement:
                            myTotal -= robNull[(int(splitLine[-2]))][k]
                    if myTotal == 0:
                        myRobNull = '0/'+str(robOrigParsToTotal[int(splitLine[-2])])
                    else:
                        myRobNull = float(myTotal)/float(robOrigParsToTotal[int(splitLine[-2])])


                if int(splitLine[-2]) not in russNull:
                    myRussNull = 'NA'
                else:
                    myTotal = russOrigParsToTotal[int(splitLine[-2])]
                    myImprovement = int(splitLine[-2])-int(splitLine[-1])
                    for k in sorted(russNull[(int(splitLine[-2]))].keys()):
                        if k < myImprovement:
                            myTotal -= russNull[(int(splitLine[-2]))][k]
                    if myTotal == 0:
                        myRussNull = '0/'+str(russOrigParsToTotal[int(splitLine[-2])])
                    else:
                        myRussNull = float(myTotal)/float(russOrigParsToTotal[int(splitLine[-2])])
                splitLine.append(myRobNull)
                splitLine.append(myRussNull)
                splitLine.append(nodeToDesc[int(splitLine[0])])
                myOutString += joiner(splitLine)+'\n'
    print("Writing with Pvals")
    open('filtering/data/combinedCatOnlyBestWithPVals.txt','w').write(myOutString)



##########################
#### HELPER FUNCTIONS ####
##########################

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

def toInt(myList):
    myReturn = []
    for k in myList:
        myReturn.append(int(k))
    return(myReturn)

def joinerU(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return '_'.join(newList)

def joinerC(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return ','.join(newList)

#########################
##### FUNCTION CALL #####
#########################

def main():
    catOnlyBest()


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit
