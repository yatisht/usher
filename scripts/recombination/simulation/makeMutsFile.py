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
import argparse

"""
Note: Newest matUtils extract -S does not always output all nodes regardless of whether they house a mutation or not.
Using commit a6f65ade7a6606ef75902ee290585c6db23aeba6 for usher
"""

##########################
###### COMMAND LINE ######
##########################

class CommandLine(object):
    """
    Handles the input arguments from the command line. Manages 
    the argument parser.
    """

    def __init__(self, inOpts=None):
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line input using argparse.
        '''
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("-s", "--samples", help="Sample mutation paths file. To generate: matUtils extract -i <input.pb> -S <sample-paths.tsv>. [REQUIRED]", default='')
        self.parser.add_argument("-l", "--leaves", help="File containing number of leaves for each node. Only nodes with at least 10 descendants will be used, "+
            "as this is required for detection. To generate: matUtils extract -i <input.pb> -L <num-leaves.tsv>. Required if --threshold > 0.", default='')
        self.parser.add_argument("-t", "--threshold", help="By default, this script includes only nodes with at least 10 descendant leaves. To use a different minimum, specify here.",default=10)
        self.parser.add_argument("-r", "--ref", help="Fasta file containing reference genome. (Default = 'wuhan.ref.fa').", default='wuhan.ref.fa')
        self.parser.add_argument("-f", "--fasta", help="Write fasta file containing all samples according to their mutation paths. This makes the script significantly slower. (Default = False)",default=False)
        
        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))
        self.args = self.parser.parse_args()     

##########################
##### MAIN FUNCTIONS #####
##########################

def getMutationsFile(myS,myL,myT,myR,myF):
    myReference = ''
    with open(myR) as f:
        for line in f:
            l = line.strip()
            if not l.startswith('>'):
                myReference += l.strip().upper()
    sys.stderr.write("Finished reading in reference.\n")

    nodeToMuts = {1:[]}
    nodeToPars = {}
    lineCounter = 0
    with open(myS) as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            myPath = stripEach(splitLine[1].split('>'))
            pastMuts = []
            for k in myPath:
                if k.endswith(')'): # for internal nodes
                    s = k.split()
                    if len(s) == 1:
                        myNode = str((s[0])[1:-1])
                        mutSplit = []
                    else:
                        myNode = str((s[1])[1:-1])
                        mutSplit = s[0].split(',')
                    nodeToPars[myNode] = len(mutSplit)
                    if not myNode in nodeToMuts:
                        nodeToMuts[myNode] = []
                        for m in pastMuts:
                            nodeToMuts[myNode].append(m)
                        for m in mutSplit:
                            nodeToMuts[myNode].append(m)
                    for m in mutSplit:
                        pastMuts.append(m)
                else: # this is the sample leaf
                    mutSplit = k.split(',')
                    nodeToMuts[str(splitLine[0])] = []
                    nodeToPars[splitLine[0]] = len(mutSplit)
                    for m in pastMuts:
                        nodeToMuts[splitLine[0]].append(m)
                    for m in mutSplit:
                        nodeToMuts[splitLine[0]].append(m)
    sys.stderr.write("Finished storing all mutations.\n")

    myOutMSA = ''
    myOutString = ''
    myOutPars = ''
    if myT > 0: # If threshold specified for leaves, read leaves file and only output nodes passing threshold
        with open(myL) as f:
            for line in f:
                splitLine = (line.strip()).split('\t')
                if splitLine[0].isdigit() and int(splitLine[1]) >= myT and int(splitLine[0]) != 1:
                    if myF != False:
                        myOutMSA += '>'+str(splitLine[0])+'\n'+makeChanges(myReference, nodeToMuts[str(splitLine[0])])+'\n'
                    myOutString += str(splitLine[0])+'\t'+','.join(nodeToMuts[str(splitLine[0])])+'\n'
                    myOutPars += str(splitLine[0])+'\t'+str(nodeToPars[str(splitLine[0])])+'\n'

    else: # Otherwise, output all nodes
        for n in nodeToMuts:
            if n != 1:
                if myF:
                    myOutMSA += '>'+str(n)+'\n'+makeChanges(myReference, nodeToMuts[n])+'\n'
                myOutString += str(n)+'\t'+','.join(nodeToMuts[n])+'\n'
                myOutPars += str(n)+'\t'+str(nodeToPars[n])+'\n'
    if myF:
        open('allNodesT'+str(myT)+'.msa.fa','w').write(myOutMSA)
    open('allNodeToMutsT'+str(myT)+'.txt','w').write(myOutString)
    open('allNodeToParsT'+str(myT)+'.txt','w').write(myOutPars)
    sys.stderr.write("Finished writing all output files.\n")


##########################
#### HELPER FUNCTIONS ####
##########################

def stripEach(myList):
    myReturn = []
    for k in myList:
        myReturn.append(k.strip())
    return(myReturn)

def makeChanges(ref, muts):
    myReturn = []
    for k in list(ref):
        myReturn.append(k)
    for k in muts:
        myCoord = int(k[1:-1])-1
        if not myReturn[myCoord] == k[0]:
            print(k, myReturn[myCoord])
        myReturn[myCoord] = str(k[-1])
        if not myReturn[myCoord] == k[-1]:
            print(k, myReturn[myCoord])
    return(''.join(myReturn))

def getPos(myInds, intLineNumToPos):
    myReturn = []
    for k in myInds:
        myReturn.append(intLineNumToPos[k])
    return(myReturn)

def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return '\t'.join(newList)

def joinerU(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return '_'.join(newList)

def joinerC(entry):
    newList = []
    for k in entry:
        newList.append('"'+str(k)+'"')
    return ','.join(newList)

#########################
##### FUNCTION CALL #####
#########################

def main(myCommandLine=None):
    """
    Initializes a CommandLine object and passes the provided 
    arguments into a new fileConverter object and calls main method.
    """
    myCommandLine = CommandLine()

    # Necessary files:
    if myCommandLine.args.samples:
        myS = myCommandLine.args.samples
    if myCommandLine.args.leaves:
        myL = myCommandLine.args.leaves
    if myCommandLine.args.threshold:
        myT = int(myCommandLine.args.threshold)
        if myT == 0:
            myL = ''
    else:
        myT = 10
    if myCommandLine.args.ref:
        myR = myCommandLine.args.ref
    else:
        myR = 'wuhan.ref.fa'
    if myCommandLine.args.fasta:
        myF = myCommandLine.args.fasta
    else:
        myF = False

    getMutationsFile(myS, myL, myT, myR, myF)


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit

























