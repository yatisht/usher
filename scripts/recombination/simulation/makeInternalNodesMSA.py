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
        self.parser.add_argument("-s", "--samples", help="Sample mutation paths file. To generate: matUtils extract -i <input.pb> -A <sample-paths.tsv>. [REQUIRED]", default='')
        self.parser.add_argument("-l", "--leaves", help="File containing number of leaves for each node. Only nodes with at least 10 descendants will be used, "+
            "as this is required for detection. To generate: matUtils extract -i <input.pb -L <num-leaves.tsv>. [REQUIRED]", default='')
        self.parser.add_argument("-t", "--threshold", help="By default, this script includes only nodes with at least 10 descendant leaves. To use a different minimum, specify here.",default=10)
        self.parser.add_argument("-r", "--ref", help="Fasta file containing reference genome. (Default = 'wuhan.ref.fa').", default='wuhan.ref.fa')
        
        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))
        self.args = self.parser.parse_args()     

##########################
##### MAIN FUNCTIONS #####
##########################

def getMutationsFile(myS, myL, myT, myR):
    myReference = ''
    with open(myR) as f:
        for line in f:
            l = line.strip()
            if not l.startswith('>'):
                myReference += l.strip().upper()

    nodeToMuts = {1:[]}
    with open(myS) as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            myPath = stripEach(splitLine[1].split('>'))
            pastMuts = []
            for k in myPath:
                if k.endswith(')'): # we don't care about mutations on sample leafs
                    s = k.split()
                    if len(s) == 1:
                        myNode = int((s[0])[1:-1])
                        mutSplit = []
                    else:
                        myNode = int((s[1])[1:-1])
                        mutSplit = s[0].split(',')
                    if not myNode in nodeToMuts:
                        nodeToMuts[myNode] = []
                        for m in pastMuts:
                            nodeToMuts[myNode].append(m)
                        for m in mutSplit:
                            nodeToMuts[myNode].append(m)
                    for m in mutSplit:
                        pastMuts.append(m)


    goodNodes = {}
    with open(myL) as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if splitLine[0].isdigit() and int(splitLine[1]) >= myT:
                goodNodes[int(splitLine[0])] = True

    myOutMSA = ''
    myOutString = ''
    for n in nodeToMuts:
        if n in goodNodes and n != 1: # don't use the root -- it has no mutations and is not useful for our purposes
            myOutMSA += '>NODE'+str(n)+'\n'+makeChanges(myReference, nodeToMuts[n])+'\n'
            myOutString += str(n)+'\t'+','.join(nodeToMuts[n])+'\n'
    open('internal_nodes.msa.fa','w').write(myOutMSA)
    open('nodeToMuts.txt','w').write(myOutString)


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
    else:
        myT = 10
    if myCommandLine.args.ref:
        myR = myCommandLine.args.ref

    getMutationsFile(myS, myL, myT, myR)


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit

























