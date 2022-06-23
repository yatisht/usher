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
    """Handles the input arguments from the command line. Manages 
    the argument parser."""

    def __init__(self, inOpts=None):
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line input using argparse.
        '''
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("-b", "--breakpoints", help="Number of breakpoints that each recombinant sample will have. Must be [1..4] (Default = 1).", default=1, type=int)
        self.parser.add_argument("-s", "--samples", help="Number of recombinant samples to create (Default = 100).", default=100, type=int)
        self.parser.add_argument("-c", "--copies", help="Number of identical copies to make for each recombinant sample (Default = 10).", default=10, type=int)
        self.parser.add_argument("-m", "--commonMutations", help="Number of mutations to add to each copy, shared by all in a set. (Default = 0).", default=0, type=int)
        self.parser.add_argument("-M", "--randomMutations", help="Number of mutations to add to each copy, randomly chosen for each copy. (Default = 0).", default=0, type=int)
        self.parser.add_argument("-t", "--threshold", help="Minimum mutational distance for acceptor/donor samples (Default = 10).", default=10, type=int)
        self.parser.add_argument("-f", "--fasta", help="Fasta file containing sequences for acceptor/donor samples. [Must include either -f or -d!]", default='')
        self.parser.add_argument("-d", "--differences", help="File containing all mutations relative to reference for each sample (allNodeToMuts.py output by makeMutsFile.py) [Must include either -f or -d!]", default='')
        self.parser.add_argument("-r", "--ref", help="Fasta file containing reference genome for use in creating VCF. (Default = 'wuhan.ref.fa').", default='wuhan.ref.fa')
        self.parser.add_argument("-S", "--separate", help="If enabled, will produce one MSA as a .fasta file for each set of recombinants to the argument directory. If not enabled, will not produce these files.", default=False)
        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))
        self.args = self.parser.parse_args()     
        if self.args.breakpoints < 1 or self.args.breakpoints > 4:
            sys.stderr.write("Please retry with between 1 and 4 breakpoints.\n")
            sys.exit(1)
        if self.args.fasta == '' and self.args.differences == '':
            sys.stderr.write("Please supply either a MSA via --fasta or a diff file via --differences.\n")
            sys.exit(1)

##########################
##### MAIN FUNCTIONS #####
##########################

def makeExamples(myS, myB, myC, myD, myF, myT, mym, myM, myR, mySep):

    ### Read in reference sequence
    posToRef = {}
    with open(myR) as f:
        for line in f:
            l = line.strip()
            if not l.startswith('>'):
                myReference = l.upper()
                for i in range(0,len(myReference)):
                    posToRef[i] = myReference[i]
            else:
                myRefName = l[1:]
    sys.stderr.write("Finished reading in reference.\n")


    ### Read in either differences file or MSA to get our node sequences
    nodeToDiffs = {}
    if myD != '':
        with open(myD) as f:
            for line in f:
                splitLine = (line.strip()).split('\t')
                if len(splitLine) == 1:
                    nodeToDiffs[splitLine[0]] = []
                else:
                    nodeToDiffs[splitLine[0]] = splitLine[1].split(',')

    sampleToSeq = {}
    if myF != '':
        with open(myF) as f:
            for line in f:
                l = line.strip()
                if l.startswith('>'):
                    mySample = l[1:]
                else:
                    sampleToSeq[mySample] = l
    sys.stderr.write("Finished reading in Diff/MSA input file.\n")

    recSampleToLog = {}
    recSampleToSeq = {}
    recSampleToDiffBetweenBps = {}
    while len(recSampleToSeq.keys()) < myS:

        ### Get samples for 2 sequences that will make up our recombinant
        if myD == '':
            samples = numpy.random.choice(list(sampleToSeq.keys()), size=2, replace=False)
            mySampleName = 'RECOMB_'+str(myB)+'_'+str(len(recSampleToSeq))+'_'+str(samples[0])+'_'+str(samples[1])
            s1 = samples[0]
            s2 = samples[1]
        elif myF == '':
            samples = numpy.random.choice(list(nodeToDiffs.keys()), size=2, replace=False)
            mySampleName = 'RECOMB_'+str(myB)+'_'+str(len(recSampleToSeq))+'_'+str(samples[0])+'_'+str(samples[1])
            s1 = samples[0]
            sampleToSeq[s1] = addMuts(myReference, nodeToDiffs[s1])
            s2 = samples[1]
            sampleToSeq[s2] = addMuts(myReference, nodeToDiffs[s2])

        ### Create recombinant sequence from our two samples
        myTotalDiff = getDiff(sampleToSeq[s1],sampleToSeq[s2], 0)
        if myB == 1:
            bps = numpy.random.choice(sorted(list(posToRef.keys())),size=1, replace=False)
            bp1 = bps[0]
            mySeq = sampleToSeq[s1][:bp1]+sampleToSeq[s2][bp1:]
            myDiff = minLen(getDiff(sampleToSeq[s1][:bp1], sampleToSeq[s2][:bp1], 0), getDiff(sampleToSeq[s1][bp1:], sampleToSeq[s2][bp1:], bp1))
            print(myDiff, myT)
        elif myB == 2:
            bps = sorted(numpy.random.choice(sorted(list(posToRef.keys())),size=2, replace=False))
            while bps[1]-bps[0] <= 1000:
                bps = sorted(numpy.random.choice(sorted(list(posToRef.keys())),size=2, replace=False))
            bp1 = bps[0]
            bp2 = bps[1]
            mySeq = sampleToSeq[s1][:bp1]+sampleToSeq[s2][bp1:bp2]+sampleToSeq[s1][bp2:]
            myDiff =  minLen(getDiff(sampleToSeq[s1][bp1:bp2], sampleToSeq[s2][bp1:bp2], bp1), getDiff(sampleToSeq[s1][:bp1]+sampleToSeq[s1][bp2:], sampleToSeq[s2][:bp1]+sampleToSeq[s2][bp2:], 0))
        elif myB == 3:
            bps = sorted(numpy.random.choice(sorted(list(posToRef.keys())),size=3, replace=False))
            while bps[1]-bps[0] <= 1000 or bps[2]-bps[1] <= 1000:
                bps = sorted(numpy.random.choice(sorted(list(posToRef.keys())),size=3, replace=False))
            bp1 = bps[0]
            bp2 = bps[1]
            bp3 = bps[2]
            mySeq = sampleToSeq[s1][:bp1]+sampleToSeq[s2][bp1:bp2]+sampleToSeq[s1][bp2:bp3]+sampleToSeq[s2][bp3:]
            diff1 = getDiff(sampleToSeq[s1][bp1:bp2]+sampleToSeq[s1][bp3:], sampleToSeq[s2][bp1:bp2]+sampleToSeq[s2][bp3:], bp1)
            diff2 = getDiff(sampleToSeq[s1][:bp1]+sampleToSeq[s1][bp2:bp3], sampleToSeq[s2][:bp1]+sampleToSeq[s2][bp2:bp3], 0)
            myDiff =  minLen(diff1, diff2)
        elif myB == 4:
            bps = sorted(numpy.random.choice(sorted(list(posToRef.keys())),size=4, replace=False))
            while bps[1]-bps[0] <= 1000 or bps[2]-bps[1] <= 1000 or bps[3]-bps[2] <= 1000:
                bps = sorted(numpy.random.choice(sorted(list(posToRef.keys())),size=4, replace=False))
            bp1 = bps[0]
            bp2 = bps[1]
            bp3 = bps[2]
            bp4 = bps[3]
            mySeq = sampleToSeq[s1][:bp1]+sampleToSeq[s2][bp1:bp2]+sampleToSeq[s1][bp2:bp3]+sampleToSeq[s2][bp3:bp4]+sampleToSeq[s1][bp4:]
            diff1 = getDiff(sampleToSeq[s1][bp1:bp2]+sampleToSeq[s1][bp3:bp4], sampleToSeq[s2][bp1:bp2]+sampleToSeq[s2][bp3:bp4], bp1)
            diff2 = getDiff(sampleToSeq[s1][:bp1]+sampleToSeq[s1][bp2:bp3]+sampleToSeq[s1][bp4:], sampleToSeq[s2][:bp1]+sampleToSeq[s2][bp2:bp3]+sampleToSeq[s2][bp4:], 0)
            myDiff =  minLen(diff1, diff2)

        ### If the differences between the parents of our recombinant is above the threshold between all adjacent breakpoint pairs, keep it
        if len(myDiff) >= myT:
            recSampleToLog[mySampleName] = [joiner(samples), joiner(bps), len(myTotalDiff), len(myDiff)]
            if mym > 0: ### Add common mutations here, prior to making copies
                myMuts = []
                for m in range(0, mym):
                    mySeq, myMut = addMut(mySeq, numpy.random.choice(sorted(list(posToRef.keys())),size=1, replace=False)[0])
                    myMuts.append(myMut)
                recSampleToLog[mySampleName].append(joiner(myMuts))
            recSampleToSeq[mySampleName] = mySeq
            recSampleToDiffBetweenBps[mySampleName] = myDiff
            sys.stderr.write("Generated "+str(len(recSampleToSeq.keys()))+" recombinant sequences.\n")
    sys.stderr.write("Finished generating recombinant sequences.\n")

    bpToHeader = {}
    bpToHeader[1] = 'recombinant_sample\tparent1\tparent2\tgenetic_distance\tmutations_in_recomb_tract\tbp1\n'
    bpToHeader[2] = 'recombinant_sample\tparent1\tparent2\tgenetic_distance\tmutations_in_recomb_tract\tbp1\tbp2\n'
    bpToHeader[3] = 'recombinant_sample\tparent1\tparent2\tgenetic_distance\tmutations_in_recomb_tract\tbp1\tbp2\tbp3\n'
    bpToHeader[4] = 'recombinant_sample\tparent1\tparent2\tgenetic_distance\tmutations_in_recomb_tract\tbp1\tbp2\tbp3\tbp4\n'

    ### Write our MSA, which we will use to create a VCF to add to the starting .pb via faToVcf and UShER
    myOutMSA = '>'+myRefName+'\n'+myReference+'\n'
    mySepMSAs = []
    myOutLog = bpToHeader[myB]
    myOutDiff = ''
    for s in recSampleToSeq:
        tempMSA = '>'+myRefName+'\n'+myReference+'\n' ### Fasta should have reference sequence first, for use with faToVcf
        for x in range(0,myC):
            if myM > 0:
                myMuts = []
                mySeq = recSampleToSeq[s]
                for m in range(0, myM):
                    mySeq, myMut = addMut(mySeq, numpy.random.choice(sorted(list(posToRef.keys())),size=1, replace=False)[0])
                    myMuts.append(myMut)
                myOutMSA += '>'+s+'_X'+str(x)+'\n'+mySeq+'\n'
                tempMSA += '>'+s+'_X'+str(x)+'\n'+mySeq+'\n'
                myOutLog += s+'_X'+str(x)+'\t'+doubleJoiner(recSampleToLog[s])+'\t'+joiner(myMuts)+'\n'
                myOutDiff += s+'_X'+str(x)+'\t'+joinerC(recSampleToDiffBetweenBps[s])+'\n'
            else:
                myOutMSA += '>'+s+'_X'+str(x)+'\n'+recSampleToSeq[s]+'\n'
                tempMSA += '>'+s+'_X'+str(x)+'\n'+recSampleToSeq[s]+'\n'
                myOutLog += s+'_X'+str(x)+'\t'+doubleJoiner(recSampleToLog[s])+'\n'
                myOutDiff += s+'_X'+str(x)+'\t'+joinerC(recSampleToDiffBetweenBps[s])+'\n'
        mySepMSAs.append(tempMSA)
    open('recombination_'+str(myB)+'_'+str(myC)+'_'+str(myM)+'.msa.fa','w').write(myOutMSA)
    open('recombination_'+str(myB)+'_'+str(myC)+'_'+str(myM)+'.log','w').write(myOutLog)
    open('recombination_'+str(myB)+'_'+str(myC)+'_'+str(myM)+'.differences.txt','w').write(myOutDiff)

    myOutFaToVcf = ''
    if mySep != False:
        if not os.path.exists(mySep):
            os.mkdir('./'+mySep)
        for i in range(1,len(mySepMSAs)+1):
            open('./'+mySep+'/recombinant_set_'+str(i)+'.fa','w').write(mySepMSAs[i-1])
            myOutFaToVcf += 'faToVcf recombinant_set_'+str(i)+'.fa recombinant_set_'+str(i)+'.vcf\n'
        open('./'+mySep+'/makeSetVCFs.sh','w').write(myOutFaToVcf)

##########################
#### HELPER FUNCTIONS ####
##########################

def addMut(seq, pos):
    #print(pos, len(seq))
    myReturn = []
    for i in range(0,len(seq)):
        if i != int(pos):
            myReturn.append(seq[i])
        else:
            if seq[i] == 'A':
                mySub = numpy.random.choice(['C','G','T'], size=1)[0]
            elif seq[i] == 'C':
                mySub = numpy.random.choice(['A','G','T'], size=1)[0]
            elif seq[i] == 'G':
                mySub = numpy.random.choice(['C','A','T'], size=1)[0]
            elif seq[i] == 'T':
                mySub = numpy.random.choice(['A','G','C'], size=1)[0]
            else:
                mySub = numpy.random.choice(['A','G','C','T'], size=1)[0]
            myMut = seq[i]+str(pos)+mySub
            myReturn.append(mySub)
    return(''.join(myReturn), myMut)

def addMuts(ref, muts):
    myReturn = []
    for k in list(ref):
        myReturn.append(k)
    if len(muts) > 0:
        for k in muts:
            myCoord = int(k[1:-1])-1
            if not myReturn[myCoord] == k[0]:
                print(k, myReturn[myCoord])
            myReturn[myCoord] = str(k[-1])
            if not myReturn[myCoord] == k[-1]:
                print(k, myReturn[myCoord])
    return(''.join(myReturn))

def getDiff(s1, s2, add):
    myReturn = []
    for i in range(0,len(s1)):
        if s1[i] != s2[i]:
            myReturn.append(add+i)
    return(myReturn)

def minLen(l1, l2):
    if len(l1) <= len(l2):
        return(l1)
    else:
        return(l2)

def doubleJoiner(myList):
    myReturn = []
    for k in myList:
        if type(k) == list:
            myReturn.append(joiner(k))
        else:
            myReturn.append(k)
    return(joiner(myReturn))

def replaceSymbols(myEntry):
    myEntry = myEntry.replace('|', '_')
    myEntry = myEntry.replace('/', '_')
    return(myEntry)

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
        newList.append(str(k))
    return(','.join(newList))

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
    else:
        myS = 100
    if myCommandLine.args.breakpoints:
        myB = myCommandLine.args.breakpoints
    else:
        myB = 1
    if myCommandLine.args.copies:
        myC = myCommandLine.args.copies
    else:
        myC = 10
    if str(myCommandLine.args.threshold):
        myT = myCommandLine.args.threshold
    else:
        myT = 10
    if myCommandLine.args.fasta:
        myF = myCommandLine.args.fasta
    else:
        myF = ''
    if myCommandLine.args.commonMutations:
        mym = myCommandLine.args.commonMutations
    else:
        mym = 0
    if myCommandLine.args.randomMutations:
        myM = myCommandLine.args.randomMutations
    else:
        myM = 0
    if myCommandLine.args.ref:
        myR = myCommandLine.args.ref
    else:
        myR = 'wuhan.ref.fa'
    if myCommandLine.args.separate:
        mySep = myCommandLine.args.separate
    else:
        mySep = False
    if myCommandLine.args.differences:
        myD = myCommandLine.args.differences
    else:
        myD = ''

    makeExamples(myS, myB, myC, myD, myF, myT, mym, myM, myR, mySep)

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit





















