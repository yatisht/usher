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

def makeSampleFiles():
    myParallelJobs = ''
    bmToDir = {'10':'S1_500/','11':'S1_500_1/','12':'S1_500_2/','13':'S1_500_3/','20':'S2_500/','21':'S2_500_1/','22':'S2_500_2/','23':'S2_500_3/','30':'S3_500/','31':'S3_500_1/','32':'S3_500_2/','33':'S3_500_3/','40':'S4_500/','41':'S4_500_1/','42':'S4_500_2/','43':'S4_500_3/'}
    for b in ['1','2','3','4']:
        for m in ['0','1','2','3']:
            l = 0
            if os.path.exists('recombination_'+b+'_1_'+m+'.log'):
                with open('recombination_'+b+'_1_'+m+'.log') as f:
                    for line in f:
                        splitLine = (line.strip()).split('\t')
                        if not splitLine[0] == 'recombinant_sample':
                            l += 1
                            # print(splitLine[0])
                            open('TEMP_SAMPLES/temp_sample_'+b+m+'_'+str(l)+'.txt','w').write(splitLine[0]+'\n')
                            if not os.path.exists('SIM_DATA/recombination_'+b+m+'_'+str(l)+'.tsv'):
                                myBashSh = '/HOME/usher/build/usher -v ../'+bmToDir[b+m]+'recombinant_set_'+str(l)+'.vcf.gz -i ../public-latest.all.masked.pb -o ../temp_'+b+m+'_'+str(l)+'.pb -T 32 2> ../LOGS/temp_'+b+m+'_'+str(l)+'.log\n'
                                myBashSh += 'grep Parsimony ../LOGS/temp_'+b+m+'_'+str(l)+'.log > ../SIM_DATA/temp_'+b+m+'_'+str(l)+'.parsimony\nmkdir ../'+b+m+'_'+str(l)+'/\n'
                                myBashSh += '/usr/bin/time -o '+b+'_'+m+'_'+str(l)+'.time -f "%E %M" /HOME/usher/build/ripples -i ../temp_'+b+m+'_'+str(l)+'.pb -s ../TEMP_SAMPLES/temp_sample_'+b+m+'_'+str(l)+'.txt -d ../'+b+m+'_'+str(l)+'/ -T 32 2> ../LOGS/findRec_'+b+m+'_'+str(l)+'.log\n'
                                myBashSh += 'mv ../'+b+m+'_'+str(l)+'/recombination.tsv ../SIM_DATA/recombination_'+b+m+'_'+str(l)+'.tsv\n'
                                myBashSh += 'mv ../'+b+m+'_'+str(l)+'/descendants.tsv ../SIM_DATA/descendants_'+b+m+'_'+str(l)+'.tsv\nrm ../temp_'+b+m+'_'+str(l)+'.pb\n'
                                open('SIM_SCRIPTS/'+b+m+'_'+str(l)+'.sh','w').write(myBashSh)
                                myParallelJobs += 'bash '+b+m+'_'+str(l)+'.sh\n'
    open('SIM_SCRIPTS/myParallelJobs.sh','w').write(myParallelJobs)



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
    makeSampleFiles()


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit