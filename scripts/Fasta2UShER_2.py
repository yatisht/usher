import glob
import os 
import argparse 
import re
from sequenceAnalyzer import FastAreader
import time



parser = argparse.ArgumentParser(description='Generates merged VCF that ignores indels and problematic sites and recognizes missing data in addition to genotypes.')

print('\n\nFor more information on problematic sites see:\n\nNicola De Maio, Landen Gozashti, Yatish Turakhia, Conor Walker, Robert Lanfear, Russell Corbett-Detig, and Nick Goldman, Issues with SARS-Cov-2 sequencing data: Updated analysis with data from 12th June 2020, Virological post 2020. https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/12 \n    and \nYatish Turakhia, Bryan Thornlow, Landen Gozashti, Angie S. Hinrichs, Jason D. Fernandes, David Haussler, and Russell Corbett-Detig, "Stability of SARS-CoV-2 Phylogenies", bioRxiv pre-print 2020.\n\n')

parser.add_argument('-inpath', nargs='?', required=True,
                     help='Path to a directory containing all input sequences you wish to consider')
parser.add_argument('-unaligned', nargs='?',const='T', default=None,
                    help='Signifies that user provided fasta files have not been aligned to the reference. If specified, fasta2vcf uses mafft to perform multiple sequence alignments.')
parser.add_argument('-auto_mask', nargs='?',const='T', default=None,
                    help='Ignore problematic sites per our masking recomendations')
parser.add_argument('-user_specified_mask', nargs='?', default=None,
                    help='Path to VCF fle containing custom masking recomendations (please ensure VCF format is consistent with https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf)')
parser.add_argument('-reference', nargs='?', required=True,
                     help='Path to reference fasta file. The reference sequence header should be identical to the reference sequence header in the input msa file')
parser.add_argument('-output', nargs='?', required=True,
                     help='Name of output vcf file')
parser.add_argument('-thread', nargs='?', required=False, default=1,
                     help='Number of threads for msa')

varDic = {}
gDic = {}
msaList = []

args = vars(parser.parse_args())


if args['unaligned'] != None:
    
    os.system('rm -r temp.fa')
    os.system('cat {0} > temp.fa'.format(' '.join(glob.glob('{0}/*'.format(args['inpath'])))))
    os.system('mafft --thread {1} --auto --keeplength --addfragments temp.fa {0} > inputMsa.fa'.format(args['reference'],args['thread']))
    msaName = 'inputMsa.fa'
    

else:
    msaName = ''.join(glob.glob('{0}/*'.format(args['inpath'])))
seqDic = {}
altCount = {}
refCount = {}


"""

Retreave problematic sites

"""
with open(args['reference'],'r') as f:
    for line in f:
        head = line.strip('>').strip()
        print(head)
        break


probDic = {}
if args['auto_mask'] != None:

    os.system('rm -r problematic_sites_sarsCov2.vcf')
    os.system('wget https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf')
    os.system('./faToVcf -maskSites=problematic_sites_sarsCov2.vcf -ref=\"{0}\" {1} {2}.vcf'.format(head,msaName,args['output'].strip()))




elif args['user_specified_mask'] != None:
     

    os.system('./faToVcf -maskSites={0} -ref=\"{1}\" {2} {3}.vcf'.format(args['user_specified_mask'],head,msaName,args['output']))

