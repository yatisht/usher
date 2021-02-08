import glob
import os
import argparse
import re
from Bio import SeqIO
import tempfile
import subprocess

parser = argparse.ArgumentParser(
    description='Generates merged VCF that ignores indels \
                 and problematic sites and \
                 recognizes missing data in addition to genotypes.')

print('\n\nFor more information on problematic sites see:\n\nNicola De Maio, '                                                                                                                                  
      'Landen Gozashti, Yatish Turakhia, Conor Walker, Robert Lanfear, '
      'Russell Corbett-Detig, and Nick Goldman, Issues with SARS-Cov-2 '
      'sequencing data: Updated analysis with data from 12th June 2020, '
      'Virological post 2020. '
      'https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/12 '
      '\n    and \nYatish Turakhia, Bryan Thornlow, Landen Gozashti, Angie '
      'S. Hinrichs, Jason D. Fernandes, David Haussler, and Russell '
      'Corbett-Detig, "Stability of SARS-CoV-2 Phylogenies", bioRxiv pre-print '
      '2020.\n\n')

# Add argparse arguments
parser.add_argument(
    '-inpath',
    nargs='?',
    required=True,
    help='Path to a directory containing all input sequences you wish to consider')
parser.add_argument(
    '-unaligned',
    nargs='?',
    const='T',
    default=None,
    help='Signifies that user provided fasta files have not been \
                    aligned to the reference. If specified, faToVcf \
                    uses mafft to perform multiple sequence alignments.')
parser.add_argument(
    '-auto_mask',
    nargs='?',
    const='T',
    default=None,
    help='Ignore problematic sites per our masking recomendations')
parser.add_argument(
    '-user_specified_mask',
    nargs='?',
    default=None,
    help='Path to VCF fle containing custom masking recomendations \
                    (please ensure VCF format is consistent with \
                    https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf)')
parser.add_argument('-reference', nargs='?', required=True,
                    help='Path to reference fasta file. \
                     The reference sequence header should be \
                     identical to the reference sequence header \
                     in the input msa file')
parser.add_argument('-output', nargs='?', required=True,
                    help='Name of output vcf file')
parser.add_argument('-thread', nargs='?', required=False, default=1,
                    help='Number of threads for msa')

varDic = {}
gDic = {}
msaList = []
headDic = {}
args = vars(parser.parse_args())
filePath = '/'.join(os.path.realpath(__file__).split('/')[0:-1])
# Check that faToVcf is in the same directory as Fasta2UShER
if os.path.isfile('{0}/faToVcf'.format(filePath)) == False:
    raise Exception("faToVcf must be in the same directory as Fasta2UShER")

"""

If the user provided unaligned sequences, perform an MSA with mafft
Otherwise, read the name of the provided msa file

"""
if args['unaligned'] is not None:

    # create temporary file to harbor input sequences to MSA
    temp = tempfile.NamedTemporaryFile(mode='wt', delete=False)
    for name in glob.glob('{0}/*'.format(args['inpath'])):
        fasta_sequences = SeqIO.parse(open(name), 'fasta')
        for fasta in fasta_sequences:
            header, sequence = fasta.id, str(fasta.seq)
            if header not in headDic:
                headDic[header] = ''
                temp.write('>{0}\n{1}\n'.format(header, sequence))
    # Perform MSA
    subprocess.call(
        'mafft --thread {1} --auto --keeplength --addfragments {3}  {0} > {2}/inputMsa.fa'.format(
            args['reference'],
            args['thread'],
            filePath,
            temp.name),
        shell=True)
    msaName = '{0}/inputMsa.fa'.format(filePath)
    temp.close()

else:
    msaName = ''.join(glob.glob('{0}/*'.format(args['inpath'])))
seqDic = {}
altCount = {}
refCount = {}


"""

Read in reference

"""
with open(args['reference'], 'r') as f:
    for line in f:
        head = line.strip('>').strip()
        break


"""

Run faToVcf

"""

probDic = {}

# Retrieve recomended problematic sites if the user specied and run faToVcf

if args['auto_mask'] is not None:

    subprocess.call(
        'wget https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf',
        shell=True)  # download masking recomendations
    # if the user specified ".vcf" in the output file name, do not add ".vcf"
    # extension
    if '.vcf' in args['output']:
        subprocess.call(
            '{3}/faToVcf -maskSites={3}/problematic_sites_sarsCov2.vcf -ref=\"{0}\" {1} {2}'.format(
                head, msaName, args['output'].strip()), shell=True)
    # else add ".vcf" extension
    else:  
        subprocess.call(
            '{3}/faToVcf -maskSites={3}/problematic_sites_sarsCov2.vcf -ref=\"{0}\" {1} {2}.vcf'.format(
                head, msaName, args['output'].strip()), shell=True)


# Retrieve user provided problematic sites if the user specied and run
# faToVcf

elif args['user_specified_mask'] is not None:
    # if the user specified ".vcf" in the output file name, do not add ".vcf"
    # extension
    if '.vcf' in args['output']:
        subprocess.call(
            '{3}/faToVcf -maskSites={3}/{0} -ref=\"{1}\" {2} {3}'.format(
                args['user_specified_mask'],
                head,
                msaName,
                args['output'],
                filePath),
            shell=True)
    # else add ".vcf" extension
    else:  
        subprocess.call('{3}/faToVcf -maskSites={3}/{0} -ref=\"{1}\" {2} {3}.vcf'.format(
            args['user_specified_mask'], head, msaName, args['output'], filePath), shell=True)


# Else run faToVcf without masking

else:

    # if the user specified ".vcf" in the output file name, do not add ".vcf"
    # extention
    if '.vcf' in args['output']:
        subprocess.call('{3}/faToVcf  -ref=\"{0}\" {1} {2}'.format(head,
                                                                   msaName, args['output'], filePath), shell=True)
    # else add ".vcf" extention
    else:
        subprocess.call('{3}/faToVcf  -ref=\"{0}\" {1} {2}.vcf'.format(head,
                                                                       msaName, args['output'], filePath), shell=True)

# remove masking recomendations file if it was downloaded
if os.path.isfile('{0}/problematic_sites_sarsCov2.vcf'.format(filePath)):

    subprocess.call(
        'rm -r {0}/problematic_sites_sarsCov2.vcf'.format(filePath),
        shell=True)
