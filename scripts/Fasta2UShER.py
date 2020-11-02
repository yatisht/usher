import glob
import os
import argparse
import re
from Bio import SeqIO
import tempfile


parser = argparse.ArgumentParser(
    description='Generates merged VCF that ignores indels and problematic sites and \
                                              recognizes missing data in addition to genotypes.')

print('\n\nFor more information on problematic sites see:\n\nNicola De Maio, Landen Gozashti, Yatish Turakhia, Conor Walker, Robert Lanfear, Russell Corbett-Detig, and Nick Goldman, Issues with SARS-Cov-2 sequencing data: Updated analysis with data from 12th June 2020, Virological post 2020. https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/12 \n    and \nYatish Turakhia, Bryan Thornlow, Landen Gozashti, Angie S. Hinrichs, Jason D. Fernandes, David Haussler, and Russell Corbett-Detig, "Stability of SARS-CoV-2 Phylogenies", bioRxiv pre-print 2020.\n\n')

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
    help='Signifies that user provided fasta files have not been aligned to the reference. \
                    If specified, fasta2vcf uses mafft to perform multiple sequence alignments.')
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
                     The reference sequence header should be identical to the reference sequence header in the input msa file')
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
# Check that fastaToVcf is in the same directory as Fasta2UShER
if os.path.isfile('{0}/faToVcf'.format(filePath)) == False:
    raise Exception("faToVcf must be in the same directory as Fasta2UShER")

# Change permissions to faToVcf
os.system('chmod 777 {0}/faToVcf'.format(filePath))


"""

If the user provided unaligned sequences, perform an MSA with mafft
Otherwise, read the name of the provided msa file

"""
if args['unaligned'] is not None:

    temp = tempfile.NamedTemporaryFile(mode='wt', delete=False)
    for name in glob.glob('{0}/*'.format(args['inpath'])):
        fasta_sequences = SeqIO.parse(open(name), 'fasta')
        for fasta in fasta_sequences:
            header, sequence = fasta.id, str(fasta.seq)
            if header not in headDic:
                headDic[header] = ''
                temp.write('>{0}\n{1}\n'.format(header, sequence))
    os.system(
        'mafft --thread {1} --auto --keeplength --addfragments {3}  {0} > {2}/inputMsa.fa'.format(
            args['reference'],
            args['thread'],
            filePath,
            temp.name))
    msaName = 'inputMsa.fa'
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

Run fastaToVcf

"""

probDic = {}

# Retrieve recomended problematic sites if the user specied and run fastaToVcf

if args['auto_mask'] is not None:

    os.system('rm -r {0}/problematic_sites_sarsCov2.vcf'.format(filePath))
    os.system('wget https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf')
    if '.vcf' in args['output']:
        os.system(
            '{3}/faToVcf -maskSites={3}/problematic_sites_sarsCov2.vcf -ref=\"{0}\" {1} {2}'.format(
                head, msaName, args['output'].strip()))
    else:
        os.system(
            '{3}/faToVcf -maskSites={3}/problematic_sites_sarsCov2.vcf -ref=\"{0}\" {1} {2}.vcf'.format(
                head, msaName, args['output'].strip()))


# Retrieve user provided problematic sites if the user specied and run
# fastaToVcf

elif args['user_specified_mask'] is not None:
    if '.vcf' in args['output']:
        os.system('{3}/faToVcf -maskSites={3}/{0} -ref=\"{1}\" {2} {3}'.format(
            args['user_specified_mask'], head, msaName, args['output'], filePath))

    else:

        os.system('{3}/faToVcf -maskSites={3}/{0} -ref=\"{1}\" {2} {3}.vcf'.format(
            args['user_specified_mask'], head, msaName, args['output'], filePath))


# Else run fastaToVcf without masking

else:

    if '.vcf' in args['output']:
        os.system('{3}/faToVcf  -ref=\"{0}\" {1} {2}'.format(head,
                                                             msaName, args['output'], filePath))

    else:

        os.system('{3}/faToVcf  -ref=\"{0}\" {1} {2}.vcf'.format(head,
                                                                 msaName, args['output'], filePath))


if os.path.isfile('{0}/problematic_sites_sarsCov2.vcf'.format(filePath)):

    os.system('rm -r {0}/problematic_sites_sarsCov2.vcf'.format(filePath))
