import glob
import os 
import argparse 
import re
from sequenceAnalyzer import FastAreader


parser = argparse.ArgumentParser(description='Generates merged VCF that ignores indels and problematic sites and recognizes missing data in addition to genotypes.')

print('\n\nFor more information on problematic sites see:\n\nNicola De Maio, Landen Gozashti, Yatish Turakhia, Conor Walker, Robert Lanfear, Russell Corbett-Detig, and Nick Goldman, Issues with SARS-Cov-2 sequencing data: Updated analysis with data from 12th June 2020, Virological post 2020. https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/12 \n    and \nYatish Turakhia, Bryan Thornlow, Landen Gozashti, Angie S. Hinrichs, Jason D. Fernandes, David Haussler, and Russell Corbett-Detig, "Stability of SARS-CoV-2 Phylogenies", bioRxiv pre-print 2020.\n\n')

parser.add_argument('-inpath', nargs='?', required=True,
                     help='Path to a directory containing all input multiple sequence alignments you wish to consider')
parser.add_argument('-unaligned', nargs='?',const='T', default=None,
                    help='Signifies that user provided fasta files have not been aligned to the reference. If specified, fasta2vcf uses mafft to perform multiple sequence alignments.')
parser.add_argument('-mask_problematic_sites', nargs='?',const='T', default=None,
                    help='Ignore problematic sites per masking recomendations')
parser.add_argument('-outfile', nargs='?', required=True,
                     help='Name of output file')

varDic = {}
gDic = {}
msaList = []

args = vars(parser.parse_args())
myReaderRef= FastAreader('MN908947.3.fa')
for header, sequence in myReaderRef.readFasta():
    reference = sequence.upper()


if args['unaligned'] != None:
    
    headCount = 0
    for name in glob.glob('{0}/*'.format(args['inpath'])):
        myReaderMSA= FastAreader(name)
        for header, sequence in myReaderMSA.readFasta():
            if header != 'MN908947.3':
                headCount += 1
                with open('{0}/temp.fa'.format(args['inpath']),'w') as f:
                    f.write('>{0}\n{1}'.format(header,sequence))
                os.system('cat MN908947.3.fa {0}/temp.fa > {0}/temp_premsa.fa'.format(args['inpath']))
                os.system('mafft {0}/temp_premsa.fa  > {0}/{1}_{2}.fa'.format(args['inpath'],header.replace('/','_').replace(' ','_').replace('\t','_'),str(headCount)))
                os.system('rm -r {0}/temp.fa'.format(args['inpath']))
                os.system('rm -r {0}/temp_premsa.fa'.format(args['inpath']))

                msaList.append('{0}/{1}_{2}.fa'.format(args['inpath'],header.replace('/','_').replace(' ','_').replace('\t','_'),str(headCount)))
    

else:
    for name in glob.glob('{0}/*'.format(args['inpath'])):
        msaList.append('{0}'.format(name))

seqDic = {}
altCount = {}
refCount = {}


"""

Retreave problematic sites

"""
probDic = {}
if args['mask_problematic_sites'] != None:

    os.system('rm -r problematic_sites_sarsCov2.vcf')
    os.system('wget https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf')



    with open('problematic_sites_sarsCov2.vcf', 'r') as f:
        for line in f:
            if '#' in line:
                pass
            else:
                cols = line.split('\t')

                if cols[6].strip() == 'mask':
                    probDic[str(cols[1].strip())] = ''


"""

Identify variants and missing data in input alignments
Ignore problematic sites

"""

for aln in msaList:
    seqDic = {}
    myReaderMSA= FastAreader(aln)


    for header, sequence in myReaderMSA.readFasta():
        
        if header == 'MN908947.3' or sequence.replace('-','').upper() == reference:
            refHead = header
            refSeq = sequence.upper()
        seqDic[header] = sequence.upper()
        myTestSeq = sequence.upper()

    for i in range(len(refSeq)):
        ref = refSeq[i]
        if i not in refCount:
            refCount[i] = 0
        if i not in altCount:
            altCount[i] = 0
        if i not in gDic:
            gDic[i] = {}
        if i not in varDic:
            varDic[i] = []
     


        for header in seqDic:
            if seqDic[header][i] == ref or str(i+1) in probDic: #or ref == '-':
                refCount[i] += 1
                gDic[i][header] = 0
            elif seqDic[header][i].upper() == 'N' or seqDic[header][i] == '-' or ref == '-':
                #varDic[i][header] = '.'
                gDic[i][header] = '.'
               # print(i,header)
            else:
                if seqDic[header][i] not in varDic[i]:
                    varDic[i].append(seqDic[header][i])
                    
                altCount[i] += 1

                gDic[i][header] = varDic[i].index(seqDic[header][i])+1

"""

Write merged VCF

"""

with open('{0}'.format(args['outfile']),'w') as f:
    f.write('##fileformat=VCFv4.2\n##source_20200913.1=vcf-merge(r953) {0}\n##sourceFiles_20200913.1={0}\n##nanopolish_window=MN908947.3:1-29902\n##INFO=<ID=TotalAltCalls,Number=1,Type=Integer,Description="The number of alternate allele calls for the variant">\n'.format(','.join(msaList)))
    f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{0}\n'.format('\t'.join(gDic[0].keys())))
    for site in gDic:
        if varDic[site]:
            stringList = [str(x) for x in gDic[site].values()]
            f.write('{0}\t{1}\t.\t{2}\t{3}\t.\t.\tAC={4}\tGT\t{5}\n'.format(refHead,str(site+1),refSeq[site].upper(),','.join(varDic[site]).upper(),altCount[site],'\t'.join(stringList)))
    

            



