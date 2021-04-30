#!/usr/bin/env python3

import sys, os, re, argparse, textwrap

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\

Compares two VCF files and outputs info if alleles don't match.
IMPORTANT: The first VCF is expected to use only TCGAN, the second
can have IUPAC codes. The second file is also allowed to contain
positions that do not occur in the first (but not vice versa).

Samples MUST be in the same order in both files, use vcfSort first.

To vastly reduce run time, please split files so that they 
contain about 100k samples each, for instance like so:

for i in $(seq 100001 100000 1300000); do
    # 1 through 9 are VCF info fields
    let k=$i+99999
    cut -f1-9,$i-$k file1.vcf  > splitVcf/vcf.tcga.$i 
    cut -f1-9,$i-$k file2.vcf  > splitVcf/vcf.iupac.$i 
done
The code will run in about an hour on each file pair.

        '''))

group = parser.add_argument_group('required arguments')
group.add_argument('tcgafile', type=str, help='VCF with TCGA')
group.add_argument('iupacfile', type=str, help='VCF with IUPAC')

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

def getToStart(fh):
    '''Loop through file with readline until first non comment line'''
    for line in fh:
        if line.startswith('#CHROM'):
            # samples start in field 10
            snames = line.strip().split('\t')[9:]
        if line.startswith('#'):
            continue
        else:
            return snames,line, fh

def alleleComp(sampleIds, apos, aref, alts, bref, blts, iupac, asamples, bsamples):
    '''Compares samples between two sets a and b, where a may contain iupac alleles'''
    problemCases = 0
    for i in range(len(asamples)):
        # if one of the positions has N, skip (N matches everything by definition)
        if asamples[i] == '.' or bsamples[i] == '.':
            #print("position {} sample number {} has tcga allele {} and ambig allele {}".format(apos, i, asamples[i],  
            #  blts[int(bsamples[i])]), file=sys.stderr)
            continue
        allele = alts[int(asamples[i])]
        bllele = blts[int(bsamples[i])]
        if allele != bllele:
            if not bllele in iupac[allele]:
                problemCases +=1
                print("Disagreement at position {} in sample {}: {} vs {} (iupac: {})".format(apos, 
                  sampleIds[i], bllele, allele, iupac[allele]))
    return problemCases

iupac = {
'A': ['A'],
'T': ['T'],
'C': ['C'],
'G': ['G'],
'R': ['A', 'G'],
'Y': ['C', 'T'],
'S': ['G', 'C'],
'W': ['A', 'T'],
'K': ['G', 'T'],
'M': ['A', 'C'],
'B': ['C', 'G', 'T'],
'D': ['A', 'G', 'T'],
'H': ['A', 'C', 'T'],
'V': ['A', 'C', 'G'],
}

# Main

problems = 0
with open(args.tcgafile, 'r') as before:
  with open(args.iupacfile, 'r') as after:
     # files may not have same number of header lines
     sampleIds, bline, before = getToStart(before)
     for aline in after:
        if aline.startswith('#'):
            continue
        afields = aline.strip().split("\t")
        bfields = bline.strip().split("\t")
        # first nine fields are info
        apos, aref, altstr = afields[1], afields[3], afields[4]
        # prepend ref allele to list
        alts = [aref] + altstr.split(',')
        # same for the tcga file
        bpos, bref, bltstr = bfields[1], bfields[3], bfields[4]
        blts = [bref] + bltstr.split(',')
        # we're out of sync if one of the files has a position the other does not
        # we expect this only in the iupac file (the after set)
        # so when that happens we need to skip the line and NOT MOVE ON in the before set
        if bpos != apos:
#            print("Files out of sync at pos {} {}".format(apos, bpos), file=sys.stderr)
            continue
        # now loop through samples
        problems += alleleComp(sampleIds, apos, aref, alts, bref, blts, iupac, afields[9:], bfields[9:])

        # read another line
        bline = before.readline()
        
before.close
after.close
print("Number of problems found: {}".format(problems))

