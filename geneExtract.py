#!/usr/bin/env python

import os
import sys
from pybedtools import BedTool
import pysam
import subprocess
import argparse


#pass argument from command
parser = argparse.ArgumentParser(description='Extract target genes and generate bloom filter')
parser.add_argument('-p', '--prefix', dest='prefix', type=str, 
                    help='Provide output file path. '
                    'Default path: %(default)s', default=os.getcwd())
parser.add_argument('-g', '--gene', dest='gene', metavar='GeneNames_file', type=file, required=True,
                    help='Provide absolute path of text file containing target names')
parser.add_argument('-gtf', dest='gtf', metavar='GTF_file', type=file, required=True,
                    help='Provide absolute path of gtf_file')
parser.add_argument('-ref', dest='ref', metavar='Ref Genome', type=file, required=True,
                    help='Provide absolute path of the reference file')
args = parser.parse_args()

print 'Prefix', args.prefix
#generate set contains gene name without duplicate
def getGeneName():
    input = open(args.gene.name, 'r')
    geneSet = set()
    for line in input:
        line = line.replace('\n', '')
        geneSet.add(line)
    print geneSet
    return geneSet

#define dictionary containg the information
#go through gtf file extract the longest
#sequence for gene in geneSet
def geneCoordinate(geneSet):
    seqdic = {}
    log = []
    for gene in geneSet:
        for feature in BedTool(args.gtf.name):
            if feature[2] == 'exon':
                if feature.attrs.has_key('gene_id'):
                    geneName = feature.attrs['gene_id']
                    if geneName == gene:
                        print 'Extracted sequence: ', geneName
                        sys.stdout.flush()
                        startP = int(feature.start) + 1
                        endP = int(feature.stop)
                        exon = [feature[0], startP, endP]
                        print 'Coordinate', exon
                        sys.stdout.flush()
                        if gene in seqdic:
                            if startP < seqdic[gene][1]:
                                seqdic[gene][1] = startP
                            if endP > seqdic[gene][2]:
                                seqdic[gene][2] = endP
                                #update the starting and end point
                        else:
                            #add new key
                            seqdic[gene] = exon
                        #end if
                    #end if
                #end if
            #end if
        #end for
    #end for

    #print error message, genes failed to find
    for gene in geneSet:
        if gene not in seqdic:    
            log.append(gene)
    if not log:                
        print 'Successfully found all genes'
    else:
        print '>>> Fialed to find following genes {} <<<'.format(log)
        print 'Please double check your GTF file or input gene name'

    #print dictionary
    for key, list in seqdic.iteritems():
        print 'gene coordinate: ', key, list

    #return dictionary
    return seqdic
#end of function

#extract genes from reference genome and write it to output file
def extractGene(seqdic):
    target_fa = os.path.join(args.prefix, 'targetGene.fasta')
    with open(target_fa, 'w') as fi:    
        ref_fasta = pysam.Fastafile(args.ref.name)
        for key, list in seqdic.iteritems():
            print list[0], list[1], list[2]
            fi.write('{}\n'.format('>'+key))
            fi.write('{}\n'.format(ref_fasta.fetch(list[0], list[1], list[2])))
            print 'Write {} to {} file'.format(key, fi.name)
        return fi.name

#make bloom filter out of targeted genes
def makeBloomFilter(input_file):
    target_bf = os.path.join(args.prefix, 'targetGene')
    subprocess.call(['samtools', 'faidx', input_file])
    print 'Run samtools'
    cmd_pass = ['biobloommaker', '-p', 'targetGene', '-o', args.prefix, input_file]
    subprocess.call(cmd_pass)
    print 'Run biobloommaker'
    return target_bf

#main function
if __name__ == "__main__":
    target_fa = os.path.join(args.prefix, 'targetGene.fasta')
    if not os.path.isfile(target_fa):
        print 'Step one: Extract Gene Name from input file {}'.format(args.gene.name)
        geneSet = getGeneName()
        print 'Step two: Get gene coordinate from GTF file {}'.format(args.gtf.name)
        seqdic = geneCoordinate(geneSet)
        print ('Step three: Get gene sequnce from '
               'provided reference genome {}'.format(args.ref.name))
        target_fa = extractGene(seqdic)
    else:
        print ('>>>File {} already exists. To overwrite '
               'please modify output name and prefix<<<'.format(target_fa))
    print 'Step four: Generate bloom filter from {} file.'.format(target_fa)
    target_bf = makeBloomFilter(target_fa)
    print 'Successfully generate bloom filter {}'.format(target_bf)
