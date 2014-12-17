#!/usr/bin/env python

import sys
import argparse
import os
import gzip

#pass argument from command
parser = argparse.ArgumentParser(description='Write Interlaced Fastq '
                                 'Standard Output Into Splited Files')
parser.add_argument('-p', '--prefix', dest='prefix', type=str, 
                    help='Provide output file path and prefix  '
                    'Format: /path/prefix_1, /path/prefix_2 '
                    'Default path: %(default)s',
                    default=os.path.join(os.getcwd(), 'matchedReads'))
parser.add_argument('-g', '--gzip', dest='gzip',action='store_true',
                    help='Specify gzip output file. '
                    'If not specified, files will not be gzip.')
parser.add_argument('-f', '--force', dest='overwrite', action='store_true',
                    help='Overwriter existing file')
args = parser.parse_args()

outpath = os.path.dirname(args.prefix)
outname = os.path.basename(args.prefix)

#check if prefix exists
if not os.path.exists(outpath):
    os.makedirs(outpath)
    print 'Make output directory ', outpath

print 'File path: ', outpath                    

def checkExist(file1, file2):
    if(os.path.isfile(file1) or os.path.isfile(file2)):
        if(args.overwrite):
            print 'File {}, {} exist, will overwrite'.format(file1, file2)
        else:
            sys.exit('File {}, {} exist.\nPlease use "-f" to overwrite or '
                     'change output file name or path'.format(file1, file2))
        #end if
    #end if
   
if(args.gzip):
    #define output files  
    name1 = os.path.join(outpath, outname + '_1.fq.gz')
    name2 = os.path.join(outpath, outname + '_2.fq.gz')
    checkExist(name1, name2)
    f1 = gzip.open(name1, 'w')
    f2 = gzip.open(name2, 'w')
else:
    #define output files  
    name1 = os.path.join(outpath, outname + '_1.fq')
    name2 = os.path.join(outpath, outname + '_2.fq')
    checkExist(name1, name2)
    f1 = open(name1, 'w')
    f2 = open(name2, 'w')

#define flag variables
buff = ""
counter = 0
flipflot = 1

#define debug flags
recs1 = 0
recs2 = 0

try:
    for line in sys.stdin:
        buff += line
        counter += 1
        if(counter==4):
            if(flipflot==1):
                f1.write(buff)
                recs1 += 1
            else:
                f2.write(buff)
                recs2 += 1
            counter = 0
            buff = ""
            flipflot *= -1
    f1.close()
    f2.close()
except KeyboardInterrupt:
    raise

if(recs1 != recs2):
    print 'The number of processed records does not match - check input data!'
print ("Processed {} records for pair file one "
       "and {} records for pair file two.".format(recs1, recs2))
#end of file
