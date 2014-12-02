#!/usr/bin/env python

import sys
import argparse
import os

parser = argparse.ArgumentParser(description = 'Write Stdout into seperate outputs')
parser.add_argument('-p', '--prefix', dest='prefix', type=str, 
                    help='Provide output file path\n '
                    'Format: /prefix/split_1.fq /prefix/split_2.fq\n '
                    'Default:%(default)s', default=os.getcwd())

args = parser.parse_args()

print 'File prefix: ', args.prefix

#define output files  
name1 = os.path.join(args.prefix, 'split_1.fq')
name2 = os.path.join(args.prefix, 'split_2.fq')
f1 = open(name1, 'w')
f2 = open(name2, 'w')

o#define flag variables
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

print ("Processed {} records for pair file one "
       "and {} records for pair file two.".format(recs1, recs2))
if(recs1 != recs2):
    print 'The number of processed records does not match - check input data!'
#end of file
