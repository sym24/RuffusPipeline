#A Targeted Clinical Pipeline Built in Ruffus

##Required installation:
	Python 2.7.8	
	Python Ruffus Library
	TransAbyss
	BioBloomTool
	Gmap
##Input files:
Full path of sequencing libraries contained in a text file. If the reads are in paired-end mode, they should be 		ulternating.
```
	For example:
	/abs/path/to/file/library1_1.fq
	/abs/path/library1_2.fq
	/abs/path/library2_1.fq
	/abs/path/library2_2.fq
```	
###Usage for fastq_split.py
Added this file to split fastq output from BioBloomTool 2.0.7, which will output matched read pairs into stdout with command '-d [FIlTERNAME]'\n 
Help display:
```
usage: fastq_split.py [-h] [-p PREFIX] [-g] [-f]

Write Interlaced Fastq Standard Output Into Splited Files

optional arguments:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        Provide output file path and prefix Format:
                        /path/name_split_1, /path/name_split_2 Default path:
                        /home/ymingsun/work/project_pipeline/ruffus
  -g, --gzip            Specify gzip output file. If not specified, files will
                        not be gzip.
  -f, --force           Overwriter existing file
```
