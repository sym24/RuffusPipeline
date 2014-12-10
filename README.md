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
Added this file to split fastq output from BioBloomTool 2.0.7, which will output matched read pairs into stdout with command '-d [FIlTERNAME]' 

Help display:

```
usage: fastq_split.py [-h] [-p PREFIX] [-g] [-f]

Write Interlaced Fastq Standard Output Into Splited Files

optional arguments:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        Provide output file path and prefix Format:
                        /path/prefix_1, /path/prefix_2 Default path: /home/ymi
                        ngsun/work/project_pipeline/ruffus/matchedReads
  -g, --gzip            Specify gzip output file. If not specified, files will
                        not be gzip.
  -f, --force           Overwriter existing file
```
###Usage for geneExtract.py
geneExtract.py will take interested gene names and find out the longest contig coordinates using GTF file. It will then extract sequence from reference genome and build bloom filter out of targeted genes. 

Help display:

```
usage: geneExtract.py [-h] [-p PREFIX] -g GeneNames_file -gtf GTF_file -ref
                      Ref Genome

Extract target genes and generate bloom filter

optional arguments:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        Provide output file path. Default path:
                        /projects/btl/pipeline/geneExtract
  -g GeneNames_file, --gene GeneNames_file
                        Provide absolute path of text file containing target
                        names
  -gtf GTF_file         Provide absolute path of gtf_file
  -ref Ref Genome       Provide absolute path of the reference file
```

