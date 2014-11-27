#!/usr/bin/env python

########import libraries########
import os
import sys
from ruffus import *
from subprocess import Popen, PIPE
import argparse
import subprocess
import time

########start timing########
startTime = time.time()

########get software path########
try:
    GMAP_PATH=os.environ['GMAP_PATH']
except KeyError as detail:
    print "Software installation path {} does not exist!".format(detail)
    raise KeyError(key)

#set up default library
def defaultwd():
    """Return the default output directory."""
    return os.path.join(os.getcwd(), "ruffus_pipeline")

parser = argparse.ArgumentParser(description='A clinical pipeline '
                         'written in Ruffus.', prog='AML Pipeline')
#specify global variables used in all commands
global_group = parser.add_argument_group("Global Variable")
global_group.add_argument('--threads', dest='threads', 
                          help='Number of threads to use. Default [%(default)s]', 
                          metavar='INT', type=int, default=1)
global_group.add_argument('--workdir', dest='workdir', 
                          help='Work directory. Program will store running results ' 
                          'under specified work directory. Default [%(default)s]', 
                          metavar='PATH', type=str, default=defaultwd())
#variables specify input options    
input_group = parser.add_argument_group("Input Options")
input_group.add_argument('--intxt', dest='intxt', metavar='.txt FILE', 
                         help='Required. Specify input reads list using '
                         'a txt file with specific name "<library_name>.in".', 
                         type=file, required=True)
input_group.add_argument('--pfilter', dest='pfilter', metavar='STR', 
                         help='The name of the filter created from '
                         'targeted gene sequence. Default: %(default)s', 
                         type=str, default='genetargets')
input_group.add_argument('--target', dest='target', metavar='.fa FILE', 
                         help='Required. BioBloommaker creates a bf and '
                         'txt file from a list of fasta files', 
                         type=file, nargs='+', required=True)
input_group.add_argument('--klist', dest='klist', metavar='INT', 
                         help='Input a list of key values. Default %(default)s', 
                         type=int, nargs='+', default=[35, 55, 75])
input_group.add_argument('--SS', dest='ss', 
                         help='Input reads are strand-specific', 
                         action='store_true', default=False)
input_group.add_argument('--se', dest='single', 
                         help='Single-end read files. '
                         'Default: Take paired-end read files', 
                         action='store_true', default=False)
#varibales specify output options
output_group = parser.add_argument_group("Output Options")
output_group.add_argument('--out', dest='outfile', 
                           help='Output file name. Defaut [%(default)s]', 
                           metavar='PATH', type=str, default='output')
output_group.add_argument('--cleanup', dest='cleanup', choices=[0,1,2,3], 
                          help='Level of clean up intermediate files. '
                          'Default [%(default)s]', type=int, default=3)

args = parser.parse_args()

###let us try to keep variables private to avoid mess up names

output_name = os.path.basename(args.intxt.name).split('.')[0]
print 'Output name configured: ', output_name
print 'Working directory: ', args.workdir

#make a list containing target gene files
targets_list = [args.target[i].name for i in range(len(args.target))]
print targets_list

def defworkdir():
    """Check if work directory exists. 
    If not, creat work directory
    """
    if not os.path.exists(args.workdir):
        os.makedirs(args.workdir)
#end function

def formworkdir():
    """Format work directory input by user."""
    if args.workdir.endswith("/"):
        return args.workdir
    else:
        return args.workdir + '/'
  
def filter_maker_outlist():
    """Generate biobloommaker output list."""
    outlist=[]
    filepath = formworkdir() + args.pfilter
    outlist.append(filepath+'.bf')
    outlist.append(filepath+'.txt')
    print 'bloom filter maker output list: ', outlist
    return outlist
#end function

def filter_cater_outlist():
    """Generate biobloomcategorizer output list."""
    cater_1 = formworkdir() + output_name + '_' + args.pfilter + '_1.fq' 
    cater_2 = formworkdir() + output_name + '_' + args.pfilter + '_2.fq'
    print 'bloomfilter categorizer output list: ', [cater_1, cater_2]
    return [cater_1, cater_2]
#end funtion

def assemble_outlist():
    """Generate transabyss assembly output list."""
    outlist = []
    for k in args.klist:
        direct = formworkdir() + 'k' + str(k)
        if not os.path.exists(direct):
            os.makedirs(direct)
        name = direct + '/' + output_name + '-final.fa'
        outlist.append(name)
    #end for
    print 'transabyss assembly output list: ', outlist
    return outlist
#end function

def sort_read_list(input_name):
    """Sort input files into two lists 
    containing corresponding paired reads files
    """
    print 'input file name is: ', input_name
    
    files_list = args.intxt.read().splitlines()
    files_list = filter(None, files_list)
    print 'files from txt: ', files_list
    
    file_1 = []
    file_2 = []
    file_1 = [files_list[i] for i in range(len(files_list)) if i % 2 == 0]
    file_2 = [files_list[i] for i in range(len(files_list)) if i % 2 == 1]
    file_1.insert(0, 'zcat')
    file_2.insert(0, 'zcat')
    return file_1, file_2
#end function

def filter_cater_run(input_name, bf_filter):
    """Functions to run biobloomcategorizer"""
    #get concantonate files 
    print 'parameter passed: ', input_name
    file_1, file_2 = sort_read_list(input_name)
    
    #specify output directory
    input_name = os.path.basename(args.intxt.name)
    outfile_dir = args.workdir + '/' + input_name.split('.')[0] 
    print 'categorizer output directory: ', outfile_dir

    input_fd1, output_fd1 = os.pipe()
    Popen(file_1, shell=False, stdout=output_fd1)
    os.close(output_fd1)
    print 'input_fd1:', str(input_fd1)

    input_fd2, output_fd2 = os.pipe()
    Popen(file_2, shell=False, stdout=output_fd2)
    os.close(output_fd2)
    print 'input_fd2:', str(input_fd2)

    if(args.single):
    #need to modify the following command follow single end reads. Should 
    #replace the last part with a file_list
        p = Popen(['biobloomcategorizer', '--fq', '-p', str(outfile_dir), 
                   '-t', str(args.threads), '-f', str(bf_filter), 
                   '/dev/fd/'+str(input_fd1), '/dev/fd/'+str(input_fd2)], 
                  shell=False)
        p.wait()
    else:
	p = Popen(['biobloomcategorizer', '--fq', '-p', str(outfile_dir), 
		   '-t', str(args.threads), '-e', '-f', str(bf_filter), 
		   '/dev/fd/'+str(input_fd1), '/dev/fd/'+str(input_fd2)], 
		  shell=False) 
	p.wait()
	#subprocess.popen must take string as argument. That is why str(prefix)
    print 'finish running biobloomcategrizer'
#end function

def abyss_assemble(input_file, output_file):
    """Functions to run transabyss assembly"""
    print input_file
    outlist = 'flag'
    cmd_pass = ['transabyss', input_file[0], input_file[1], 
		'--name', str(output_name), 
		'--threads', str(args.threads), 
		'--cleanup', str(args.cleanup)
		]
    print cmd_pass
    if args.ss:
	cmd_pass.append('--SS')
	print '--SS'

    if args.single:
	cmd_pass.insert(1, '--se')
	print '--se'
    else:
	cmd_pass.insert(1, '--pe')
	print '--pe'
    
    for outfile, k in zip(output_file, args.klist):
	outdir = os.path.dirname(outfile)
	print outdir	    
	cmd_copy = cmd_pass[:]
	cmd_copy.append('--outdir')
	cmd_copy.append(outdir)
	cmd_copy.append('-k')
	cmd_copy.append(str(k))
	print cmd_copy
	subprocess.call(cmd_copy)      ###subprocss  
    #endfor
#end funtion

@transform(targets_list, formatter(".fa$"), 
	   "{path[0]}/{basename[0]}.fa.fai")
def run_samtools(input_file, output_file):
    """Use Samtools to index target genes. """
    print 'input_files = ', input_file
    print 'output_files = ', output_file 
    subprocess.call(['/projects/btl/arch/samtools', 'faidx', input_file])
    print '--successfully genereate indexed file'
#end of samtools

@follows(run_samtools)
@files(targets_list, 
       filter_maker_outlist(), args.pfilter
       )
def biobloom_filter(input_file, output_file, prefix):
    """Run biobloommaker to make bloomfilter out of target genes"""
    print 'input_files = ', input_file
    print 'output_files = ', output_file
    cmd_pass = ['biobloommaker', '-p', prefix, '-o', args.workdir]
    for filters in input_file:
         cmd_pass.append(filters)
    print 'command passed: ', cmd_pass
    subprocess.call(cmd_pass)
#end of biobloommaker

###make directory for each individual library to store out put value 
#before running biobloomcategorizer
@transform(biobloom_filter, 
           formatter(),
           filter_cater_outlist()
           )
def biobloom_cater(input_file, output_file):
    """Run biobloom catergrizer to fish out 
    matching reads from libraries
    """
    #remove empty strings in the list
    print 'input_files = ', input_file
    print 'output_parameters = ', output_file
    #call function to run biobloom categorizer
    filter_cater_run(args.intxt.name, input_file[0])
#end of biobloomcategorizer

@transform(biobloom_cater,
           formatter(".+/(?P<LIBN>\w+)_(?P<STR>\w+)_1.fq"),
           assemble_outlist(),
           )      
def assemble_transAbyss(input_file, output_file):
    """Run transabyss assembly multiple times based on k-mer sizes"""
    print 'input_files = ', input_file
    print 'output_parameters = ', output_file 
    #call function to run transabyss assembly
    abyss_assemble(input_file, output_file)
#end of transabyss assembly

@mkdir(assemble_transAbyss, 
       formatter(),
       "{subpath[0][1]}/merged"
       )
@transform(assemble_transAbyss, 
           formatter(".+/(?P<LIBN>\w+)-final.fa"), 
           "{subpath[0][1]}/merged/{LIBN[0]}_"
           "{subdir[0][0]}_{subdir[1][0]}_{subdir[2][0]}-merged.fa",
           args.klist
           )  
def merge_transAbyss(input_file, output_file, k_values):   
    """Run transabyss merge to merge multiple k-mer assemblies"""
    cmd_pass = ['transabyss-merge', '--mink', str(min(k_values)), 
                '--maxk', str(max(k_values)), '--out', output_file, 
                '--threads', str(args.threads)
                ]
    for k in k_values:
        if('--prefixes' not in cmd_pass):
            cmd_pass.append('--prefixes')
        cmd_pass.append('.' + str(k))
        print cmd_pass
    #end for
    if(args.ss):
        cmd_pass.append('--SS')
        print 'Strand specific'
    for files in input_file:
        if(files not in cmd_pass):
            cmd_pass.append(files)
    #end for
    print cmd_pass
    subprocess.call(cmd_pass)
    print 'input_files = ', input_file
    print 'output_parameters = ', output_file 
    print 'k_values = ', k_values
#end of function

@mkdir(merge_transAbyss, formatter(), "{subpath[0][1]}/alignment")
@transform(merge_transAbyss, formatter(), 
           "{subpath[0][1]}/alignment/c2g.bam" 
           )
def gmap_alignment(input_contigs, output_file):
    """Align contiges to genome using gmap"""
    cmd = ("{} -d hg19 "
    "-D /projects/btl/arch/gmapdb_sarray/hg19/ {} "
    "-t {} -f samse -O -x 10 -n 0 | "
    "grep -v chrM | /projects/btl/arch/samtools "
    "view -bhS - -o {}").format(GMAP_PATH, input_contigs, str(args.threads), output_file)
    print 'command pass to gamp: ', cmd
    subprocess.call(['/bin/sh', '-c', cmd])
    
    print 'input contigs: ', input_contigs
    print 'output file directory: ', output_file

#startTime = time.time()
totalTime = time.time() - startTime
print "Programming Running time: ", totalTime
#pipeline_run(multiprocess=2, touch_files_only=False,
#             history_file=".history_file=.ruffus.command2", checksum_level=1)

pipeline_run(multiprocess=2, touch_files_only=False,
             history_file=".history_file=.ruffus.command2", checksum_level=1)

#pipeline_printout(sys.stdout, gmap_alignment, checksum_level=1, history_file=".history_file=.ruffus.command2", verbose=6, verbose_abbreviated_path=0)
#CHECKSUM_level = 0     # only rerun when the file timestamps are out of date (classic mode)
#CHECKSUM_level = 1     # Default: also rerun when the history shows a job as being out of date


