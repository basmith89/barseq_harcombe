#!/usr/bin/python
import subprocess
import argparse,textwrap 
import re
import os
import glob

#######It might be useful to add a qsub submission portion to run on a cluster

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description = textwrap.dedent('''\
                                 Author: Brian A. Smith
                                 University of Minnesota
                                 smithba@umn.edu
		
		This is a driver for running multiple samples with 
		the Barseq script MultiCodes.pl.  It requires a file 
		with a tab separated file containing 
		[sample name]/[fastq file names]/[sequencing indexes].  

		Note: unzipped fastq files are required for Multicodes.pl
            	 '''))

parser.add_argument("-i", "--info_file", required = True,
           help = "File containing filename prefixes and indexes required")
parser.add_argument("-p", "--program_path", required = True,
           help = "Provide the path to MultiCodes.pl")
parser.add_argument("-f", "--fastq_dir", required = False,
           help = "Provide the directory where fastq files are stored if not executing within that directory")
parser.add_argument("-c", "--create_out_dir", required = False,
           help = "Creates an output directory to store data")

args = parser.parse_args()

#input processing for user friendliness
if args.fastq_dir:
	if args.fastq_dir.endswith("/"):
		fastq_dir = args.fastq_dir
	elif not args.fastq_dir.endswith("/"):
		fastq_dir = args.fastq_dir + "/"
else:
	fastq_dir = ""

if args.create_out_dir:
	if args.create_out_dir.endswith("/"):
		os.makedirs(args.create_out_dir)
		out_dir = args.create_out_dir
	elif not args.create_out_dir.endswith("/"):
		os.makedirs(args.create_out_dir)
		out_dir = args.create_out_dir + "/"
else:
	out_dir = ""

#Reading input file and running MultiCodes.pl
with open(args.info_file, 'r') as infile:
	for line in infile:
		line = line.strip("\n")
		line = line.split("\t")
		#debug print statements
		#print(line[0])
		#print(line[1])
		file_prefix = line[0]
		file_name = line[1]
		seq_index = line[2]
		#file_name = glob.glob(file_prefix + "*_R1_*.fastq")
		#print(file_name)
		multicodes_cmd = "perl " + args.program_path + " -out " + out_dir + file_prefix + " -index " + seq_index + " < " + fastq_dir + file_name
		#multicodes_cmd = "perl " + args.program_path + " -out " + out_dir + file_prefix + " -index " + seq_index + " < " + file_prefix
		print(multicodes_cmd)
		exit_status = subprocess.call(multicodes_cmd, shell=True)
		if exit_status is 1:
			print("Job " + seq_index + " failed to run".format(multicodes_cmd))
print("Driver is done.")