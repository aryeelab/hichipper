import click
import os
import os.path
import sys
import shutil
import yaml
import shutil
import random
import string
import time
import re
from pkg_resources import get_distribution
from subprocess import call, check_call
from .hichipperHelp import *
from .hichipperProjectClass import *
from .hicproHelper import *

@click.command()
@click.version_option()

# User output / mode
@click.argument('mode')
@click.option('--out', "-o", default="hichipper_out", required=True, help='Output directory name; must not be already existing [Required]')
@click.option('--keep-temp-files', "-z", is_flag=True, help='Keep temporary files?')

# User input
@click.option('--input-vi', '-ii', default = "", help='Comma-separted list of interactions files for loop calling; option valid only in `call` mode')
@click.option('--restriction-frags', '-rf', default = "", help='Filepath to restriction fragment files; will overwrite specification of this file when a .yaml is supplied for mode')
@click.option('--peaks', '-p', default = "", help='Either 1 of 4 peak logic strings or a valid filepath to a .bed (or similary formatted) file; defers to what is in the .yaml')

# Essential options
@click.option('--keep-samples', "-k", default="ALL", help='Comma separated list of sample names to keep; ALL (special string) by default')
@click.option('--ignore-samples', "-x", default="NONE", help='Comma separated list of sample names to ignore; NONE (special string) by default')
@click.option('--read-length', "-l", default="75", help='Length of reads from sequencing runs; default = 75')

# Loop Distance options
@click.option('--min-dist', "-mi", default="5000", help='Minimum distance for loop calls; default = 5000')
@click.option('--max-dist', "-ma", default="2000000", help='Maximum distance for loop calls; default = 2000000')

# MACS2 Configurations
@click.option('--macs2-string', default="-q 0.01 --extsize 147 --nomodel", help='String of arguments to pass to MACS2; only is called when peaks are set to be called; default = "-q 0.01 --extsize 147 --nomodel"')
@click.option('--macs2-genome', default="hs", help='Argument to pass to the -g variable in MACS2 (mm for mouse genome; hs for human genome); default = "hs"')

# Loop anchor options
@click.option('--peak-pad', "-pp", default="500", help='Peak padding width (applied on both left and right); default = 500')
@click.option('--merge-gap', "-mg", default="500", help='Merge nearby peaks (after all padding is complete; default = 500')
@click.option('--no-merge', "-nm", is_flag=True, help='Completely skip anchor merging; will affect summary statistics. Not recommended unless understood what is happening.')

@click.option('--skip-resfrag-pad', is_flag=True, help='Skip restriction fragment aware padding')
@click.option('--skip-background-correction', is_flag=True, help='Skip restriction fragment aware background correction?')

# External Dependencies
@click.option('--make-ucsc', "-mu", is_flag=True, help='Make additional output files that can support viewing in UCSC genome browser; requires tabix and bgzip; does the same thing as --make-washu.')
@click.option('--make-washu', "-mw", is_flag=True, help='Make additional output files that can support viewing in WashU genome browser; requires tabix and bgzip; does the same thing as --make-ucsc.')

@click.option('--basic-qc', is_flag=True, help='Create a simple QC report without Pandoc')
@click.option('--skip-diffloop', is_flag=True, help='Skip analyses in diffloop (e.g. Mango loop calling; .rds generation)')

# Software options
@click.option('--bedtools-path', default = "", help='Path to bedtools; by default, assumes that bedtools is in PATH')
@click.option('--macs2-path',    default = "", help='Path to macs2; by default, assumes that macs2 is in PATH')
@click.option('--tabix-path',    default = "", help='Path to samtools; by default, assumes that tabix is in PATH')
@click.option('--bgzip-path',    default = "", help='Path to macs2; by default, assumes that bgzip is in PATH')
@click.option('--r-path',        default = "", help='Path to R; by default, assumes that R is in PATH')


def main(mode, out, keep_temp_files, 
	input_vi, restriction_frags, peaks,
	keep_samples, ignore_samples, read_length,
	min_dist, max_dist,
	macs2_string, macs2_genome,
	peak_pad, merge_gap, no_merge, skip_resfrag_pad, skip_background_correction,
	basic_qc, skip_diffloop, make_ucsc, make_washu,
	bedtools_path,  macs2_path, tabix_path, bgzip_path, r_path):
	
	
	"""
	hichipper: a preprocessing and QC pipeline for HiChIP data. \n
	(c) Aryee Lab, 2019 \n
	See https://hichipper.readthedocs.io for more details.\n
	
	hichipper mode: [call, *.yaml]
	^ either specify the word `call` and feed in a valid interactions file
	OR specify the .yaml format for options to be parsed from a manifest file (see documentation)
	"""
	
	# Staples
	__version__ = get_distribution('hichipper').version
	
	#-------------------------------------------
	# Initial verification for external software
	#-------------------------------------------
	if(make_ucsc or make_washu):
		tabix = get_software_path("tabix", tabix_path)
		bgzip = get_software_path("bgzip", bgzip_path)
	else:
		tabix = ""
		bgzip = ""
	bedtools = get_software_path("bedtools", bedtools_path)
	
	# See if R is necessary
	if( not (skip_background_correction and skip_resfrag_pad and basic_qc and skip_diffloop)):
		Rscript = get_software_path("R", r_path) + "script"
	else: 
		Rscript = ""
	
	halfLength = int(float(read_length)/2)
	if not halfLength> 0:
		sys.exit('ERROR: Specify appropriate read length > 1; QUITTING')
	
	#------------------------------
	# Handle initial QC reporting
	#------------------------------
	if os.path.exists(out):
		sys.exit("ERROR: Output path (%s) already exists; remove it or specify a new location." % out)
	os.mkdir(out)
	logf = open(out + "/" + out + ".hichipper.log", 'w')
	click.echo(gettime() + "Starting hichipper pipeline v%s" % __version__, file = logf)
	click.echo(gettime() + "Starting hichipper pipeline v%s" % __version__)
	
	# Handle directories
	script_dir = os.path.dirname(os.path.realpath(__file__))
	outfolder = os.path.abspath(out)  
	cwd = os.getcwd()
	click.echo(gettime() + "Executed from: %s" % cwd, logf)
	click.echo(gettime() + "Output folder: %s" % outfolder, logf) 
	
	click.echo(gettime() + "Parsing user parameters")
	
	 # Check for UCSC/WASHU Specification
 	if make_ucsc or make_washu:
 		ucscoutput = "true"
 	else:
 		ucscoutput = "false"
	
	# Check no merge specification
	if no_merge:
 		no_merge_str = "true"
 	else:
 		no_merge_str = "false"
	
	#------------------------------
	# If it is a manifest file, handle it as such; otherwise check for the loop call mode
	#------------------------------
	if mode.endswith(('.yaml', '.yml')):
		m = parse_manifest(mode)
		click.echo(gettime() + ".yaml file detected")
		click.echo(gettime() + "Parsed manifest as follows: ", logf)
		click.echo(m, logf)
	elif(mode == "call"):
		click.echo(gettime() + "Direct loop call option detected.")
		if not os.path.exists(out):
			sys.exit("ERROR: Peaks file (%s) cannot be found; with the `call` option, one must supply peaks" % out)
	else:
		sys.exit("ERROR: Mode option (%s) is invalid. Choose either 'call' or specify a valid path to a .yaml file." % mode)

	# Project
	p = hichipperProject(script_dir, mode, out, peaks, restriction_frags,
		skip_resfrag_pad, skip_background_correction)
	peaks = p.peaks
		
	peakopts = ["COMBINED,ALL", "EACH,ALL", "COMBINED,SELF", "EACH,SELF"]
	if(p.peaks in peakopts):
		macs2 = get_software_path("macs2", macs2_path)
	elif not os.path.isfile(peaks):
		sys.exit('ERROR: Could not identify the ' + peaks + ' file; correctly specify file location or use special variable {COMBINED,EACH},{ALL,SELF} for peaks')

	if(p.go == "yaml"):
		
		# Get samples
		samples = samplesHelper(keep_samples, ignore_samples, p.hicprooutput, logf)
		click.echo(gettime() + "Determined that the following samples are good to go: ", logf)
		click.echo(samples, logf)

		# Set up peaks files if necessary; moreover, specify vector of peaks / sample
		click.echo(gettime() + "User defined peaks specification: " + peaks, logf)
		peakfilespersample = peakHelper(p.peaks, p.hicprooutput, p.resfrags, halfLength, peak_pad, out, samples,
			Rscript, skip_resfrag_pad, skip_background_correction,
			logf, macs2_string, macs2_genome, script_dir)
		logf.close()

		# Call putative interactions
		for i in range(len(samples)):
			hichipperRun = os.path.join(script_dir, 'interactionsCall.sh')	
			cmd = ['bash', hichipperRun, cwd, out, p.hicprooutput, samples[i], peakfilespersample[i], min_dist, max_dist, merge_gap, str(halfLength), ucscoutput, no_merge_str]        
			call(cmd)
			if not os.path.isfile(out + "/" + samples[i] + ".stat"):
				sys.exit('ERROR: something failed at the individual sample level; check the .log file for more info')

	else:
		# do the new implementation for `call`
		print(input_vi)
		if(os.path.isfile(input_vi)):
			click.echo(gettime() + "Verified valid interations file: %s" % input_vi, logf) 
		else:
			sys.exit('ERROR: in `call` mode, specify `--input-vi`')
		
		samples = ["one"]
		hicprooutput = ""
		
		# Most of this isn't needed but we'll call it for simplicity with the other parameters
		peakfilespersample = peakHelper(p.peaks, hicprooutput, p.resfrags, halfLength, peak_pad, out, samples,
			Rscript, skip_resfrag_pad, skip_background_correction,
			logf, macs2_string, macs2_genome, script_dir)
		click.echo(gettime() + "Pulling interaction PETs from valid interactions file (rather than full HiC-pro output): " + peaks, logf)
		logf.close()
		hichipperRunFrags = os.path.join(script_dir, 'interactionsCall_inputFrags.sh')	
		i = 0
		cmd = ['bash', hichipperRunFrags, cwd, out, input_vi, samples[i], peakfilespersample[i], min_dist, max_dist, merge_gap, str(halfLength), ucscoutput, no_merge_str]        
		call(cmd)
		if not os.path.isfile(out + "/" + samples[i] + ".stat"):
			sys.exit('ERROR: something failed at the individual sample level; check the .log file for more info')
	

	# Back to Python
	logf = open(out + "/" + out + ".hichipper.log", 'a')
		
	# QC Report			
	if basic_qc:
		click.echo(gettime() + "Skipping QC report generation since --skip-qc was specified", logf)
	else:
		click.echo(gettime() + "Creating QC report", logf)
		cftg = ' '.join(samples)
		call(['cp', os.path.join(script_dir, 'qcReport_make.Rmd'), out + '/qcReport_make.Rmd'])
		call(['cp', os.path.join(script_dir, 'qcReport.R'), out + '/qcReport.R'])
		cmd = [Rscript, out + '/qcReport.R', script_dir, str(out), cwd, __version__ , cftg] 		
		click.echo(gettime())
		click.echo(cmd)
		call(cmd)
		
	# diffloop work
	if skip_diffloop:
		click.echo(gettime() + "Skipping diffloop analyses since --skip-diffloop was specified", logf)
	else:
		click.echo(gettime() + "Creating .rds and .mango files", logf)
		cftg = ' '.join(samples)
		cmd = [Rscript, os.path.join(script_dir, 'diffloop_work.R'), cwd, out, cftg] 
		call(cmd)	
	
	# Temporary File Management
	if keep_temp_files:
		click.echo(gettime() + "Temporary files not deleted since --keep-temp-files was specified", logf)
	else:
		click.echo(gettime() + "Deleting temporary files", logf)
		files = os.listdir(outfolder)
		for file in files:
			if file.endswith(".tmp"):
				os.remove(os.path.join(outfolder,file))
			elif "_temporary_" in file:
				os.remove(os.path.join(outfolder,file))
			elif "qcReport" in file:
				os.remove(os.path.join(outfolder,file))
	click.echo(gettime() + "Done", logf)
	
	logf.close()
