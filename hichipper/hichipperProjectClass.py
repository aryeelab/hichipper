import click
import os
import os.path
import sys
import random
import string
import itertools
import time
import platform
import yaml

from .hichipperHelp import *

class hichipperProject():
	def __init__(self, script_dir, mode, out, peaks, restriction_frags,
		skip_resfrag_pad, skip_background_correction):
		
		#-----------------------------------------
		# Potentially parse .yaml file if supplied
		#-----------------------------------------
		if mode.endswith(('.yaml', '.yml')):
			m = parse_manifest(mode)
			self.peaks = m['peaks'][0]
			self.resfrags = m['resfrags'][0]
			self.hicprooutput = m['hicpro_output'][0]
			self.go = "yaml"
		else:
			self.go = "call"
			self.peaks = peaks
			self.resfrags = restriction_frags
		
		#------------------------------
		# Determine if restriction fragment calling is happening or not
		#------------------------------
		if not (skip_resfrag_pad or skip_background_correction):
			if not os.path.isfile(self.resfrags):
				sys.exit('ERROR: Could not find the restriction fragment file ' + self.resfrags + '; either correctly specify file in .yaml; use the --skip-resfrag-pad and --skip-background-correction flags; or supply them directly with the `--restriction-frags` flag')
	
		#----------------------------------
		# Assign straightforward attributes
		#----------------------------------
		self.script_dir = script_dir
		self.mode = mode
		self.out = out
		
		
	#--------------------------------------------------------------------------------
	# Define a method to dump the object as a .yaml/dictionary for use in other files
	#--------------------------------------------------------------------------------
	def __iter__(self):
		
		yield 'script_dir', self.script_dir
		yield 'mode', self.mode
		yield 'output', self.output
		yield 'bamfile', self.bamfile

		yield 'cluster', self.cluster
		yield 'jobs', self.jobs
		yield 'minimum_barcode_fragments', self.minimum_barcode_fragments
		yield 'minimum_cell_fragments', self.minimum_cell_fragments
		
		yield 'extract_mito', self.extract_mito
		yield 'tssFile', self.tssFile
		yield 'blacklistFile', self.blacklistFile
		yield 'bedtoolsGenomeFile', self.bedtoolsGenomeFile
		yield 'R', self.R
		
		yield 'barcode_tag', self.barcode_tag
		yield 'bam_name', self.bam_name
		
		yield 'bowtie2', self.bowtie2
		yield 'bowtie2_index', self.bowtie2_index
		
	