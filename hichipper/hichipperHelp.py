import glob
import itertools
import time
from whichcraft import which
import re
import os
import sys
import csv
import yaml

'''
This script contains a set of helper functions for various hichipper analysis calls
'''

def string_hamming_distance(str1, str2):
	'''
	Fast hamming distance over 2 strings known to be of same length.
	In information theory, the Hamming distance between two strings of equal 
	length is the number of positions at which the corresponding symbols 
	are different.
	eg "karolin" and "kathrin" is 3.
	'''
	return sum(itertools.imap(operator.ne, str1, str2))

def make_folder(folder):
	"""
	Function to only make a given folder if it does not already exist
	"""
	if not os.path.exists(folder):
		os.makedirs(folder)

def verify_file(filename):
	"""
	Ensure that file can both be read and exists
	"""
	try:
		extension = filename.split('.')[-1]
		if extension == "gz":
			fp = gzip.open(filename, 'rt')
		else:
			fp = open(filename)
		one = fp.readline()
		fp.close()
	except IOError as err:
		sys.exit(gettime() + "Error reading the file {0}: {1}".format(filename, err))
	return(filename)


def get_software_path(tool, abs_path):
	'''
	Function takes a tool name and a possible absolute path
	specified by the user input and returns the absolute path
	to the tool
	'''
	tool_path = which(tool)
	if(str(tool_path) == "None"):
		if(abs_path == ""):
			sys.exit("ERROR: cannot find "+tool+" in environment; add it to user PATH environment or specify executable using a flag.")
	if(abs_path != ""):
		if(os.path.isfile(abs_file)):
			tool_path = abs_path
	return(tool_path)

def rev_comp(seq):
	'''
	Fast Reverse Compliment
	'''  
	tbl = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
	return ''.join(tbl[s] for s in seq[::-1])

def gettime(): 
	'''
	Matches `date` in unix
	'''
	return(time.strftime("%a ") + time.strftime("%b ") + time.strftime("%d ") + time.strftime("%X ") + 
		time.strftime("%Z ") + time.strftime("%Y")+ ": ")
		

def findIdx(list1, list2):
	'''
	Return the indices of list1 in list2
	'''
	return [i for i, x in enumerate(list1) if x in list2]

def check_R_packages(required_packages, R_path):
	'''
	Determines whether or not R packages are properly installed
	'''
	installed_packages = os.popen(R_path + ''' -e "installed.packages()" | awk '{print $1}' | sort | uniq''').read().strip().split("\n")
	if(not set(required_packages) < set(installed_packages)):
		sys.exit("ERROR: cannot find the following R package: " + str(set(required_packages) - set(installed_packages)) + "\n" + 
			"Install it in your R console and then try rerunning proatac (but there may be other missing dependencies).")
		
def parse_manifest(manifest):
    samples = []
    if manifest.endswith(('.yaml', '.yml')):
        with open(manifest, 'r') as f: 
            m = yaml.load(f)
        return m
    else:
        click.echo(gettime() + "Please specify a valid .yaml file for samples")

# Simple intersection function
def intersect(a, b):
	return list(set(a) & set(b))

#Grab Sample Names
def get_subdirectories(dir):
	return [name for name in os.listdir(dir)
		if os.path.isdir(os.path.join(dir, name))]

# https://stackoverflow.com/questions/1006289/how-to-find-out-the-number-of-cpus-using-python	
def available_cpu_count():
	'''
	Number of available virtual or physical CPUs on this system, i.e.
	user/real as output by time(1) when called with an optimally scaling
	userspace-only program
	'''
	# cpuset
	# cpuset may restrict the number of *available* processors
	try:
		m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$',
					  open('/proc/self/status').read())
		if m:
			res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
			if res > 0:
				return res
	except IOError:
		pass

	# Python 2.6+
	try:
		import multiprocessing
		return multiprocessing.cpu_count()
	except (ImportError, NotImplementedError):
		pass

	# http://code.google.com/p/psutil/
	try:
		import psutil
		return psutil.cpu_count()   # psutil.NUM_CPUS on old versions
	except (ImportError, AttributeError):
		pass

	# POSIX
	try:
		res = int(os.sysconf('SC_NPROCESSORS_ONLN'))

		if res > 0:
			return res
	except (AttributeError, ValueError):
		pass

	# Windows
	try:
		res = int(os.environ['NUMBER_OF_PROCESSORS'])

		if res > 0:
			return res
	except (KeyError, ValueError):
		pass

	# jython
	try:
		from java.lang import Runtime
		runtime = Runtime.getRuntime()
		res = runtime.availableProcessors()
		if res > 0:
			return res
	except ImportError:
		pass

	# BSD
	try:
		sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'],
								  stdout=subprocess.PIPE)
		scStdout = sysctl.communicate()[0]
		res = int(scStdout)

		if res > 0:
			return res
	except (OSError, ValueError):
		pass

	# Linux
	try:
		res = open('/proc/cpuinfo').read().count('processor\t:')

		if res > 0:
			return res
	except IOError:
		pass

	# Solaris
	try:
		pseudoDevices = os.listdir('/devices/pseudo/')
		res = 0
		for pd in pseudoDevices:
			if re.match(r'^cpuid@[0-9]+$', pd):
				res += 1

		if res > 0:
			return res
	except OSError:
		pass

	# Other UNIXes (heuristic)
	try:
		try:
			dmesg = open('/var/run/dmesg.boot').read()
		except IOError:
			dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
			dmesg = dmesgProcess.communicate()[0]

		res = 0
		while '\ncpu' + str(res) + ':' in dmesg:
			res += 1

		if res > 0:
			return res
	except OSError:
		pass

	raise Exception('Can not determine number of CPUs on this system')

# Make sure that files are present
def verify_sample_hichipper_old(sample, boo, hicprooutput):
	bt = glob.glob(hicprooutput + '/bowtie_results/bwt2/' + sample + '/*.pairstat')
	hc1 = glob.glob(hicprooutput + '/hic_results/data/' + sample + '/*.RSstat')
	hc2 = glob.glob(hicprooutput + '/hic_results/data/' + sample + '/*.*Pairs')
	hc2 = [a for a in hc2 if re.search("DEPairs|SCPairs|validPairs|SinglePairs|DumpPairs",a) is not None]
	hc3 = glob.glob(hicprooutput + '/hic_results/data/' + sample + '/*[_,.]allValidPairs')
	
	if(boo):
		return((len(bt) > 0) and (len(hc1) > 0) and (len(hc2) % 5) == 0 and (len(hc3) > 0))
	
