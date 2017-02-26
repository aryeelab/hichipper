import click
import os
import os.path
import sys
import shutil
import yaml
import shutil
import random
import string
from pkg_resources import get_distribution
from subprocess import call, check_call

def parse_manifest(manifest):
    samples = []
    if manifest.endswith(('.yaml', '.yml')):
        with open(manifest, 'r') as f: 
            m = yaml.load(f)
        return m
    else:
        click.echo("Please specify a valid .yaml file for samples")

#Grab Sample Names
def get_subdirectories(dir):
    return [name for name in os.listdir(dir)
            if os.path.isdir(os.path.join(dir, name))]

@click.command()
@click.option('--out', default="hichipper_out", required=True, help='Output directory name')
@click.option('--min-dist', default="5000", help='Minimum distance ; default = 5000')
@click.option('--max-dist', default="2000000", help='Peak padding width (applied on both left and right); default = 2000000')
@click.option('--macs2-string', default="-q 0.01 --extsize 147 --nomodel", help='String of arguments to pass to MACS2; only is called when peaks in .yaml is set to NONE; default = "-q 0.01 --extsize 147 --nomodel"')
@click.option('--peak-pad', default="1000", help='Peak padding width (applied on both left and right); default = 1500')
@click.option('--keep-temp-files', is_flag=True, help='Keep temporary files?')
@click.option('--skip-resfrag', is_flag=True, help='Skip restriction fragment padding?')
@click.option('--skip-qc', is_flag=True, help='Skip QC report generation?')
@click.option('--skip-diffloop', is_flag=True, help='Skip QC report generation?')
@click.option('--ignore-samples', default="NONE", help='Comma separated list of sample names to ignore')
@click.option('--read-length', default="75", help='Length of reads from sequencing runs')
@click.argument('manifest')
@click.version_option()

def main(manifest, out, min_dist, max_dist, macs2_string, peak_pad, keep_temp_files, skip_resfrag, skip_qc, skip_diffloop, read_length, ignore_samples):
	"""A preprocessing and QC pipeline for HiChIP data."""
	__version__ = get_distribution('hichipper').version
	
	click.echo("Starting hichipper pipeline v%s" % __version__)
	if os.path.exists(out):
		sys.exit("ERROR: Output path (%s) already exists." % out)
	os.mkdir(out)        
	script_dir = os.path.dirname(os.path.realpath(__file__))
	outfolder = os.path.abspath(out)  
	cwd = os.getcwd()
	click.echo("Executed from: %s" % cwd)
	click.echo("Output folder: %s" % outfolder) 
	m = []
	m = parse_manifest(manifest)
	click.echo(m)
	
	peaks = m['peaks'][0]
	resfrags = m['resfrags'][0]
	hicprooutput = m['hicpro_output'][0]
	
	# Determine if restriction fragment calling is happening or not
	if not skip_resfrag:
		if not os.path.isfile(resfrags):
			sys.exit('ERROR: Could not find the restriction fragment file ' + resfrags + '; either correctly specify file in .yaml or use the --skip-resfrag flag')
		
	# Determine samples from user-defined input
	bwt_samples = os.popen('ls ' + hicprooutput + '/bowtie_results/bwt2').read().strip().split("\n")
	hic_samples = os.popen('ls ' + hicprooutput + '/hic_results/data').read().strip().split("\n")
	
	def intersect(a, b):
		return list(set(a) & set(b))

	samples = intersect(bwt_samples, hic_samples)
	
	if(ignore_samples != "NONE"):
		igslist = ignore_samples.split(",")
		for byesample in igslist:
			click.echo("Attempting to remove " + byesample + " from processing")
			if byesample in samples: samples.remove(byesample)
	
	if not len(samples) > 0:
		sys.exit('ERROR: Could not import any samples from the .yaml specification; quitting')
	
	halfLength = int(float(read_length)/2)
	if not halfLength> 0:
		sys.exit('ERROR: Specify appropriate read length > 1')
	
	
	# Make sure that files are present
	def verify_sample(sample, boo):
		bt = os.popen('ls ' + hicprooutput + '/bowtie_results/bwt2/' + sample + "/*.pairstat").read().strip().split("\n")
		hc1 = os.popen('ls ' + hicprooutput + '/hic_results/data/' + sample + "/*.RSstat").read().strip().split("\n")
		hc2 = os.popen('ls ' + hicprooutput + '/hic_results/data/' + sample + "/*Pairs").read().strip().split("\n")
		if(boo):
			return((len(bt) > 0) and (len(hc1) > 0) and (len(hc2) % 5) == 0)
	
	for sample in samples:
		if not verify_sample(sample, True):
			sys.exit('ERROR: Missing hic_results or bowtie_results files for ' + sample + "; either exclude the sample or manage the file architecture/input")
	
	click.echo("Determined that the following samples are good to go: ")
	click.echo(samples)

	# Set up peaks files if necessary; moreover, specify vector of peaks / sample
	peakfilespersample = []
	click.echo("User defined peaks specification: " + peaks)
	
	# Call peaks integrating over all samples
	macs2Error = 'ERROR: macs2 peak calling was not successful; make sure macs2 is in the environment and that *DEPairs/*SCPairs files exist for samples. Finally, note that only the narrowPeak file is being considered here, which may be a problem if you tried to do broad peak calling.'
	if(peaks == "ALL"):
		click.echo("Calling one set of peaks from HiChIP self-ligation reads across all samples.")
		for sample in samples:
		    os.system('cat ' + hicprooutput + '/hic_results/data/' + sample + "/*DEPairs >> " + out + "/allDEandSCpairs.Pairs")
		    os.system('cat ' + hicprooutput + '/hic_results/data/' + sample + "/*SCPairs >> " + out + "/allDEandSCpairs.Pairs")
		os.system('''awk -v RL=''' + str(halfLength)+r''' '{print $2"\t"$3-RL"\t"$3+RL}' '''+ out + '''/allDEandSCpairs.Pairs | awk '$2 > 0 {print $0}' > ''' + out + '''/allDEandSCpairs.bed''')
		os.system('''awk -v RL=''' + str(halfLength)+r''' '{print $5"\t"$6-RL"\t"$6+RL}' '''+ out + '''/allDEandSCpairs.Pairs | awk '$2 > 0 {print $0}' >> ''' + out + '''/allDEandSCpairs.bed''')
		macs2cmd = 'macs2 callpeak -t ' + out + '/allDEandSCpairs.bed ' + macs2_string + " -f BED -n " + out + '/allSamples'
		click.echo( "macs2 command: " + macs2cmd)
		os.system(macs2cmd)
		if not os.path.isfile(out + '/allSamples_peaks.narrowPeak'):
			sys.exit(macs2Error)
		peakfilespersample = [out + '/allSamples_peaks.narrowPeak'] * len(samples)
	
	# Call peaks from within each sample
	elif(peaks == "EACH"):
		click.echo("Calling peaks from self-ligation reads for each sample")
		for sample in samples:
		    os.system('cat ' + hicprooutput + '/hic_results/data/' + sample + "/*DEPairs >> " + out + "/" + sample +"DEandSCpairs.Pairs")
		    os.system('cat ' + hicprooutput + '/hic_results/data/' + sample + "/*SCPairs >> " + out + "/" + sample +"DEandSCpairs.Pairs")
		    os.system('''awk -v RL=''' + str(halfLength)+r''' '{print $2"\t"$3-RL"\t"$3+RL}' '''+ out + '''/''' + sample +'''DEandSCpairs.Pairs | awk '$2 > 0 {print $0}' >  ''' + out + '''/''' + sample +'''DEandSCpairs.bed''')
		    os.system('''awk -v RL=''' + str(halfLength)+r''' '{print $5"\t"$6-RL"\t"$6+RL}' '''+ out + '''/''' + sample +'''DEandSCpairs.Pairs | awk '$2 > 0 {print $0}' >> ''' + out + '''/''' + sample +'''DEandSCpairs.bed''')
		    macs2cmd = 'macs2 callpeak -t ' + out + '/' + sample + 'DEandSCpairs.bed ' + macs2_string + " -f BED -n " + out + '/' + sample
		    click.echo( "macs2 command: " + macs2cmd)
		    os.system(macs2cmd)
		    if not os.path.isfile(out + '/' + sample + '_peaks.narrowPeak'):
				sys.exit(macs2Error)
		    peakfilespersample.append(out + '/' + sample + '_peaks.narrowPeak')
	
	# Check to see if it's even a valid file
	elif not os.path.isfile(peaks):
		sys.exit('ERROR: Could not identify the ' + peaks + ' file; correctly specify file location or use special variables ALL or EACH in .yaml')
	
	# Peaks are apparently determined or user-supplied...
	else:
		click.echo("Using user-defined peaks file " + peaks +  " for analysis.")
		peakfilespersample = [peaks] * len(samples)
	
	# Do the initial peak padding
	if(int(peak_pad) > 0):
		for peakFile in list(set(peakfilespersample)):
			os.system('''awk -v RL=''' + peak_pad + r''' '{print $1"\t"$2-RL"\t"$2+RL}' '''+ peakFile + ''' | awk '$2 > 0 {print $0}' >  ''' + peakFile + '''_pad.bed''')
		peakfilespersample = [peakFile + "_pad.bed" for peakFile in peakfilespersample]
	else:
		click.echo("No initial padding as --peak-pad set to 0 or less")
		
	# Restriction fragment aware padding	
	if skip_resfrag:
		click.echo("Skipping restriction fragment-aware padding based on user flag")
	else:
		click.echo("Performing restriction fragment-aware padding")
		for peakFile in list(set(peakfilespersample)):
			call(['Rscript', os.path.join(script_dir, 'resFragAnchorProcess.R'), resfrags, peakFile])
		peakfilespersample = [peakFile + "rf" for peakFile in peakfilespersample]
 
	for i in range(len(samples)):
		hichipperRun = os.path.join(script_dir, 'hichipper.sh')	
		cmd = ['bash', hichipperRun, cwd, out, sample['name'], sample['read1'], sample['read2'], peak_pad, merge_gap, min_qual, min_dist, max_dist, macs2_string]        
	call(cmd)
	statout = out + "/" + sample['name'] + ".stat"
	
	if not os.path.isfile(statout):
		sys.exit('ERROR: Could not execute the pipeline successfully; check the .log file for more info')
		
	# QC Report			
	if skip_qc:
		click.echo("Skipping QC report generation since --skip-qc was specified")
	else:
		click.echo("Creating QC report")
		cftg = ' '.join(sample_names)
		cmd = ['Rscript', os.path.join(script_dir, 'qcReport.R'), script_dir, out, cwd, __version__ , cftg] 
		#click.echo("    Executing: %s" % " ".join(cmd))
		call(cmd)
		
	# diffloop work
	if skip_diffloop:
		click.echo("Skipping diffloop analyses since --skip-diffloop was specified")
	else:
		click.echo("Creating .rds and .mango files")
		cftg = ' '.join(sample_names)
		cmd = ['Rscript', os.path.join(script_dir, 'diffloop_work.R'), cwd, out, cftg] 
		#click.echo("    Executing: %s" % " ".join(cmd))
		call(cmd)	
	
	# Temporary File Management
	if keep_temp_files:
		click.echo("Temporary files not deleted since --keep-temp-files was specified")
	else:
		click.echo("Deleting temporary files")
		files = os.listdir(outfolder)
		for file in files:
			if file.endswith(".tmp"):
				os.remove(os.path.join(outfolder,file))
			elif "_temporary_" in file:
				os.remove(os.path.join(outfolder,file))
	click.echo("Done")