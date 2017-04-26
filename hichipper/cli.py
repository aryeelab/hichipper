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
from pkg_resources import get_distribution
from subprocess import call, check_call

def parse_manifest(manifest):
    samples = []
    if manifest.endswith(('.yaml', '.yml')):
        with open(manifest, 'r') as f: 
            m = yaml.load(f)
        return m
    else:
        click.echo(gettime() + "Please specify a valid .yaml file for samples")

#Grab Sample Names
def get_subdirectories(dir):
    return [name for name in os.listdir(dir)
            if os.path.isdir(os.path.join(dir, name))]

@click.command()
@click.option('--out', default="hichipper_out", required=True, help='Output directory name; must not be already existing [Required]')
@click.option('--min-dist', default="5000", help='Minimum distance ; default = 5000')
@click.option('--max-dist', default="2000000", help='Peak padding width (applied on both left and right); default = 2000000')
@click.option('--macs2-string', default="-q 0.01 --extsize 147 --nomodel", help='String of arguments to pass to MACS2; only is called when peaks are set to be called; default = "-q 0.01 --extsize 147 --nomodel"')
@click.option('--macs2-genome', default="hs", help='Argument to pass to the -g variable in MACS2 (mm for mouse genome; hs for human genome); default = "hs"')
@click.option('--peak-pad', default="500", help='Peak padding width (applied on both left and right); default = 500')
@click.option('--merge-gap', default="500", help='Merge nearby peaks (after all padding is complete; default = 500')
@click.option('--keep-temp-files', is_flag=True, help='Keep temporary files?')
@click.option('--skip-background-correction', is_flag=True, help='Skip restriction fragment aware background correction?')
@click.option('--skip-resfrag-pad', is_flag=True, help='Skip restriction fragment aware padding')
@click.option('--skip-qc', is_flag=True, help='Skip QC report generation?')
@click.option('--skip-diffloop', is_flag=True, help='Skip analyses in diffloop (e.g. Mango loop calling; .rds generation)')
@click.option('--make-ucsc', is_flag=True, help='Make additional output files that can support viewing in UCSC genome browser; requires tabix and htslib tools.')
@click.option('--keep-samples', default="ALL", help='Comma separated list of sample names to keep; ALL (special string) by default')
@click.option('--ignore-samples', default="NONE", help='Comma separated list of sample names to ignore; NONE (special string) by default')
@click.option('--read-length', default="75", help='Length of reads from sequencing runs; default = 75')
@click.argument('manifest')
@click.version_option()

def main(manifest, out, min_dist, max_dist, macs2_string, macs2_genome, peak_pad, merge_gap, keep_temp_files, skip_background_correction, skip_resfrag_pad, skip_qc, skip_diffloop, read_length, keep_samples, ignore_samples, make_ucsc):
	"""A preprocessing and QC pipeline for HiChIP data."""
	__version__ = get_distribution('hichipper').version
	if os.path.exists(out):
		sys.exit("ERROR: Output path (%s) already exists." % out)
	
	def gettime(): # Matches `date` in Linux
		return(time.strftime("%a ") + time.strftime("%b ") + time.strftime("%d ") + time.strftime("%X ") + time.strftime("%Z ") + time.strftime("%Y")+ ": ")
	
	os.mkdir(out)
	logf = open(out + "/" + out + ".hichipper.log", 'w')
	click.echo(gettime() + "Starting hichipper pipeline v%s" % __version__, file = logf)
	click.echo(gettime() + "Starting hichipper pipeline v%s" % __version__)
	click.echo(gettime() + "Parsing user parameters")
	script_dir = os.path.dirname(os.path.realpath(__file__))
	outfolder = os.path.abspath(out)  
	cwd = os.getcwd()
	click.echo(gettime() + "Executed from: %s" % cwd, logf)
	click.echo(gettime() + "Output folder: %s" % outfolder, logf) 
	m = parse_manifest(manifest)
	click.echo(gettime() + "Parsed manifest as follows: ", logf)
	click.echo(m, logf)
	
	peaks = m['peaks'][0]
	resfrags = m['resfrags'][0]
	hicprooutput = m['hicpro_output'][0]
	
	# Determine if restriction fragment calling is happening or not
	if not (skip_resfrag_pad or skip_background_correction):
		if not os.path.isfile(resfrags):
			sys.exit('ERROR: Could not find the restriction fragment file ' + resfrags + '; either correctly specify file in .yaml or use the --skip-resfrag-pad and --skip-background-correction flags')
		
	# Determine samples from user-defined input
	bwt_samples = os.popen('ls ' + hicprooutput + '/bowtie_results/bwt2').read().strip().split("\n")
	hic_samples = os.popen('ls ' + hicprooutput + '/hic_results/data').read().strip().split("\n")
	
	def intersect(a, b):
		return list(set(a) & set(b))

	samples = intersect(bwt_samples, hic_samples)
	
	if(keep_samples != "ALL"):
		keeplist = keep_samples.split(",")
		click.echo(gettime() + "Intersecting detected samples with user-retained ones: " + keep_samples)
		samples = intersect(samples, keeplist)
		
	if(ignore_samples != "NONE"):
		igslist = ignore_samples.split(",")
		for byesample in igslist:
			click.echo(gettime() + "Attempting to remove " + byesample + " from processing", logf)
			if byesample in samples: samples.remove(byesample)
	
	if not len(samples) > 0:
		sys.exit('ERROR: Could not import any samples from the user specification; check flags, logs and .yaml configuration; quitting')
	
	halfLength = int(float(read_length)/2)
	if not halfLength> 0:
		sys.exit('ERROR: Specify appropriate read length > 1')
	
	
	# Make sure that files are present
	def verify_sample(sample, boo):
		bt = os.popen('ls ' + hicprooutput + '/bowtie_results/bwt2/' + sample + "/*.pairstat").read().strip().split("\n")
		hc1 = os.popen('ls ' + hicprooutput + '/hic_results/data/' + sample + "/*.RSstat").read().strip().split("\n")
		hc2 = os.popen('ls ' + hicprooutput + '/hic_results/data/' + sample + "/*.*Pairs").read().strip().split("\n")
		hc3 = os.popen('ls ' + hicprooutput + '/hic_results/data/' + sample + "/*_allValidPairs").read().strip().split("\n")
		
		if(boo):
			return((len(bt) > 0) and (len(hc1) > 0) and (len(hc2) % 5) == 0 and (len(hc3) > 0))
	
	for sample in samples:
		if not verify_sample(sample, True):
			sys.exit('ERROR: Missing hic_results or bowtie_results files for ' + sample + "; either exclude the sample or manage the file architecture/input")
	
	click.echo(gettime() + "Determined that the following samples are good to go: ", logf)
	click.echo(samples, logf)

	# Set up peaks files if necessary; moreover, specify vector of peaks / sample
	peakfilespersample = []
	click.echo(gettime() + "User defined peaks specification: " + peaks, logf)
	
	# Call peaks integrating over all samples for self ligation reads
	macs2Error = 'ERROR: macs2 peak calling was not successful; make sure macs2 is in the environment and that *Pairs files exist for samples. Finally, note that only the narrowPeak file is being considered here, which may be a problem if you tried to do broad peak calling.'
	if(peaks == "COMBINED,SELF"):
		click.echo(gettime() + "Calling one set of peaks from HiChIP self-ligation reads across all samples.", logf)
		for sample in samples:
		    os.system('cat ' + hicprooutput + '/hic_results/data/' + sample + "/*DEPairs >> " + out + "/allDEandSCpairs.Pairs.tmp")
		    os.system('cat ' + hicprooutput + '/hic_results/data/' + sample + "/*SCPairs >> " + out + "/allDEandSCpairs.Pairs.tmp")
		os.system('''awk -v RL=''' + str(halfLength)+r''' '{print $2"\t"$3-RL"\t"$3+RL}' '''+ out + '''/allDEandSCpairs.Pairs.tmp | awk '$2 > 0 {print $0}' > ''' + out + '''/allDEandSCpairs.bed.tmp''')
		os.system('''awk -v RL=''' + str(halfLength)+r''' '{print $5"\t"$6-RL"\t"$6+RL}' '''+ out + '''/allDEandSCpairs.Pairs.tmp | awk '$2 > 0 {print $0}' >> ''' + out + '''/allDEandSCpairs.bed.tmp''')
		macs2cmd = 'macs2 callpeak -t ' + out + '/allDEandSCpairs.bed.tmp ' + macs2_string + " -g " + macs2_genome + " -B -f BED --verbose 0 -n " + out + '/allSamples_temporary'
		click.echo(gettime() +  "macs2 command: " + macs2cmd, logf)
		os.system(macs2cmd)
		if not os.path.isfile(out + '/allSamples_temporary_peaks.narrowPeak'):
			sys.exit(macs2Error)
		
		# Run hichipper modified background correction for all samples combined using self ligation reads
		if (not skip_background_correction):
			click.echo(gettime() + "Performing hichipper-modified background peak calling")
			call(['Rscript', os.path.join(script_dir, 'lambdaProcess.R'), resfrags,  os.path.join(out, "allSamples_temporary_treat_pileup.bdg"), os.path.join(out,"allSamples_temporary_control_lambda.bdg"), out])
			click.echo(gettime() + "Modified background signal; performing peak calling")
			click.echo(gettime() +  "MACS2 TEXT OUTPUT INCOMING")
			os.system("macs2 bdgcmp -t " + os.path.join(out,"adjustedTreatment.bdg.tmp") + " -c " + os.path.join(out, "adjustedBackground.bdg.tmp") + " -m qpois -o " + os.path.join(out,"hichipper_qvalue.bdg.tmp"))
			os.system("macs2 bdgpeakcall -i " + os.path.join(out,"hichipper_qvalue.bdg.tmp") + " -c 2 -l 147 -g 100 -o " + os.path.join(out,"hichipper_peaks.bed"))
			os.system("mv " + os.path.join(out,"hichipper_peaks.bed") + " " + os.path.join(out,"allSamples_temporary_hichipperPeaks.bed"))
			peakfilespersample = [out + '/allSamples_temporary_hichipperPeaks.bed'] * len(samples)
		else:
			peakfilespersample = [out + '/allSamples_temporary_peaks.narrowPeak'] * len(samples)
	
	# Call peaks from within each sample using self-ligation reads
	elif(peaks == "EACH,SELF"):
		click.echo(gettime() + "Calling peaks from self-ligation reads for each sample", logf)
		for sample in samples:
		    os.system('cat ' + hicprooutput + '/hic_results/data/' + sample + "/*DEPairs >> " + out + "/" + sample +"DEandSCpairs.Pairs.tmp")
		    os.system('cat ' + hicprooutput + '/hic_results/data/' + sample + "/*SCPairs >> " + out + "/" + sample +"DEandSCpairs.Pairs.tmp")
		    os.system('''awk -v RL=''' + str(halfLength)+r''' '{print $2"\t"$3-RL"\t"$3+RL}' '''+ out + '''/''' + sample +'''DEandSCpairs.Pairs.tmp | awk '$2 > 0 {print $0}' >  ''' + out + '''/''' + sample +'''DEandSCpairs.bed.tmp''')
		    os.system('''awk -v RL=''' + str(halfLength)+r''' '{print $5"\t"$6-RL"\t"$6+RL}' '''+ out + '''/''' + sample +'''DEandSCpairs.Pairs.tmp | awk '$2 > 0 {print $0}' >> ''' + out + '''/''' + sample +'''DEandSCpairs.bed.tmp''')
		    macs2cmd = 'macs2 callpeak -t ' + out + '/' + sample + 'DEandSCpairs.bed.tmp ' + macs2_string + " -g " + macs2_genome + " -B -f BED --verbose 0 -n " + out + '/' + sample + "_temporary"
		    click.echo(gettime() +  "macs2 command: " + macs2cmd, logf)
		    os.system(macs2cmd)
		    if not os.path.isfile(out + '/' + sample + '_temporary_peaks.narrowPeak'):
				sys.exit(macs2Error)
			
		    # Run hichipper modified background correction for each sample using self-ligation reads
		    if (not skip_background_correction):
				click.echo(gettime() + "Performing hichipper-modified background peak calling")
				os.system('Rscript ' + os.path.join(script_dir, 'lambdaProcess.R') + " " + resfrags + " " + out + '/' + sample + "_temporary_treat_pileup.bdg " + out + '/' + sample + "_temporary_control_lambda.bdg " + out)
				click.echo(gettime() + "Modified background signal; performing peak calling")
				click.echo(gettime() +  "MACS2 TEXT OUTPUT INCOMING")
				os.system("macs2 bdgcmp -t " + os.path.join(out,"adjustedTreatment.bdg.tmp") + " -c " + os.path.join(out, "adjustedBackground.bdg.tmp") + " -m qpois -o " + os.path.join(out,"hichipper_qvalue.bdg.tmp"))
				os.system("macs2 bdgpeakcall -i " + os.path.join(out,"hichipper_qvalue.bdg.tmp") + " -c 2 -l 147 -g 100 -o " + os.path.join(out,"hichipper_peaks.bed"))
				os.system("mv " + os.path.join(out,"hichipper_peaks.bed") + " " + out + "/" + sample + "_temporary_hichipperPeaks.bed")
				peakfilespersample.append(out + '/' + sample + '_temporary_hichipperPeaks.bed')
		    else:
				peakfilespersample.append(out + '/' + sample + '_temporary_peaks.narrowPeak')
		    
	
	# Call peaks integrating over all samples using  all reads
	elif(peaks == "COMBINED,ALL"):
		click.echo(gettime() + "Calling one set of peaks from all HiChIP reads across all samples.", logf)
		for sample in samples:
		    os.system('cat ' + hicprooutput + '/hic_results/data/' + sample + "/*Pairs >> " + out + "/all.Pairs.tmp")
		os.system('''awk -v RL=''' + str(halfLength)+r''' '{print $2"\t"$3-RL"\t"$3+RL}' '''+ out + '''/all.Pairs.tmp | awk '$2 > 0 {print $0}' > ''' + out + '''/allpairs.bed.tmp''')
		os.system('''awk -v RL=''' + str(halfLength)+r''' '{print $5"\t"$6-RL"\t"$6+RL}' '''+ out + '''/all.Pairs.tmp | awk '$2 > 0 {print $0}' >> ''' + out + '''/allpairs.bed.tmp''')
		macs2cmd = 'macs2 callpeak -t ' + out + '/allpairs.bed.tmp ' + macs2_string + " -g " + macs2_genome + " -B -f BED --verbose 0 -n " + out + '/allSamples_temporary'
		click.echo(gettime() +  "macs2 command: " + macs2cmd, logf)
		os.system(macs2cmd)
		if not os.path.isfile(out + '/allSamples_temporary_peaks.narrowPeak'):
			sys.exit(macs2Error)
			
		# Run hichipper modified background correction for all samples combined using all
		if (not skip_background_correction):
			click.echo(gettime() + "Performing hichipper-modified background peak calling")
			call(['Rscript', os.path.join(script_dir, 'lambdaProcess.R'), resfrags,  os.path.join(out, "allSamples_temporary_treat_pileup.bdg"), os.path.join(out,"allSamples_temporary_control_lambda.bdg"), out])
			click.echo(gettime() + "Modified background signal; performing peak calling")
			click.echo(gettime() +  "MACS2 TEXT OUTPUT INCOMING")
			os.system("macs2 bdgcmp -t " + os.path.join(out,"adjustedTreatment.bdg.tmp") + " -c " + os.path.join(out, "adjustedBackground.bdg.tmp") + " -m qpois -o " + os.path.join(out,"hichipper_qvalue.bdg.tmp"))
			os.system("macs2 bdgpeakcall -i " + os.path.join(out,"hichipper_qvalue.bdg.tmp") + " -c 2 -l 147 -g 100 -o " + os.path.join(out,"hichipper_peaks.bed"))
			os.system("mv " + os.path.join(out,"hichipper_peaks.bed") + " " + os.path.join(out,"allSamples_temporary_hichipperPeaks.bed"))
			peakfilespersample = [out + '/allSamples_temporary_hichipperPeaks.bed'] * len(samples)
		else:
			peakfilespersample = [out + '/allSamples_temporary_peaks.narrowPeak'] * len(samples)
	
	# Call peaks from within each sample using self-ligation reads
	elif(peaks == "EACH,ALL"):
		click.echo(gettime() + "Calling peaks from all reads for each sample", logf)
		for sample in samples:
		    os.system('cat ' + hicprooutput + '/hic_results/data/' + sample + "/*Pairs >> " + out + "/" + sample +".all.Pairs.tmp")
		    os.system('''awk -v RL=''' + str(halfLength)+r''' '{print $2"\t"$3-RL"\t"$3+RL}' '''+ out + '''/''' + sample +'''.all.Pairs.tmp | awk '$2 > 0 {print $0}' >  ''' + out + '''/''' + sample +'''.all.Pairs.bed.tmp''')
		    os.system('''awk -v RL=''' + str(halfLength)+r''' '{print $5"\t"$6-RL"\t"$6+RL}' '''+ out + '''/''' + sample +'''.all.Pairs.tmp | awk '$2 > 0 {print $0}' >> ''' + out + '''/''' + sample +'''.all.Pairs.bed.tmp''')
		    macs2cmd = 'macs2 callpeak -t ' + out + '/' + sample + '.all.Pairs.bed.tmp ' + macs2_string + " -g " + macs2_genome + " -B -f BED --verbose 0 -n " + out + '/' + sample + "_temporary"
		    click.echo(gettime() +  "macs2 command: " + macs2cmd, logf)
		    os.system(macs2cmd)
		    if not os.path.isfile(out + '/' + sample + '_temporary_peaks.narrowPeak'):
				sys.exit(macs2Error)

		    # Run hichipper modified background correction for each sample using self-ligation reads
		    if (not skip_background_correction):
				click.echo(gettime() + "Performing hichipper-modified background peak calling")
				os.system('Rscript ' + os.path.join(script_dir, 'lambdaProcess.R') + " " + resfrags + " " + out + '/' + sample + "_temporary_treat_pileup.bdg " + out + '/' + sample + "_temporary_control_lambda.bdg " + out)
				click.echo(gettime() + "Modified background signal; performing peak calling")
				click.echo(gettime() +  "MACS2 TEXT OUTPUT INCOMING")
				os.system("macs2 bdgcmp -t " + os.path.join(out,"adjustedTreatment.bdg.tmp") + " -c " + os.path.join(out, "adjustedBackground.bdg.tmp") + " -m qpois -o " + os.path.join(out,"hichipper_qvalue.bdg.tmp"))
				os.system("macs2 bdgpeakcall -i " + os.path.join(out,"hichipper_qvalue.bdg.tmp") + " -c 2 -l 147 -g 100 -o " + os.path.join(out,"hichipper_peaks.bed"))
				os.system("mv " + os.path.join(out,"hichipper_peaks.bed") + " " + out + "/" + sample + "_temporary_hichipperPeaks.bed")
				peakfilespersample.append(out + '/' + sample + '_temporary_hichipperPeaks.bed')
		    else:
				peakfilespersample.append(out + '/' + sample + '_temporary_peaks.narrowPeak')
		    
	
	# Check to see if it's even a valid file
	elif not os.path.isfile(peaks):
		sys.exit('ERROR: Could not identify the ' + peaks + ' file; correctly specify file location or use special variable {COMBINED,EACH},{ALL,SELF} for peaks in .yaml')
	
	# Peaks are apparently determined or user-supplied...
	else:
		click.echo(gettime() + "Using user-defined peaks file " + peaks +  " for analysis.", logf)
		os.system("cp " + peaks + " " + out + "/userSuppliedPeaks.bed.tmp")
		peakfilespersample = [out + "/userSuppliedPeaks.bed.tmp"] * len(samples)
	
	# Do the initial peak padding
	if(int(peak_pad) > 0):
		for peakFile in list(set(peakfilespersample)):
			os.system('''awk -v RL=''' + peak_pad + r''' '{print $1"\t"$2-RL"\t"$3+RL}' '''+ peakFile + ''' | awk '$2 > 0 {print $0}' >  ''' + peakFile + '''_pad.bed.tmp''')
		peakfilespersample = [peakFile + "_pad.bed.tmp" for peakFile in peakfilespersample]
	else:
		click.echo(gettime() + "No initial padding as --peak-pad set to 0 or less", logf)
		
	# Restriction fragment aware padding	
	if skip_resfrag_pad:
		click.echo(gettime() + "Skipping restriction fragment-aware padding based on user flag", logf)
	else:
		click.echo(gettime() + "Performing restriction fragment-aware padding", logf)
		for peakFile in list(set(peakfilespersample)):
			call(['Rscript', os.path.join(script_dir, 'resFragAnchorProcess.R'), resfrags, peakFile])
		peakfilespersample = [peakFile + "rf.tmp" for peakFile in peakfilespersample]
 	logf.close()
 	
 	# Check for UCSC Specification
 	if make_ucsc:
 		ucscoutput = "true"
 	else:
 		ucscoutput = "false"
 	
 	
 	# Call putative interactions
	for i in range(len(samples)):
		hichipperRun = os.path.join(script_dir, 'interactionsCall.sh')	
		cmd = ['bash', hichipperRun, cwd, out, hicprooutput, samples[i], peakfilespersample[i], min_dist, max_dist, merge_gap, str(halfLength), ucscoutput]        
		call(cmd)
		if not os.path.isfile(out + "/" + samples[i] + ".stat"):
			sys.exit('ERROR: something failed at the individual sample level; check the .log file for more info')
	
	# Back to Python
	logf = open(out + "/" + out + ".hichipper.log", 'a')
		
	# QC Report			
	if skip_qc:
		click.echo(gettime() + "Skipping QC report generation since --skip-qc was specified", logf)
	else:
		click.echo(gettime() + "Creating QC report", logf)
		cftg = ' '.join(samples)
		call(['cp', os.path.join(script_dir, 'qcReport_make.Rmd'), out + '/qcReport_make.Rmd'])
		call(['cp', os.path.join(script_dir, 'qcReport.R'), out + '/qcReport.R'])
		cmd = ['Rscript', out + '/qcReport.R', script_dir, str(out), cwd, __version__ , cftg] 		
		click.echo(gettime())
		click.echo(cmd)
		call(cmd)
		
	# diffloop work
	if skip_diffloop:
		click.echo(gettime() + "Skipping diffloop analyses since --skip-diffloop was specified", logf)
	else:
		click.echo(gettime() + "Creating .rds and .mango files", logf)
		cftg = ' '.join(samples)
		cmd = ['Rscript', os.path.join(script_dir, 'diffloop_work.R'), cwd, out, cftg] 
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