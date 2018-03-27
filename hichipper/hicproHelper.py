import os
import click
from .hichipperHelp import *
from subprocess import call, check_call

#------------------------------------------------------------------
# Orient peaks based on any number of conditions; kind of dizzying
#------------------------------------------------------------------
def peakHelper(peaks, hicprooutput, resfrags, halfLength, peak_pad, out, samples,
	Rscript, skip_resfrag_pad, skip_background_correction,
	logf, macs2_string, macs2_genome, script_dir):
	
	peakfilespersample = []
	# Call peaks integrating over all samples for self ligation reads
	macs2Error = 'ERROR: macs2 peak calling was not successful; make sure macs2 is in the environment and that *Pairs files exist for samples. Finally, note that only the narrowPeak file is being considered here, which may be a problem if you tried to do broad peak calling.'
	if(peaks == "COMBINED,SELF"):
		click.echo(gettime() + "Calling one set of peaks from HiChIP self-ligation reads across all samples.", logf)
		for sample in samples:
			os.system('cat ' + hicprooutput + '/hic_results/data/' + sample + "/*DEPairs >> " + out + "/allDEandSCpairs.Pairs.tmp")
			os.system('cat ' + hicprooutput + '/hic_results/data/' + sample + "/*SCPairs >> " + out + "/allDEandSCpairs.Pairs.tmp")
		os.system('''awk -v RL=''' + str(halfLength)+r''' '{print $2"\t"$3-RL"\t"$3+RL}' '''+ out + '''/allDEandSCpairs.Pairs.tmp | awk '$2 > 0 {print $0}' > ''' + out + '''/allDEandSCpairs.bed.tmp''')
		os.system('''awk -v RL=''' + str(halfLength)+r''' '{print $5"\t"$6-RL"\t"$6+RL}' '''+ out + '''/allDEandSCpairs.Pairs.tmp | awk '$2 > 0 {print $0}' >> ''' + out + '''/allDEandSCpairs.bed.tmp''')
		macs2cmd = 'macs2 callpeak -t ' + out + '/allDEandSCpairs.bed.tmp --keep-dup all ' + macs2_string + " -g " + macs2_genome + " -B -f BED --verbose 0 -n " + out + '/allSamples_temporary'
		click.echo(gettime() +  "macs2 command: " + macs2cmd, logf)
		os.system(macs2cmd)
		if not os.path.isfile(out + '/allSamples_temporary_peaks.narrowPeak'):
			sys.exit(macs2Error)
	
		# Run hichipper modified background correction for all samples combined using self ligation reads
		if (not skip_background_correction):
			click.echo(gettime() + "Performing hichipper-modified background peak calling")
			call([Rscript, os.path.join(script_dir, 'lambdaProcess.R'), resfrags,  os.path.join(out, "allSamples_temporary_treat_pileup.bdg"), os.path.join(out,"allSamples_temporary_control_lambda.bdg"), out])
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
			macs2cmd = 'macs2 callpeak -t ' + out + '/' + sample + 'DEandSCpairs.bed.tmp  --keep-dup all ' + macs2_string + " -g " + macs2_genome + " -B -f BED --verbose 0 -n " + out + '/' + sample + "_temporary"
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
		macs2cmd = 'macs2 callpeak -t ' + out + '/allpairs.bed.tmp  --keep-dup all ' + macs2_string + " -g " + macs2_genome + " -B -f BED --verbose 0 -n " + out + '/allSamples_temporary'
		click.echo(gettime() +  "macs2 command: " + macs2cmd, logf)
		os.system(macs2cmd)
		if not os.path.isfile(out + '/allSamples_temporary_peaks.narrowPeak'):
			sys.exit(macs2Error)
		
		# Run hichipper modified background correction for all samples combined using all
		if (not skip_background_correction):
			click.echo(gettime() + "Performing hichipper-modified background peak calling")
			call([Rscript, os.path.join(script_dir, 'lambdaProcess.R'), resfrags,  os.path.join(out, "allSamples_temporary_treat_pileup.bdg"), os.path.join(out,"allSamples_temporary_control_lambda.bdg"), out])
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
			macs2cmd = 'macs2 callpeak -t ' + out + '/' + sample + '.all.Pairs.bed.tmp  --keep-dup all ' + macs2_string + " -g " + macs2_genome + " -B -f BED --verbose 0 -n " + out + '/' + sample + "_temporary"
			click.echo(gettime() +  "macs2 command: " + macs2cmd, logf)
			os.system(macs2cmd)
			if not os.path.isfile(out + '/' + sample + '_temporary_peaks.narrowPeak'):
				sys.exit(macs2Error)

			# Run hichipper modified background correction for each sample using self-ligation reads
			if (not skip_background_correction):
				click.echo(gettime() + "Performing hichipper-modified background peak calling")
				os.system(Rscript + ' ' + os.path.join(script_dir, 'lambdaProcess.R') + " " + resfrags + " " + out + '/' + sample + "_temporary_treat_pileup.bdg " + out + '/' + sample + "_temporary_control_lambda.bdg " + out)
				click.echo(gettime() + "Modified background signal; performing peak calling")
				click.echo(gettime() +  "MACS2 TEXT OUTPUT INCOMING")
				os.system("macs2 bdgcmp -t " + os.path.join(out,"adjustedTreatment.bdg.tmp") + " -c " + os.path.join(out, "adjustedBackground.bdg.tmp") + " -m qpois -o " + os.path.join(out,"hichipper_qvalue.bdg.tmp"))
				os.system("macs2 bdgpeakcall -i " + os.path.join(out,"hichipper_qvalue.bdg.tmp") + " -c 2 -l 147 -g 100 -o " + os.path.join(out,"hichipper_peaks.bed"))
				os.system("mv " + os.path.join(out,"hichipper_peaks.bed") + " " + out + "/" + sample + "_temporary_hichipperPeaks.bed")
				peakfilespersample.append(out + '/' + sample + '_temporary_hichipperPeaks.bed')
			else:
				peakfilespersample.append(out + '/' + sample + '_temporary_peaks.narrowPeak')
		

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
			call([Rscript, os.path.join(script_dir, 'resFragAnchorProcess.R'), resfrags, peakFile])
		peakfilespersample = [peakFile + "rf.tmp" for peakFile in peakfilespersample]
	return(peakfilespersample)
	
	
#------------------------------------------------------------------
# Determine samples from user-defined input for .yaml configuration
#------------------------------------------------------------------
def samplesHelper(keep_samples, ignore_samples, hicprooutput, logf):
	bwt_samples = os.popen('ls ' + hicprooutput + '/bowtie_results/bwt2').read().strip().split("\n")
	hic_samples = os.popen('ls ' + hicprooutput + '/hic_results/data').read().strip().split("\n")

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

	for sample in samples:
		if not verify_sample_hichipper_old(sample, True, hicprooutput):
			sys.exit('ERROR: Missing hic_results or bowtie_results files for ' + sample + "; either exclude the sample or manage the file architecture/input")
	return(samples)