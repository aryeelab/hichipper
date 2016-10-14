import click
import os
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
        sample_names = m['samples'].keys()
        sample_names.sort()
        for sample_name in sample_names:
            runs = m['samples'][sample_name]
            read1 = []
            read2 = []
            for run in runs:
                bam1, bam2 = run.split(" ")
                read1.append(bam1)
                read2.append(bam2)
            d = {
                'name': sample_name, 
                'read1': ','.join(read1),
                'read2': ','.join(read2)
                }
            samples.append(d)
        return samples
    else:
        click.echo("Please specify a valid .yaml file for samples")

#Grab Sample Names
def get_subdirectories(dir):
    return [name for name in os.listdir(dir)
            if os.path.isdir(os.path.join(dir, name))]

@click.command()
@click.option('--out', default=".", required=True, help='Output directory name')
@click.option('--peak-pad', default="1500", help='Peak padding width (applied on both left and right); default = 1500')
@click.option('--merge-gap', default="1500", help='Max gap size for merging peaks; default = 1500')
@click.option('--keep-temp-files', is_flag=True, help='Keep temporary files?')
@click.option('--skip-R', is_flag=True, help='Skip QC report generation? (Requires R)')
@click.argument('manifest')

def main(manifest, out, peak_pad, merge_gap, keep_temp_files, skip_R):
    """A preprocessing and QC pipeline for HiChIP data."""
    __version__ = get_distribution('dnaloop').version
    click.echo("Starting hichipper pipeline v%s" % __version__)
    if os.path.exists(out):
        sys.exit("ERROR: Output path (%s) already exists." % out)
    os.mkdir(out)        
    os.mkdir(os.path.join(out, 'log'))
    with open(os.path.join(out, 'log', 'VERSION.txt'), 'w') as f: 
        f.write(__version__ + '\n')
    script_dir = os.path.dirname(os.path.realpath(__file__))
    out = os.path.abspath(out)    
    click.echo("Output folder: %s" % out) 
    # Preprocess individual samples
    samples = parse_manifest(manifest)
    i = 0
    for sample in samples:
        i += 1
        click.echo("\nProcessing sample %d of %d: %s" % (i, len(samples), sample['name']))
        click.echo("    Read 1: %s" % sample['read1']) 
        click.echo("    Read 2: %s" % sample['read2'])    
        hichipperRun = os.path.join(script_dir, 'hichipper.sh')
        if allow_missing_linker:
            allow_missing_linker_opt = "allow_missing_linker=yes"
        else:
            allow_missing_linker_opt = "allow_missing_linker=no"
        
        cmd = [hichipperRun, os.path.join(out, 'samples', sample['name']), merge_gap, sample['read1'], sample['read2'], bwa_mode, allow_missing_linker_opt, linkers]
        click.echo("    Executing: %s" % " ".join(cmd))
        call(cmd)

    
    if skip_R:
        click.echo("Skipping QC report generation since --skip-R was specified")
    else:
        click.echo("Creating QC report")
        cmd = ['Rscript', os.path.join(script_dir, 'qcReport.R'), out] 
        call(cmd)
    if keep_temp_files:
        click.echo("Temporary files not deleted since --keep-temp-files was specified")
    else:
        click.echo("Deleting temporary files")
    click.echo("Done")



