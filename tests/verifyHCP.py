import glob
import sys
import re
import os

from optparse import OptionParser

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Script to verify HiC-Pro architecture for a given sample"
opts = OptionParser(usage=usage)
opts.add_option("--hicpro", help="Path to HiC-Pro directory")
opts.add_option("--sample", help="sample name within the HiC-Pro directory")

options, arguments = opts.parse_args()

hicprooutput = options.hicpro
sample = options.sample

bt = glob.glob(hicprooutput + '/bowtie_results/bwt2/' + sample + '/*.pairstat')
hc1 = glob.glob(hicprooutput + '/hic_results/data/' + sample + '/*.RSstat')
hc2 = glob.glob(hicprooutput + '/hic_results/data/' + sample + '/*.*Pairs')
hc2 = [a for a in hc2 if re.search("DEPairs|SCPairs|validPairs|SinglePairs|DumpPairs",a) is not None]
hc3 = glob.glob(hicprooutput + '/hic_results/data/' + sample + '/*_allValidPairs')

print(bt)
print(hc1)
print(hc2)
print(hc3)

print(len(bt))
print(len(hc1))
print(len(hc2))
print(len(hc3))