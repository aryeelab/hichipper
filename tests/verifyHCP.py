import sys
import re
import os

from optparse import OptionParser

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Script to verify HiC-Pro architecture for a given sample"
opts = OptionParser(usage=usage)
opts.add_option("--hicpro", help="Path to HiC-Pro directory")
opts.add_option("--sample", help="Path to HiC-Pro directory")

options, arguments = opts.parse_args()

hicprooutput = options.hicpro
sample = options.sample

bt = os.popen('ls ' + hicprooutput + '/bowtie_results/bwt2/' + sample + "/*.pairstat").read().strip().split("\n")
hc1 = os.popen('ls ' + hicprooutput + '/hic_results/data/' + sample + "/*.RSstat").read().strip().split("\n")
hc2 = os.popen('ls ' + hicprooutput + '/hic_results/data/' + sample + "/*.*Pairs").read().strip().split("\n")
hc2 = [a for a in hc2 if re.search("DEPairs|SCPairs|validPairs|SinglePairs|DumpPairs",a) is not None]
hc3 = os.popen('ls ' + hicprooutput + '/hic_results/data/' + sample + "/*_allValidPairs").read().strip().split("\n")

print(bt)
print(hc1)
print(hc2)
print(hc3)

print(len(bt))
print(len(hc1))
print(len(hc2))
print(len(hc3))