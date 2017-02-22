# hichipper
This package is maintained by [Caleb Lareau](mailto:caleblareau@g.harvard.edu) in the [Aryee Lab](http://aryee.mgh.harvard.edu/). Source code is made freely available here and a packaged install version is provided through [PyPi](https://pypi.python.org/pypi/hichipper/).

[![Build Status](https://travis-ci.org/aryeelab/hichipper.svg?branch=master)](https://travis-ci.org/aryeelab/hichipper) [![PyPI version](https://badge.fury.io/py/hichipper.svg)](https://badge.fury.io/py/hichipper) [![MIT Licence](https://badges.frapsoft.com/os/mit/mit.svg?v=103)](https://opensource.org/licenses/mit-license.php) 

## About<a name="about"></a>

The **hichipper** package implements our data processing and quality control pipeline for 
[HiChIP](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3999.html) data.
This package takes output from a [HiC-Pro](https://github.com/nservant/HiC-Pro) run and a sample manifest file (`.yaml`) that coordinates optional high-quality peaks (identified through ChIP-Seq) and restriction fragment locations (see [folder here](RestrictionFragmentFiles))
as input and produces output that can be used to 1) determine library quality, 2) identify and characterize DNA loops and 3) interactively visualize loops. Loops are assigned strength and confidence metrics that can be used to evaluate samples individually or for differential analysis in downstream tools like [diffloop](http://github.com/aryeelab/diffloop). We have used the library QC metrics with as few as 1 million reads, enabling library quality to be assessed through shallow (and cheap) sequencing before performing a full depth sequencing run.

A graphical overview showing how **hichipper** integrates with other tools in the analysis of raw HiChIP data is shown in the overview figure below. Detailed descriptions of the different branches of output from **hichipper** are discussed at the bottom of this guide. 
![big1](media/Big1.png)
![big2](media/Big2.png)<br>
A higher resolution [slide of this image](media/Big.pptx) is in the [media](media) folder.


## Table of Contents<a name="toc"></a>
- [About](#about)
- [Table of Contents](#toc)
- [Workflow Overview](#ugo)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Simple Usage Example](#sue)
- [More typical example](#moe)
- [Output](#output)
- [Configurations](#configuration)
- [Parameter Explanations](#pe)
- [Quality Control reports](#qcr)
- [Interactive visualization of loops](#viz)
- [Analyzing loops in the R](#loops)

## Workflow Overview<a name="ugo"></a>
A simple graphical guide to processing HiChIP data is shown below. The role of **hichipper**
is to import aligned read files from (e.g. [HiC-Pro](https://github.com/nservant/HiC-Pro))
as well as a sample `.yaml` file and produce user-friendly output. 
 
![hichipper_overview](media/Overview.png)


## Dependencies<a name="dependencies"></a>

The following dependencies need to be installed before running **hichipper**: [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html), OpenSSL, libcurl, and libxml2

On an Ubuntu system these can be installed with:
```
apt-get install bedtools libssl-dev libcurl4-openssl-dev libxml2-dev
```

Additionally, `R` must be available in the environment as well as [pandoc](http://pandoc.org) and a few packages that can be downloaded running the following in an 'R' environment:

```
install_pkgs <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
}
install_pkgs(c("DT", "devtools", "foreach", "ggplot2", "knitr", "networkD3", "readr", "reshape2"))
devtools::install_github("aryeelab/diffloop")
```

Convenient [pandoc binaries](https://s3.amazonaws.com/rstudio-buildtools/pandoc-1.12.4.2.zip) for Linux, Mac and Windows are available for download from RStudio.

## Installation<a name="installation"></a>

To install **hichipper** given the dependencies above, simply run:

```
pip install hichipper
```

## Simple usage example<a name="sue"></a>

The example below uses the test dataset bundled with the **hichipper** package source code. Download the package and change to the test directory with:

```
git clone https://github.com/aryeelab/hichipper.git
cd hichipper/tests
```

1. Create a sample description file:
  
  Sample description files can be created with the `.yaml` format. 

  **Processing `.yaml` format**
   
   Example [yaml](https://en.wikipedia.org/wiki/YAML) format sample description file:
   
   ```
   FIX ME
   ```
  Note: This file is available as `example.yaml` in the `hichipper/tests` directory.
  
  In this example, we call loops from two GM12878 samples using just chromosome 22.  
  
  
2. Run the pipeline:
	```
	hichipper --out output1 example.yaml
	```

## More typical example<a name="moe"></a>
While the example above references files that are part of the **hichipper** distribution,
our experience using this tool in conjunction with HiC-Pro suggests that a file hierarchy
like the following may be more typical. 

```
ls -LR
```
yields (a slightly modified version of)
```
<pre>

./hicpro/bowtie_results:
bwt2  bwt2_global  bwt2_local

./hicpro/bowtie_results/bwt2:
SRR3467175  SRR3467176  SRR3467177  SRR3467178

./hicpro/bowtie_results/bwt2/SRR3467175:
SRR3467175_1_hg19.bwt2merged.bam  SRR3467175_2_hg19.bwt2merged.bam  SRR3467175_hg19.bwt2pairs.bam
SRR3467175_1_hg19.mapstat         SRR3467175_2_hg19.mapstat         <b>SRR3467175_hg19.bwt2pairs.pairstat</b>

./hicpro/bowtie_results/bwt2/SRR3467176:
SRR3467176_1_hg19.bwt2merged.bam  SRR3467176_2_hg19.bwt2merged.bam  SRR3467176_hg19.bwt2pairs.bam
SRR3467176_1_hg19.mapstat         SRR3467176_2_hg19.mapstat         <b>SRR3467176_hg19.bwt2pairs.pairstat</b>

./hicpro/bowtie_results/bwt2/SRR3467177:
SRR3467177_1_hg19.bwt2merged.bam  SRR3467177_2_hg19.bwt2merged.bam  SRR3467177_hg19.bwt2pairs.bam
SRR3467177_1_hg19.mapstat         SRR3467177_2_hg19.mapstat         <b>SRR3467177_hg19.bwt2pairs.pairstat</b>

./hicpro/bowtie_results/bwt2/SRR3467178:
SRR3467178_1_hg19.bwt2merged.bam  SRR3467178_2_hg19.bwt2merged.bam  SRR3467178_hg19.bwt2pairs.bam
SRR3467178_1_hg19.mapstat         SRR3467178_2_hg19.mapstat         <b>SRR3467178_hg19.bwt2pairs.pairstat</b>

./hicpro/hic_results:
data

./hicpro/hic_results/data:
SRR3467175  SRR3467176  SRR3467177  SRR3467178

./hicpro/hic_results/data/SRR3467175:
<b>SRR3467175_hg19.bwt2pairs.DEPairs    SRR3467175_hg19.bwt2pairs.RSstat   SRR3467175_hg19.bwt2pairs.SinglePairs
SRR3467175_hg19.bwt2pairs.DumpPairs  SRR3467175_hg19.bwt2pairs.SCPairs  SRR3467175_hg19.bwt2pairs.validPairs</b>

./hicpro/hic_results/data/SRR3467176:
<b>SRR3467176_hg19.bwt2pairs.DEPairs    SRR3467176_hg19.bwt2pairs.RSstat   SRR3467176_hg19.bwt2pairs.SinglePairs
SRR3467176_hg19.bwt2pairs.DumpPairs  SRR3467176_hg19.bwt2pairs.SCPairs  SRR3467176_hg19.bwt2pairs.validPairs</b>

./hicpro/hic_results/data/SRR3467177:
<b>SRR3467177_hg19.bwt2pairs.DEPairs    SRR3467177_hg19.bwt2pairs.RSstat   SRR3467177_hg19.bwt2pairs.SinglePairs
SRR3467177_hg19.bwt2pairs.DumpPairs  SRR3467177_hg19.bwt2pairs.SCPairs  SRR3467177_hg19.bwt2pairs.validPairs</b>

./hicpro/hic_results/data/SRR3467178:
<b>SRR3467178_hg19.bwt2pairs.DEPairs    SRR3467178_hg19.bwt2pairs.RSstat   SRR3467178_hg19.bwt2pairs.SinglePairs
SRR3467178_hg19.bwt2pairs.DumpPairs  SRR3467178_hg19.bwt2pairs.SCPairs  SRR3467178_hg19.bwt2pairs.validPairs</b>

...
</pre>

```
where files denoted in **bold** are assumed to exist. Typically, an analysis folder may look like so:

```
fastq/
|-- S75
|  |-- SRR3467175_1.fastq.gz
|  |-- SRR3467175_2.fastq.gz
|-- S76
|  |-- SRR3467176_1.fastq.gz
|  |-- SRR3467176_2.fastq.gz
|-- S77
|  |-- SRR3467177_1.fastq.gz
|  |-- SRR3467177_2.fastq.gz
|-- S78
|  |-- SRR3467178_1.fastq.gz
|  |-- SRR3467178_2.fastq.gz
hicpro/
|-- HiCPro_step1_hic.sh
|-- bowtie_results/
|  |-- bwt2/
|  |  |-- SRR3467175
|  |  |  |-- SRR3467175_hg19.bwt2pairs.pairstat
|  |  |-- SRR3467176
|  |  |  |-- SRR3467176_hg19.bwt2pairs.pairstat
|  |  |-- SRR3467177
|  |  |  |-- SRR3467177_hg19.bwt2pairs.pairstat
|  |  |-- SRR3467178
|  |  |  |-- SRR3467178_hg19.bwt2pairs.pairstat
|-- hic_results/
|  |-- data/
|  |  |-- SRR3467175
|  |  |  |-- SRR3467175*Pairs # 6 Files
|  |  |-- SRR3467176
|  |  |  |-- SRR3467176*Pairs # 6 Files
|  |  |-- SRR3467177
|  |  |  |-- SRR3467177*Pairs # 6 Files
|  |  |-- SRR3467178
|  |  |  |-- SRR3467178*Pairs # 6 Files
GM12878_SMC3_ChIPSeq.narrowPeak
hg19_MboI_resfrag.bed.gz
config.yaml
config-hicpro-mboi-ext12.txt
```
where the results in the `hicpro` directory could have been obtained by running: 
```
HiC-Pro -i fastq/ -o hicpro/ -c config-hicpro-mboi-ext12.txt -p
```
and subsequently executing the resulting `HiCPro_step1_hic.sh`.  Thus, the `config.yaml` file
needed for **hichipper** when executed from the current working directory would look like this:

```
peaks:
  - GM12878_SMC3_ChIPSeq.narrowPeak
resfrags:
  - hg19_MboI_resfrag.bed.gz
hicpro_output:
  - hicpro
```

And could be executed running this command: 

```
hichipper --out GM12878 config.yaml
```

would yield the default output from **hichipper**. 

## Output<a name="output"></a>

### Per-run output files
Each time the user runs **hichipper**, a `*.hichipper.log` file containing information pertaining to the 
flow of the software execution is placed in the `out` directory (specified by the `--out` flag). Unless
otherwise specified, a file ending in `hichipper-qcReport.html` provides an interactive quality control report for all samples. 

### Per-sample output files
Per sample, six (yes, 6, but don't worry-- there's lots of redundancy) output files are created. They are:

1. `*.stat` Key summary statistics that show the number of PETs meeting certain criteria

2. `*.inter.loop_counts.bedpe` Interchromosomal looping between anchor loci. 

3. `*.intra.loop_counts.bedpe` Intrachromosomal looping between **all**  anchor loci

4. `*.filt.intra.loop_counts.bedpe` Intrachromosomal looping between anchor loci where loops meet min/max distance requirements.

5. `*interactions.all.mango` The same set of loops as 4 but with per-loop FDR measures from the loop proximity bias correction algorithm originally implemented in [Mango](https://github.com/dphansti/mango) and presented in the same format. 

6. `*.rds` The same set of loops as 4 but in an R binary compressed format of a `loops()` S4 object from [diffloop](http://bioconductor.org/packages/release/bioc/html/diffloop.html). Can immediately be imported for interactive visualization in [DNAlandscapeR](https://dnalandscaper.aryeelab.org).

So, outputs 4, 5, and 6 are identical except in presentation. These data are a subset of those presented in 3. Interchromosomal interactions from 2 are often discarded by other preprocessing pipelines, but they may hold value. 
If the `qcReport` is generated, then the `.stat` file won't tell you anything new. However, if `R` is not installed on your machine, this will be a useful file for assessing the quality of your library.  

## Configurations<a name="configuration"></a>
Running
```
hichipper --help
```
shows the parameters that can be used in this software package as reproduced below.

```
Usage: hichipper [OPTIONS] MANIFEST

  A preprocessing and QC pipeline for HiChIP data.

Options:
  --out TEXT           Output directory name  [required]
  --min-dist TEXT      Minimum distance ; default = 5000
  --max-dist TEXT      Peak padding width (applied on both left and right);
                       default = 2000000
  --macs2-string TEXT  String of arguments to pass to MACS2; default = "-p
                       0.01 --nomodel"
  --peak-pad TEXT      Peak padding width (applied on both left and right);
                       default = 1500
  --merge-gap TEXT     Max gap size for merging peaks; default = 1500
  --min-qual TEXT      Minimum quality for read; default = 30
  --read-length TEXT   Length of reads from experiment; default = 75
  --keep-temp-files    Keep temporary files?
  --skip-qc            Skip QC report generation? (Requires R + dependent
                       packages (see README))
  --skip-diffloop      Skip diffloop processing of loops? (Requires R +
                       [diffloop](http://bioconductor.org/packages/release/bioc/html/diffloop.html))
  --version            Show the version and exit.
  --help               Show this message and exit.
```
 
Running
```
hichipper --version
```  
will show the version of this package currently installed. 

```
hichipper, version 0.4.0
```
Check the badge up top to see if a newer version is available or try directly through `pip`:

```
pip install hichipper --upgrade
```

Unless these flags are supplied, the pipeline will attempt to run. Minimally sufficient parameters include
the `--out` flag and a `.yaml` file as shown in the example executions. Below are some explanations of the
additional parameters than can be configured when executing the pipeline. 

## Parameter explanations<a name="pe"></a>

Most of the parameter options are fairly straight forward. Running `hichipper --version` or `hichipper --help`
doesn't run the tool but supplies the information noted above. Otherwise, the default run mode requires 
a `.yaml` file supplied in addition to the `--out` parameter, which specifies the output directory of the run. 
Users can decide to customize final output by using boolean flags  `--keep-temp-files`, `--skip-qc`, and  `--skip-diffloop`. 
If the `R` dependencies noted above could not be successfully installed, these last two flags may be useful. Only 
reads that pass the `--min-qual` flag (default = `30`) will be retained in downstream analyses. 

The figure below shows a cartoon of loops, PETs, and anchors and how the various parameters pertain to padding anchors. 

![param](media/param1.png)

As noted in orange, defined peaks are automatically padded by some integer width from the `--peak-pad` flag. By default, 
this pad extends 1500 base pairs in either direction. Padding the peaks boosts the number of PETs that can be mapped to loops. 
For example, `PET 1` would not be considered in loop since the left end of the read does not overlap with the called anchor. However,
it does overlap with the padded peak, so it is retained. When two peaks are close to one another, they may be merged using the 
`--merge-gap` command. As noted in purple, the padded peaks `B` and `C` are sufficiently close to be merged into a single anchor. 
Note that this can lead to some PETs becoming self-ligation (e.g. `1` and `3`). Note, the `--merge-gap` command is equivalent to running 
[bedtools merge -d](http://bedtools.readthedocs.io/en/latest/content/tools/merge.html) on the padded anchors. 

We compared various parameter settings [for the same sample here](https://cdn.rawgit.com/aryeelab/hichipper/master/qcReports/Parameters/peakPlay.hichipper.html). Each sample was processed with a peak pad and merge gap of 250, 500, 1000, and 1500. By default, we've set these parameters at 1500 to mirror those established in ChIA-PET preprocessing. However, the strong retention of reads near the called anchor loci suggest that using smaller parameter values (i.e. 250 or 500 bp) may be optimal for HiChIP analyses to maximize resolution of loop contact loci. 

The `dist` or distance between two peaks is noted in black as the center of two peaks. The `--min-dist` flag is the smallest
and `--max-dist` is the largest integer number that ensures this distance falls between to be considered in a loop. These defaults
are 5Kb and 2Mb as smaller reads are likely self-ligations whereas larger reads are unlikely to be biologically real loops.

Finally, the `--macs2-string` allows users to directly configure the `macs2` call when defining anchors. By default, we use
a lenient threshold for anchors and use a model-free call. 

## Quality control reports
In the [qcReports folder](qcReports), we collect the `.html` QC report files associated with text annotations
from the experiments performed in the [original HiChIP manuscript](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3999.html)
as well as other reports generated by anonymous collaborators that demonstrate libraries that did not prepare well, likely due
to poor in situ ligation. To determine the quality of a new HiChIP library, we recommend comparing the vital statistics and 
interactive tables and figures between existing libraries. 


## Finding differences<a name="loops"></a>
Have you generated a bunch of HiChIP samples and want to see what's different between them? Check out
the [diffloop vignette](https://rpubs.com/caleblareau/diffloop_vignette) for an example analysis
comparing loops from ChIA-PET (the cranky uncle of HiChIP) between K562 and MCF-7. Installation
instructions for this package are shown in the [dependencies](#dependencies) section. 

## Interactive visualization of loops<a name="viz"></a>
One you've (hopefully) assessed that your samples look good, now go visualize them! One option
is to link the `.bedpe` file to the [WashU Genome Browser](http://epigenomegateway.wustl.edu/). Another
option is to upload the `.rds` to our genome topology browser, [DNAlandscapeR](http://dnalandscaper.aryeelab.org). Navigate
to the **Guide** tab to get a sense of how the browser works and ultimately add your sample(s) to a local user session
using the **Import** tab. Note: the browser currently supports hg19/hg37 and mm9 genome builds. 

## Questions/comments/feedback
are always welcomed. Email [Caleb](mailto:caleblareau@g.harvard.edu) anytime! 
