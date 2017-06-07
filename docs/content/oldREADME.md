<p align="left">
  <br><br><br>
  <img src="media/hichipper_logo.png" width="50%"/>
</p>

This package is maintained by [Caleb Lareau](mailto:caleblareau@g.harvard.edu) in the [Aryee Lab](http://aryee.mgh.harvard.edu/). Source code is made freely available here and a packaged install version is provided through [PyPi](https://pypi.python.org/pypi/hichipper/).

[![Build Status](https://travis-ci.org/aryeelab/hichipper.svg?branch=master)](https://travis-ci.org/aryeelab/hichipper)
[![PyPI version](https://img.shields.io/badge/pypi-0.6.0-brightgreen.svg)](https://pypi.python.org/pypi/hichipper/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## About<a name="about"></a>

The **hichipper** package implements our data processing and quality control pipeline for 
[HiChIP](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3999.html) data.
This package takes output from a [HiC-Pro](https://github.com/nservant/HiC-Pro) run and a sample manifest file (`.yaml`) that coordinates optional high-quality peaks (identified through ChIP-Seq) and restriction fragment locations (see [folder here](RestrictionFragmentFiles))
as input and produces output that can be used to 1) determine library quality, 2) identify and characterize DNA loops and 3) interactively visualize loops. Loops are assigned strength and confidence metrics that can be used to evaluate samples individually or for differential analysis in
downstream tools like [diffloop](http://github.com/aryeelab/diffloop).
We have used the library QC metrics with as few as 1 million reads, enabling library quality to be assessed through shallow (and cheap) sequencing before performing a full depth sequencing run.

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
- [Peaks](#HPC)
- [HiChIP-Specific Bias Correction](#hcsbc)
- [Configurations](#configuration)
- [Restriction-fragment aware padding](#rfap)
- [Parameter Explanations](#pe)
- [User parameter recommendations](#ur)
- [Quality Control reports](#qcr)
- [Interactive visualization of loops](#viz)
- [Visualization in UCSC](#vizUCSC)
- [Analyzing loops in the R](#loops)

## Workflow Overview<a name="ugo"></a>
A simple graphical guide to processing HiChIP data is shown below. The role of **hichipper**
is to import aligned read files from (e.g. [HiC-Pro](https://github.com/nservant/HiC-Pro))
as well as location of restriction fragment files
([available here](https://github.com/aryeelab/hichipper/tree/master/RestrictionFragmentFiles)) coordinated through a
`.yaml` configuration file and produce user-friendly output. 

In particular, **hichipper** allows users to pre-supply their own set of gold-standard peaks (e.g. from ChIP-Seq)
or call peaks directly from HiChIP data using a novel background detection algorithm. In either case, interactions
and chromatin loops can be called using a restriction fragment-aware approach that substantially increases read density in loops. 
 
![hichipper_overview](media/Overview.png)

A higher resolution [slide of this image](media/Overview.pptx) is in the [media](media) folder.




## Output

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

6. `*.rds` The same set of loops as 4 but in an R binary compressed format of a `loops()` S4 object from [diffloop](http://bioconductor.org/packages/release/bioc/html/diffloop.html). Can
immediately be imported for interactive visualization in [DNAlandscapeR](https://dnalandscaper.aryeelab.org).

So, outputs 4, 5, and 6 are identical except in presentation. These data are a subset of those presented in 3. Interchromosomal interactions from 2 are often discarded by other preprocessing pipelines, but they may hold value. 
If the `qcReport` is generated, then the `.stat` file won't tell you anything new. However, if `R` is not installed on your machine, this will be a useful file for assessing the quality of your library.  

## Peaks <a name="HPC"></a>

To call peaks from HiChIP data directly, **hichipper** aggregates read density from either all samples or each sample individually. Additionally,
users can specify whether all read density is used or if only self-ligation reads are used. To specify these options, put the appropriate 
string of the form `{COMBINED,EACH},{ALL,SELF}` in the `peaks` slot of the `.yaml`.

For example, to replicate the peak calling performed in Mumbach *et al.*, one would use the following `.yaml`: 

```
peaks:
  - COMBINED,SELF
resfrags:
  - hg19_MboI_resfrag.bed.gz
hicpro_output:
  - hicpro
```

Alternatively, we can call peaks from the HiChIP data for each sample individually using all reads using this specification--

```
peaks:
  - EACH,ALL
resfrags:
  - hg19_MboI_resfrag.bed.gz
hicpro_output:
  - hicpro
```

The figure below shows all options for peak specification in **hichipper** including every option
for inferring peaks which are noted in the table.

![peakParam](media/peakParamPng.png)

Alternatively, users can pre-specify a set of peaks to used. In this case, a "connectome" will be inferred between the peaks
specified in the `.bed` file. Of note, pre-specified peaks will still be padded either by fixed amounts or to the edges of the restriction 
fragment pads (or both) unless the user specifies these flags differently (see below). 

```
peaks:
  - predeterminedPeaks.bed
resfrags:
  - hg19_MboI_resfrag.bed.gz
hicpro_output:
  - hicpro
```

Note: the input of pre-determined peaks does not have to explicitly be a `.bed` file. Rather, any file name is acceptable so long
as the first three columns indicate appropriate genomic loci as if it were a `.bed` file. For example, `.narrowPeak` files from 
macs2 should be fine. 

## HiChIP-Specific Bias Correction<a name="hcsbc"></a>

A key difference of HiChIP data compared to ChIA-PET, ChIP-Seq, and related
immunoprecipitation assays is the a notable bias where a greater read density
accumulates near the motif used in the restriction enzyme digestion. The image below
shows the ratio of the treatment to the background (the statistic used in macs2 to call peaks)
as a function of distance to the nearest restriction fragment locus. Note the plot below-- 

![hichipPlot](media/hichip_bias.png)

A more detailed description of this bias and our analysis is contained in this [writeup.](media/hichipper_supplement.pdf)

## Quality control reports
In the [qcReports folder](qcReports), we collect the `.html` QC report files associated with text annotations
from the experiments performed in the [original HiChIP manuscript](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3999.html)
as well as other reports generated by anonymous collaborators that demonstrate libraries that did not prepare well, likely due
to poor in situ ligation. To determine the quality of a new HiChIP library, we recommend comparing the vital statistics and 
interactive tables and figures between existing libraries. 

## Finding differences<a name="loops"></a>
Have you generated a bunch of HiChIP samples and want to see what's different between them? Check out
the [diffloop vignette](https://rpubs.com/caleblareau/diffloop_vignette) for an example analysis
comparing loops from ChIA-PET (a similar 3C method to HiChIP) between K562 and MCF-7. Installation
instructions for this package are shown in the [dependencies](#dependencies) section. 

## Interactive visualization of loops<a name="viz"></a>
One you've (hopefully) assessed that your samples look good, now go visualize them! One option
is to link the `.bedpe` file to the [WashU Genome Browser](http://epigenomegateway.wustl.edu/). Another
option is to upload the `.rds` to our genome topology browser, [DNAlandscapeR](http://dnalandscaper.aryeelab.org). Navigate
to the **Guide** tab to get a sense of how the browser works and ultimately add your sample(s) to a local user session
using the **Import** tab. Note: the browser currently supports hg19/hg37 and mm9 genome builds. 

## Visualization in UCSC<a name="vizUCSC"></a>
Users can specify the `--make-ucsc` flag to produce output that can be imported into UCSC.
See [this discussion](https://groups.google.com/a/soe.ucsc.edu/forum/#!topic/genome/kE2pIZUvfnA) for an overview of the format.
In order to produce this output, **hichipper** needs access to the [htslib](http://www.htslib.org/download/) suite of tools
in the computational environment. You can see if you have these dependences available (namely, `tabix` and `bgzip`) by making
sure the following works:

```
tabix --version
```

Specifying this flag will create the additional files `*.txt.gz` and `.txt.gz.tbi`, which can be used to make a UCSC track. 
(Shout out to Gary for helping us with this!)

