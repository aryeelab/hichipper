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


## Finding differences<a name="loops"></a>
Have you generated a bunch of HiChIP samples and want to see what's different between them? Check out
the [diffloop vignette](https://rpubs.com/caleblareau/diffloop_vignette) for an example analysis
comparing loops from ChIA-PET (a similar 3C method to HiChIP) between K562 and MCF-7. Installation
instructions for this package are shown in the [dependencies](#dependencies) section. 
