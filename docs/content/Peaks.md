# Peaks with hichipper

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

<img src="https://github.com/aryeelab/hichipper/blob/master/docs/content/media/peakParamPng.png?raw=True" style="width:100%;"><br>

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
`macs2` should be fine. 

# Multiple ChIP-Seq peaks as input

As raised in this [issue](https://github.com/aryeelab/hichipper/issues/18), if you have multiple samples and multiple ChIP-Seq or related high-quality peak definitions to be used as an input, the way to do this is to create two or more `.yaml` files, each one specifying its own bed file of peaks. Then, execute `hichipper` such that you restrict the analysis to the sample you want per bed file using the `--keep-samples` or `--ignore-samples` flags. Thanks to user **sb5169** for bringing this up. 

# HiChIP-Specific Bias Correction

A key difference of HiChIP data compared to ChIA-PET, ChIP-Seq, and related
immunoprecipitation assays is the a notable bias where a greater read density
accumulates near the motif used in the restriction enzyme digestion. The image below
shows the ratio of the treatment to the background (the statistic used in macs2 to call peaks)
as a function of distance to the nearest restriction fragment locus. Note the plot below-- 

<img src="https://github.com/aryeelab/hichipper/blob/master/docs/content/media/hichip_bias.png?raw=True" style="width:100%;"><br>

A more detailed description of this bias and our analysis is contained in this [writeup](https://github.com/aryeelab/hichipper/blob/master/docs/content/media/hichipper_supplement.pdf).

<br><br>
