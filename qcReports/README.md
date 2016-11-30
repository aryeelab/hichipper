# QC Reports
In this folder hierarchy, we hub all QC Reports from HiChIP that we can get our hands on to facilitate comparisons. We've pre-categorized current reports into libraries that look good and bad. Additionally, we've done some analyses to determine some noteworthy features of the HiChIP experiment/sequencing run (see below). Have more reports to contribute? We'd greatly appreciate it if you submitted a pull request or emailed [Caleb](mailto:caleblareau@g.harvard.edu) with your report!  

When planning a HiChIP experiment and sequencing run, we consider three key factors: 1) cellular input, 2) read length, and 3) read count.  

## Variable Cell Input
[Compiled HTML here](https://cdn.rawgit.com/aryeelab/hichipper/master/qcReports/Noteworthy/cellInput.hichipper.html)

While each sample looks mostly successful, we note the lower proportion of long range interactions in the 5 million cell count relative 
to the other samples. As the 1 million cell sample looks slightly better, we cannot determine, at least from these samples,
that reduced cellular input has a direct effect on library quality.  

## Variable Read Length (75bp versus 100bp)
[Compiled HTML here](https://cdn.rawgit.com/aryeelab/hichipper/master/qcReports/Noteworthy/readLength.hichipper.html)

The [original HiChIP paper](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3999.html) only sequenced using 75 bp PE reads for all samples. To infer the impact of read length, we sequencing a high quality library at 100 bp PE reads and then trimmed these reads to 75 bp. 
From our analyses, this causes a 5-6% reduction in mappable interactions and ~12-13% reduction in reads in loops. Depending on the availability of sequencing technology, we suggest using 75 bp reads with more depth may be a more economical experimental strategy. 

## Downsampled (Variable Read Depth)
[Compiled HTML here](https://cdn.rawgit.com/aryeelab/hichipper/master/qcReports/Noteworthy/downsampled.hichipper.html)

Three samples were created from random downsampling of the GM12878 HiChIP experiment with
1 million, 500,000 and 250,000 reads. As the proportion of long range interactions remains
relatively constant, we suggest that a relatively small number of reads may be sufficient 
to determine library quality before performing a deeper sequencing run. This comparison is just
a reorganization of data in the Good/Mumbach_etal/ folder. The original SRR files shown here include
SRR3467183 (10 million), SRR3467185 (5 million), and SRR3467187 (1 million). 