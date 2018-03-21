# Visualizing loops

## Visualization in UCSC / WashU Epigenome Browser<a name="vizUCSC"></a>
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

#### For WashU Epigenome Browser

We recommend building a WashU track hub to facilitate data i/o and batch visualization. 
[This Rscript]() may help in creating the .json file for visualization.

## Old method of visualizing loops<a name="viz"></a>
We've publicly distributed an R/Shiny application for visualizing DNA loops called [DNAlandscapeR](https://molpath.shinyapps.io/DNAlandscapeR).
**While DNAlandscapeR is not actively maintained**, the latest build is stable and should facilitate visualizing `.rds`
files. The [code is all made publicly available here](https://github.com/aryeelab/dnalandscaper).

To visualize loops in this browser, navigate to the **Guide** tab to get a sense of how the browser works and ultimately add your sample(s) to a local user session
using the **Import** tab. Note: the browser currently supports hg19/hg37 and mm9 genome builds. 


<br><br>