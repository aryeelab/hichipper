# Visualizing loops

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

<br><br>