# Noteworthy changes in hichipper

## Split .fastq files as input to HiC-Pro

As of version `0.5.3` of **hichipper**, users should be able to input split `.fastq` files into HiC-Pro
and have **hichipper** function properly. No extra user flags are needed for this functionality. Thanks
to our early users for helping us figure this out. 

## A note on duplicate PETs in loops

In certain versions (`0.4.4` to `0.5.3`) of **hichipper**, duplicates were not being filtered out 
by default. These duplicates were potentially inflating the number of PETs mapping to loops only. 

In version `0.7.1`, peaks are now called using all reads without again filtering duplicates
as this is done during the HiC-Pro step. 

<br><br>