# Noteworthy changes in hichipper

## Split .fastq files as input to HiC-Pro

As of version `0.5.3` of **hichipper**, users should be able to input split `.fastq` files into HiC-Pro
and have **hichipper** function properly. No extra user flags are needed for this functionality. Thanks
to our early users for helping us figure this out. 

## A note on duplicates

In certain versions (`0.4.4` to `0.5.3`) of **hichipper**, duplicates were not being filtered out 
by default. 
