# Output

## Per-run output files
Each time the user runs **hichipper**, a `*.hichipper.log` file containing information pertaining to the 
flow of the software execution is placed in the `out` directory (specified by the `--out` flag). Unless
otherwise specified, a file ending in `hichipper-qcReport.html` provides an interactive quality control report for all samples. 

## Per-sample output files
Per sample, six (yes, 6, but don't worry-- there's lots of redundancy) output files are created. They are:

1. `*.stat` Key summary statistics that show the number of PETs meeting certain criteria

2. `*.inter.loop_counts.bedpe` Interchromosomal looping between anchor loci. 

3. `*.intra.loop_counts.bedpe` Intrachromosomal looping between **all**  anchor loci

4. `*.filt.intra.loop_counts.bedpe` Intrachromosomal looping between anchor loci where loops meet min/max distance requirements.

5. `*interactions.all.mango` The same set of loops as 4 but with per-loop FDR measures from the loop proximity bias correction algorithm originally implemented in [Mango](https://github.com/dphansti/mango) and presented in the same format. 

6. `*.interaction.txt.gz` The same set of loops as 4 but in a compressed format for visualizing with WashU and UCSC Genome Browser tools.

So, outputs 4, 5, and 6 are identical except in presentation. These data are a subset of those presented in 3. Interchromosomal interactions from 2 are often discarded by other preprocessing pipelines, but they may hold value. 
If the `qcReport` is generated, then the `.stat` file won't tell you anything new. However, if `R` is not installed on your machine, this will be a useful file for assessing the quality of your library.  

<br><br>