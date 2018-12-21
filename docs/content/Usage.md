# Simple usage example

The example below uses the test dataset bundled with the **hichipper** package source code.
Download the package, change to the test directory, and execute the basic working example with:

```
git clone https://github.com/aryeelab/hichipper.git
cd hichipper/tests
hichipper --out output1 yaml/one.yaml 
```

Here's a more detailed description of what just happened. First, we had to create a sample description
file that specifies how peaks are to be inferred (in this example, they are pre-specified from a ChIP-Seq experiment). 
Next, one must specify the location of a restriction fragment file. Finally, a path to the HiC-Pro output folder
must be designated. These are encoded
through the `peaks`, `resfrags`, and `hicpro_output` variables that will be parsed from the `.yaml` format. 

1. Create a sample description file:
  
  Description files can be created with the `.yaml` format. 

  **Processing `.yaml` format**
   
   Example [yaml](https://en.wikipedia.org/wiki/YAML) format sample description file:
   
   ```
   peaks:
    - chipseq/GM12878_SMC3_ChIPSeq_chr22.narrowPeak
   resfrags:
    - resfrag/hg19_MboI_resfrag_chr22.bed.gz
   hicpro_output:
    - hicpro
   ```
   
  Note: This file is available as `example.yaml` in the `hichipper/tests` directory.
  
  In this example, we call loops from two GM12878 samples using just chromosome 22
  using pre-determined peaks from a ChIP-Seq file.   
  
  
2. Run the pipeline:

	```
	hichipper --out output1 yaml/one.yaml 
	```

Additional details concerning user configuration options are shown below. 

# More typical example
While the example above references files that are part of the **hichipper** distribution,
our experience using this tool in conjunction with HiC-Pro suggests that a file hierarchy
like the following may be more typical. 

```
ls -LR
```
yields (a slightly modified version of)
```
./hicpro/bowtie_results:
bwt2  bwt2_global  bwt2_local

./hicpro/bowtie_results/bwt2:
SRR3467175  SRR3467176  SRR3467177  SRR3467178

./hicpro/bowtie_results/bwt2/SRR3467175:
SRR3467175_1_hg19.bwt2merged.bam  SRR3467175_2_hg19.bwt2merged.bam  SRR3467175_hg19.bwt2pairs.bam
SRR3467175_1_hg19.mapstat         SRR3467175_2_hg19.mapstat         SRR3467175_hg19.bwt2pairs.pairstat*

./hicpro/bowtie_results/bwt2/SRR3467176:
SRR3467176_1_hg19.bwt2merged.bam  SRR3467176_2_hg19.bwt2merged.bam  SRR3467176_hg19.bwt2pairs.bam
SRR3467176_1_hg19.mapstat         SRR3467176_2_hg19.mapstat         SRR3467176_hg19.bwt2pairs.pairstat*

./hicpro/bowtie_results/bwt2/SRR3467177:
SRR3467177_1_hg19.bwt2merged.bam  SRR3467177_2_hg19.bwt2merged.bam  SRR3467177_hg19.bwt2pairs.bam
SRR3467177_1_hg19.mapstat         SRR3467177_2_hg19.mapstat         SRR3467177_hg19.bwt2pairs.pairstat*

./hicpro/bowtie_results/bwt2/SRR3467178:
SRR3467178_1_hg19.bwt2merged.bam  SRR3467178_2_hg19.bwt2merged.bam  SRR3467178_hg19.bwt2pairs.bam
SRR3467178_1_hg19.mapstat         SRR3467178_2_hg19.mapstat         SRR3467178_hg19.bwt2pairs.pairstat*

./hicpro/hic_results:
data

./hicpro/hic_results/data:
SRR3467175  SRR3467176  SRR3467177  SRR3467178

./hicpro/hic_results/data/SRR3467175:
SRR3467175_hg19.bwt2pairs.DEPairs*    SRR3467175_hg19.bwt2pairs.RSstat*   SRR3467175_hg19.bwt2pairs.SinglePairs*
SRR3467175_hg19.bwt2pairs.DumpPairs*  SRR3467175_hg19.bwt2pairs.SCPairs*  SRR3467175_hg19.bwt2pairs.validPairs*
SRR3467175.allValidPairs*

./hicpro/hic_results/data/SRR3467176:
SRR3467176_hg19.bwt2pairs.DEPairs*    SRR3467176_hg19.bwt2pairs.RSstat*   SRR3467176_hg19.bwt2pairs.SinglePairs*
SRR3467176_hg19.bwt2pairs.DumpPairs*  SRR3467176_hg19.bwt2pairs.SCPairs*  SRR3467176_hg19.bwt2pairs.validPairs*
SRR3467176.allValidPairs*

./hicpro/hic_results/data/SRR3467177:
SRR3467177_hg19.bwt2pairs.DEPairs*    SRR3467177_hg19.bwt2pairs.RSstat*   SRR3467177_hg19.bwt2pairs.SinglePairs*
SRR3467177_hg19.bwt2pairs.DumpPairs*  SRR3467177_hg19.bwt2pairs.SCPairs*  SRR3467177_hg19.bwt2pairs.validPairs*
SRR3467177.allValidPairs*

./hicpro/hic_results/data/SRR3467178:
SRR3467178_hg19.bwt2pairs.DEPairs*    SRR3467178_hg19.bwt2pairs.RSstat*   SRR3467178_hg19.bwt2pairs.SinglePairs*
SRR3467178_hg19.bwt2pairs.DumpPairs*  SRR3467178_hg19.bwt2pairs.SCPairs*  SRR3467178_hg19.bwt2pairs.validPairs*
SRR3467178.allValidPairs*

...

```
where files denoted in with an asterisk* are assumed to exist. Typically, an analysis folder may look like so:

```
fastq/
|-- SRR3467175
|  |-- SRR3467175_1.fastq.gz
|  |-- SRR3467175_2.fastq.gz
|-- SRR3467176
|  |-- SRR3467176_1.fastq.gz
|  |-- SRR3467176_2.fastq.gz
|-- SRR3467177
|  |-- SRR3467177_1.fastq.gz
|  |-- SRR3467177_2.fastq.gz
|-- SRR3467178
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
|  |  |  |-- SRR3467175*RSstat
|  |  |  |-- SRR3467175*Pairs # 5 Files
|  |  |  |-- SRR3467175.allValidPairs
|  |  |-- SRR3467176
|  |  |  |-- SRR3467176*RSstat
|  |  |  |-- SRR3467176*Pairs # 5 Files
|  |  |  |-- SRR3467175.allValidPairs
|  |  |-- SRR3467177
|  |  |  |-- SRR3467177*RSstat
|  |  |  |-- SRR3467177*Pairs # 5 Files
|  |  |  |-- SRR3467175.allValidPairs
|  |  |-- SRR3467178
|  |  |  |-- SRR3467178*RSstat
|  |  |  |-- SRR3467178*Pairs # 5 Files
|  |  |  |-- SRR3467175.allValidPairs
GM12878_SMC3_ChIPSeq.narrowPeak
hg19_MboI_resfrag.bed.gz
yaml/
|-- one.yaml
config-hicpro-mboi-ext12.txt
```
where the results in the `hicpro` directory could have been obtained by running: 
```
HiC-Pro -i fastq/ -o hicpro/ -c config-hicpro-mboi-ext12.txt -p
```
and subsequently executing the resulting `HiCPro_step1_hic.sh` and `HiCPro_step2_hic.sh`. 

Thus, the `yaml/one.yaml` file
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
