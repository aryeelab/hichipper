# Inputs for hichipper

Based on the original design of `hichipper` users needed to specify the full HiC-Pro output. 
However, in subsequent versions, different inputs are accepted; specifically, only a valid fragments file. 

While `.bam` files processed via different software (e.g. HiCUP) are not explicitly supported, one 
can easily convert `.bam` file formats into valid fragment pairs files using a combination of samtools and bedtools
commands, such as this shown below:

```
samtools sort -n HiCUP_chr22.bam  | bamToBed -i - -bedpe |
	awk 'OFS="\t" {print $7,$1,int(($2+$3)/2), $9, $4, int(($5+$6)/2), $10}' > sample_converted.allValidPairs
```

The resulting `sample_converted.allValidPairs` file can be used as input to `hichipper`, using the `--input-vi` flag.

<br><br>