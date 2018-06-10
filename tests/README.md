# hichipper tests

equivalent to running
```
hichipper --out output1 --peak-pad 1000 --skip-resfrag-pad --basic-qc --skip-diffloop yaml/one.yaml

# Now with absolute file paths
hichipper --out output1 --peak-pad 1000 --skip-resfrag-pad --basic-qc --skip-diffloop yaml/one_abs.yaml

```
and testing the output against something we think is correct

## Verify the `call` mode

```
hichipper call --out outputCallBAM --restriction-frags ../RestrictionFragmentFiles/hg19_MboI_resfrag.bed.gz --peaks chipseq/GM12878_SMC3_ChIPSeq_chr22.narrowPeak --skip-resfrag-pad --basic-qc --skip-diffloop --input-bam call_inputs/HiCUP_chr22.bam
hichipper call --out outputCallvi --restriction-frags ../RestrictionFragmentFiles/hg19_MboI_resfrag.bed.gz --peaks chipseq/GM12878_SMC3_ChIPSeq_chr22.narrowPeak --skip-resfrag-pad --basic-qc --skip-diffloop --input-bam call_inputs/dSRR3467177_allValidPairs 

```


## Other useful commands -- peak calling
```
hichipper --out combinedall --skip-diffloop yaml/example_COMBINED_ALL.yaml
hichipper --out combinedself --skip-diffloop yaml/example_COMBINED_SELF.yaml

hichipper --out eachall --skip-diffloop yaml/example_EACH_ALL.yaml
hichipper --out eachself --skip-diffloop yaml/example_EACH_SELF.yaml
```

