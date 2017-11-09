# hichipper tests

equivalent to running
```
hichipper --out output1 --peak-pad 1000 --skip-resfrag-pad --skip-qc --skip-diffloop yaml/one.yaml
```
and testing the output against something we think is correct


## Other useful commands
```
hichipper --out combinedall --skip-diffloop yaml/example_COMBINED_ALL.yaml
hichipper --out combinedself --skip-diffloop yaml/example_COMBINED_SELF.yaml

hichipper --out eachall --skip-diffloop yaml/example_EACH_ALL.yaml
hichipper --out eachself --skip-diffloop yaml/example_EACH_SELF.yaml
```

