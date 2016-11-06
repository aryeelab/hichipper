#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

scriptdir <- args[1]
outdir <- args[2]
cwd <- args[3]
version <- args[4]
samples <- unlist(strsplit(args[seq(-1,-4)], split = " "))

rmarkdown::render(paste0(scriptdir,"/qcReport_make.Rmd"), params = list(
  scriptdir = scriptdir, 
  outdir = outdir, 
  cwd = cwd, 
  samples = samples,
  version = version
), output_file=paste0(cwd, "/", outdir, "/", "hichipper-qcReport.html"))