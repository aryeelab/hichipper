# hichipper executable scripts

Below are some details of what's what in **hichipper**

### cli.py
The master coordinator file for running **hichipper** as implemented. See details 
[here](https://pythonhosted.org/pyCLI/). This is the script that is called when 
**hichipper** is evoked from the command line.

If interested in how peak files are coordinated, the `.yaml` is parsed, and how all the
rest of the scripts are stitched together, this is the place to look. 

### diffloop_work.R
R script that builds the `.rds` file and the `.mango` file. All filtered interaction 
files are read and and individually processed per sample.

### interactionsCall.sh
This is the main worker script given some set of defined anchors. Most user arguments
related to loops will be passed here. This will also generate the `.stat` file and a 
majority of the log file. 

### lambdaProcess.R
This script will perform the modified background correction implemented in **hichipper** where 
notably a smoothing spline function is fit over the mean treatment / background lambda. 

### qcReport.R
Token wrapper R script that renders the actual working script--

### qcReport_make.Rmd
Using `knitr` and related packages, this generated the `.html` file associated with 
the quality control report. It's pretty messy and peculiar so proceed with caution 
when reading it. 

### requirementsInstall.R
R script that isn't actually called at any point by **hichipper** but is useful for installing
dependencies as coordinated through the main README in github. 

### resFragAnchorProcess.R
This R script will adaptively pad a given bed file to the edges of the restriction fragments,
which should also be supplied by a user defined input. 

