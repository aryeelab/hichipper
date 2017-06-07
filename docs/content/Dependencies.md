# Required software for hichipper

The following dependencies need to be installed before running **hichipper**:
[bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html), OpenSSL, libcurl, libxml2,
and [samtools](http://www.htslib.org/download/). Depending on if you want some bonus functionality,
you may need to download additional requirements. 

Except for `bedtools` and `samtools`, these other dependencies came out of the
box with the unix/linux systems that we've used **hichipper** on. 

But just to be safe, on an Ubuntu system, all of the dependencies can be installed with:
```
apt-get install bedtools libssl-dev libcurl4-openssl-dev libxml2-dev
```

Additionally, `R` must be available in the environment as well as a reasonably recent
version of [pandoc](http://pandoc.org) and a few packages that can be downloaded running the following
in an 'R' environment. :

```
install_pkgs <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
}
install_pkgs(c("DT", "data.table", "devtools", "foreach", "ggplot2", "knitr", "networkD3", "readr", "reshape2"))

source("https://bioconductor.org/biocLite.R")
install_pkgs_bioc <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) biocLite(new.pkg, dependencies = TRUE)
}
install_pkgs_bioc(c("diffloop"))
```

or simply download [this R script](https://github.com/aryeelab/hichipper/blob/master/hichipper/requirementsInstall.R) and run:
```
Rscript requirementsInstall.R
```

Convenient [pandoc binaries](https://s3.amazonaws.com/rstudio-buildtools/pandoc-1.12.4.2.zip) for Linux,
Mac and Windows are available for download from RStudio. If you prefer to install pandoc globally
on your machine, installation instructions can be found [here](http://pandoc.org/installing.html).

<br><br>