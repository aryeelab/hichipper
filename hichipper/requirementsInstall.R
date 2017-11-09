install_pkgs <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
}
install_pkgs(c("data.table", "devtools", "foreach", "ggplot2", "knitr", "networkD3", "readr", "reshape2"))

source("https://bioconductor.org/biocLite.R")
install_pkgs_bioc <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) biocLite(new.pkg, dependencies = TRUE)
}
install_pkgs_bioc(c("diffloop"))