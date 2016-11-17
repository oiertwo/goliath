#sudo apt install r-bioc-biomart r-bioc-bsgenome -y
#http://cran.cnr.berkeley.edu/bin/linux/ubuntu/
source("https://bioconductor.org/biocLite.R")
#library(BiocInstaller)
#nbiocLite("BiocUpgrade")
#biocLite("wavClusteR")
install.packages("doParallel")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
biocLite("XML")
install.packages("devtools")
library(devtools)
install.packages("survival")
install_github(repo = 'Bioconductor-mirror/BiocGenerics', ref = "release-3.4")
install_github(repo = 'Bioconductor-mirror/S4Vectors', ref = "release-3.4")
install_github(repo = 'Bioconductor-mirror/IRanges', ref = "release-3.4")
install_github(repo = 'Bioconductor-mirror/XVector', ref = "release-3.4")
install_github(repo = 'Bioconductor-mirror/GenomicRanges', ref = "release-3.4")
install_github(repo = 'Bioconductor-mirror/Biostrings', ref = "release-3.4")
install_github(repo = 'Bioconductor-mirror/SummarizedExperiment', ref = "release-3.4")
install_github(repo = 'Bioconductor-mirror/GenomicAlignments', ref = "release-3.4")
install_github(repo = 'Bioconductor-mirror/wavClusteR', ref = "release-3.4")

install_github(repo = "FedericoComoglio/wavClusteR", ref = "devel", build_vignettes = TRUE)




