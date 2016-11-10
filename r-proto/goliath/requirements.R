#sudo apt install r-bioc-biomart r-bioc-bsgenome -y
source("https://bioconductor.org/biocLite.R")
library(BiocInstaller)
biocLite("wavClusteR")
install.packages("doParallel")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
biocLite("XML")
