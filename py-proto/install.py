#!/usr/bin/python3

HG19_LOCATION = '~/hg19/'

from subprocess import call

comands = [
    'sudo apt-get install r-bioc-biomart -y',
    'sudo apt-get install r-bioc-bsgenome -y',
    'pip3 install -r requirements.txt',
    'sudo apt-get install bowtie -y',
]

for c in comands:
    call(c.split(' '))

r_commands = [
    'source("https://bioconductor.org/biocLite.R")',
    'library(BiocInstaller)',
    'biocLite("wavClusteR")',
    'install.packages("doParallel")',
    'biocLite("BSgenome.Hsapiens.UCSC.hg19")',
    'biocLite("XML")',
]

import rpy2.robjects as robjects

for r in r_commands:
    robjects.r(r)


#hg19 for bowtie

cmd = [
    'cd {}'.format(HG19_LOCATION),
    'wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19.ebwt.zip',
    'unzip hg19.ebwt.zip',
    'sh make_hg19.sh',
]