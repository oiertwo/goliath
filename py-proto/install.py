#!/usr/bin/python3

HG19_LOCATION = '~/hg19/'

from subprocess import call
comands = [
    'sudo apt-get install r-bioc-biomart -y',
    'sudo apt-get install r-bioc-bsgenome -y',
    'sudo apt-get install python3-pip -y',
    'sudo apt-get install bowtie -y',
]

for c in comands:
    call(c.split(' '))

#PIP INSTALL (IN THE RUNNING ENV)
import pip
py_libs = [
    'rpy2',
    'cutadapt',
]

for p in py_libs:
    if pip.main(['install', p]):
        print("{} SUCCESFULY INSTALLED".format(p))

r_commands = [
    'source("https://bioconductor.org/biocLite.R")',
    'library(BiocInstaller)',
    'biocLite("wavClusteR")',
    'install.packages("foreach")'
    'install.packages("doParallel")',
    'biocLite("BSgenome.Hsapiens.UCSC.hg19")',
    'biocLite("XML")',
]

import rpy2.robjects as robjects

for r in r_commands:
    robjects.r(r)


#hg19 for bowtie

download = False

if download:
    cmd = [
        'cd {}'.format(HG19_LOCATION),
        'wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19.ebwt.zip',
        'unzip hg19.ebwt.zip',
        'sh make_hg19.sh',
    ]

else:
    cmd = [
        'mv ../hg19.ebwt.zip {}/hg19.ebwt.zip'.format(),
        'cd {}'.format(HG19_LOCATION),
        'unzip hg19.ebwt.zip',
        'sh make_hg19.sh',
        ]

for c in cmd:
    call(c.split(' '))
