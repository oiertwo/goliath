from Bio import SeqIO
import os
import sys
import itertools
import gzip
import shutil
from itertools import zip_longest

def merge2(f1, f2, f_out, compress=True):

    outfile = open(f_out, "w")
    iter1 = SeqIO.parse(open(f1), "fastq")
    iter2 = SeqIO.parse(open(f2), "fastq")

    for t in zip_longest(iter1, iter2):
        add = []
        for rec in t:
            if (rec != None):
                add.append(rec)
        SeqIO.write(add, outfile, "fastq")
    outfile.close()

    # compress
    if compress:
        with open(f_out, 'rb') as f_in, gzip.open("{}.gz".format(f_out), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def merge3(f1, f2, f3, f_out, compress=True):

    outfile = open(f_out, "w")
    iter1 = SeqIO.parse(open(f1), "fastq")
    iter2 = SeqIO.parse(open(f2), "fastq")
    iter3 = SeqIO.parse(open(f3), "fastq")


    for t in zip_longest(iter1, iter2, iter3):
        add = []
        for rec in t:
            if (rec != None):
                add.append(rec)
        SeqIO.write(add, outfile, "fastq")
    outfile.close()

    # compress
    if compress:
        with open(f_out, 'rb') as f_in, gzip.open("{}.gz".format(f_out), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)



def call_to_R(bam, output="./output"):
    cmd = [
        'library(doParallel)',
        'cl < - makeCluster(2)',
        'registerDoParallel(cl)',
        'library(wavClusteR)',
        #'filename < - system.file("{}", package="wavClusteR")'.format(bam),
        'Bam < - readSortedBam(filename="{}")'.format(bam),
        'setwd("{}")'.format(output),
        'countTable < - getAllSub(Bam, minCov=10)',
        "png('substitutions_1.png')",
        'plotSubstitutions(countTable, highlight="TC")',
        'dev.off()',
        'model < - fitMixtureModel(countTable, substitution="TC")',
        'data(model)',
        'str(model)',
        '(support < - getExpInterval(model, bayes=TRUE))',
        "png('model_densities.png')",
        '(support < - getExpInterval(model, bayes=FALSE, leftProb=0.9, rightProb= 0.9 ) )',
        'dev.off()',
        "png('substitutions_2.png')",
        'plotSubstitutions( countTable, highlight = "TC", model )',
        'dev.off()',
        'highConfSub <- getHighConfSub(countTable,support=support,substitution="TC")',
        'coverage < - coverage(Bam)',
        'clusters < - getClusters(highConfSub=highConfSub,coverage=coverage,sortedBam=Bam,method="mrn",cores=1)',
        'write.csv(as.data.frame(clusters), file = "clusters.csv")',
        'require(BSgenome.Hsapiens.UCSC.hg19)',
        'wavclusters < - filterClusters(clusters=clusters,highConfSub=highConfSub,coverage=coverage,model=model,genome=Hsapiens,refBase="T",minWidth=12)',
        'write.csv(as.data.frame(wavclusters), file = "wavclusters.csv")',
        "png('distribution.png')",
        'plotSizeDistribution(clusters=wavclusters, showCov=TRUE, col="skyblue2")',
        'dev.off()',
        "png('statistics.png')",
        'plotStatistics(clusters=wavclusters,corMethod="spearman",lower=panel.smooth)',
        'dev.off()',
    ]

    import rpy2.robjects as robjects

    for r in cmd:
        robjects.r(r)