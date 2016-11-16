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



def call_to_R():

    cmd = [
        library(doParallel)
        cl < - makeCluster(2)
        registerDoParallel(cl)

        library(wavClusteR)
        # filename <- system.file( "extdata", "example.bam", package = "wavClusteR" )
        filename < - system.file("/scratch/Marta/out.sam", package="wavClusteR")
        Bam < - readSortedBam(filename="/scratch/Marta/out.bam")
        Bam

        ## ----eval=TRUE-----------------------------------------------------------
        countTable < - getAllSub(Bam, minCov=10)
        head(countTable)

        ## ----fig.width=5, fig.height=5, fig.align='center', eval=TRUE------------
        plotSubstitutions(countTable, highlight="TC")

        ## ----eval=FALSE----------------------------------------------------------
        model <- fitMixtureModel(countTable, substitution = "TC")

        ## ----eval=TRUE-----------------------------------------------------------
        data(model)
        str(model)

        ## ----fig.width=7, fig.height=4.5, fig.align='center', eval=TRUE----------
        (support < - getExpInterval(model, bayes=TRUE))

        ## ----fig.width=7, fig.height=4.5, fig.align='center', eval=TRUE----------
        (support < - getExpInterval(model, bayes=FALSE, leftProb=0.9, rightProb=0.9))

        ## ----fig.width=7, fig.height=5, fig.align='center', eval=TRUE------------
        plotSubstitutions(countTable, highlight="TC", model)

        ## ----eval=TRUE-----------------------------------------------------------
        highConfSub < - getHighConfSub(countTable,
                                       support=support,
                                       substitution="TC")
        head(highConfSub)

        ## ----eval=TRUE-----------------------------------------------------------
        coverage < - coverage(Bam)
        coverage$chrX

    ## ----eval=TRUE-----------------------------------------------------------
    clusters < - getClusters(highConfSub=highConfSub,
                             coverage=coverage,
                             sortedBam=Bam,
                             method="mrn",
                             threshold=1,
                             cores=1)
    clusters

    ## ----eval=TRUE-----------------------------------------------------------
    clusters < - getClusters(highConfSub=highConfSub,
                             coverage=coverage,
                             sortedBam=Bam,
                             method="mrn",
                             cores=1)
    clusters

    ## ----eval=TRUE-----------------------------------------------------------
    require(BSgenome.Hsapiens.UCSC.hg19)

    wavclusters < - filterClusters(clusters=clusters,
                                   highConfSub=highConfSub,
                                   coverage=coverage,
                                   model=model,
                                   genome=Hsapiens,
                                   refBase="T",
                                   minWidth=12)

    wavclusters

    ## ----fig.width=5, fig.height=5, fig.align='center', eval=TRUE------------
    plotSizeDistribution(clusters=wavclusters, showCov=TRUE, col="skyblue2")


    ## ----eval=TRUE-----------------------------------------------------------
    sessionInfo()

    ]
