#!/usr/bin/python3
import os
import datetime

### CONFIG

#general settings
RANDOM = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M")
RANDOM = ""
DEST = "/home/ibai/output-{}".format(RANDOM)
MERGE_OUT = os.path.join(DEST,"merge_out.fastq")
BAM_FILE = os.path.join(DEST,"out.bam")
#bowtie settings
BOWTIE_HG19 = '/home/ibai/hg19' #please write full path
BOWTIE_OUT = os.path.join(DEST, 'bow_out.sam')

#cutadapt settings
CA_3ADAPT = "TGGAATTCTCGGGTGCCAAGG"
CA_LENGTH = 15
CA_OUT = os.path.join(DEST,"ca_out.fastq")
CA_PRINT_FILE = os.path.join(DEST, "trim_stats.txt")

#########

### SCRIPT

#!/usr/bin/env python


import argparse
import sys
import script as tools
from subprocess import call

def create_argparser():
    parser = argparse.ArgumentParser(
        description='a FAST CLIP pypeline made easy'
    )

    return parser


def main():
    #parser = create_argparser()

    #try:
    #    args = parser.parse_args()
    #except argparse.ArgumentError as exc:
    #    print('Error parsing arguments.')
    #    parser.error(str(exc.message))
    #    exit(-1)
    if not os.path.exists(DEST):
        os.makedirs(DEST)
    try:

        argc = len(sys.argv[:])
        print("merging... please wait, this is a good time for a coffee or a tea...")

        if argc == 3:
            print("Merging 2 files...")
            tools.merge2(sys.argv[1],sys.argv[2],MERGE_OUT)
        elif argc == 4:
            print("Merging 3 files...")
            tools.merge3(sys.argv[1],sys.argv[2],sys.argv[3],MERGE_OUT)
        elif argc == 2:
            globals()['MERGE_OUT'] =  sys.argv[1]
            print("No need to merge, only one file...")
        else:
            print("HELP:")
            print("arguments missing")
            print("structure: ")
            print("pype.py [file1] [file2] ... [fileN]")
            print("(hint) all file extension should be *.fastq")
            print("(hint 2) Maximum: 3 files")
            return


        print("NEW PROCESS: Cutadapt")

        cmd = "cutadapt -a {} -m {} -o {} {}".format(CA_3ADAPT,CA_LENGTH,CA_OUT,MERGE_OUT)
        print("CALL: " + cmd)
        call(cmd.split(' '))

        print("NEW PROCESS: Aligning... time to go for lunch")


        cmd = 'bowtie {} --threads 4 -v 2 -m 10 --best --strata {} -S {}'.format(os.path.join(BOWTIE_HG19,'hg19'),
                                                                             CA_OUT, BOWTIE_OUT)
        print("CALL: " + cmd)
        a = []
        for i in cmd.split(' '):
            if i:
                a.append(i)
        call(a)

        ## 	#ALIGN:
        ## 		sample.fastq -> sample.sam
        ## 	#CONVERT:
        ## 		samtools view -b -S sample.sam -o sample.bam
        ## 	#SORT:
        ## 		samtools sort sample.bam sample_sorted
        ## 	#INDEXING:
        ## 		samtools index sample_sorted.bam
        ## #-F 4:
        ##     leaves the unmaped out

        cmd = "samtools view -b -S {} -o {}".format(BOWTIE_OUT, BAM_FILE)
        print("CALL: " + cmd)
        call(cmd.split(" "))

        tools.call_to_R(BAM_FILE, output=DEST)

    except:
        print("ERROR:", sys.exc_info())

        print("HELP:")
        print("arguments missing????")
        print("structure: ")
        print("biomerge [file1] [file2] ... [fileN] [mergedfile]")
        print("(hint) all file extension should be *.fastq")
        print("(hint 2) Maximum: 3 files")


if __name__ == '__main__':
    main()
