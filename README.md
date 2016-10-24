# goliath
Biomedical tools 

#links referencia:

[bowtie](http://bowtie-bio.sourceforge.net/index.shtml) : software para alinear (y mas..) <br\>
```bash
zcat ${tmp_fastq1_gz} | bowtie "$tmp_genome_path" --threads 4 -v 2 -m 10 --best --strata - -S ${tmp_sampleID}.1.sam 
```
[wavClusteR](https://cran.r-project.org/web/packages/doParallel/vignettes/gettingstartedParallel.pdf) : la opcion de R 


