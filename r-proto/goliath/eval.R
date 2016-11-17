## ----eval=FALSE----------------------------------------------------------
## 	#ALIGN:
## 		sample.fastq -> sample.sam
## 	#CONVERT:
## 		samtools view -b -S sample.sam -o sample.bam
## 	#SORT:
## 		samtools sort sample.bam sample_sorted
## 	#INDEXING:
## 		samtools index sample_sorted.bam

## ----eval=TRUE-----------------------------------------------------------

library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)

library(wavClusteR)
#filename <- system.file( "extdata", "example.bam", package = "wavClusteR" )
filename <- system.file( "/scratch/Marta/out.sam", package = "wavClusteR" )
Bam <- readSortedBam(filename = "/scratch/Marta/out.bam")
Bam

## ----eval=TRUE-----------------------------------------------------------
countTable <- getAllSub( Bam, minCov = 10 )
head( countTable )

## ----fig.width=5, fig.height=5, fig.align='center', eval=TRUE------------
png('substitutions_1.png')
plotSubstitutions( countTable, highlight = "TC" )
dev.off()

## ----eval=FALSE----------------------------------------------------------
model <- fitMixtureModel(countTable, substitution = "TC")

## ----eval=TRUE-----------------------------------------------------------
data(model)
str(model)

## ----fig.width=7, fig.height=4.5, fig.align='center', eval=TRUE----------
(support <- getExpInterval( model, bayes = TRUE ) )

## ----fig.width=7, fig.height=4.5, fig.align='center', eval=TRUE----------
png('model_densities.png')
(support <- getExpInterval( model, bayes = FALSE, leftProb = 0.9, rightProb = 	0.9 ) )
dev.off()

## ----fig.width=7, fig.height=5, fig.align='center', eval=TRUE------------
png('substitutions_2.png')
plotSubstitutions( countTable, highlight = "TC", model )
dev.off()

## ----eval=TRUE-----------------------------------------------------------
highConfSub <- getHighConfSub( countTable, 
                               support = support, 
                               substitution = "TC" )
head( highConfSub )                               

## ----eval=FALSE----------------------------------------------------------
## highConfSub <- getHighConfSub( countTable,
##                                supportStart = 0.2,
##                                supportEnd = 0.7,
##                                substitution = "TC" )
## head( highConfSub )

## ----eval=TRUE-----------------------------------------------------------
coverage <- coverage( Bam )
coverage$chrX

## ----eval=TRUE-----------------------------------------------------------
clusters <- getClusters( highConfSub = highConfSub,
                         coverage = coverage,
                         sortedBam = Bam,
                         method = "mrn",
                         cores = 1 )
write.csv(as.data.frame(clusters), file = "clusters.csv")

## ----eval=TRUE-----------------------------------------------------------
require(BSgenome.Hsapiens.UCSC.hg19)

wavclusters <- filterClusters( clusters = clusters, 
                               highConfSub = highConfSub,
                               coverage = coverage, 
                               model = model, 
                               genome = Hsapiens, 
                               refBase = "T", 
                               minWidth = 12)

write.csv(as.data.frame(wavclusters), file = "wavclusters.csv")

## ----eval=FALSE----------------------------------------------------------
## exportHighConfSub( highConfSub = highConfSub,
##                    filename = "hcTC.bed",
##                    trackname = "hcTC",
##                    description = "hcTC" )

## ----eval=FALSE----------------------------------------------------------
## exportClusters( clusters = wavclusters,
##                 filename = "wavClusters.bed",
##                 trackname = "wavClusters",
##                 description = "wavClusters" )

## ----eval=FALSE----------------------------------------------------------
## exportCoverage( coverage = coverage, filename = "coverage.bigWig" )

## ----eval=FALSE----------------------------------------------------------
## txDB <- makeTxDbFromUCSC(genome = "hg19", tablename = "ensGene")

## ----eval=FALSE----------------------------------------------------------
## annotateClusters( clusters = wavclusters,
##               txDB = txDB,
##               plot = TRUE,
##               verbose = TRUE)

## ----eval=FALSE----------------------------------------------------------
## getMetaGene( clusters = wavclusters,
##              txDB = txDB,
##              upstream = 1e3,
##              downstream = 1e3,
##              nBins = 40,
##              nBinsUD = 10,
##              minLength = 1,
##              plot = TRUE,
##              verbose = TRUE )

## ----eval=FALSE----------------------------------------------------------
## getMetaTSS( sortedBam = Bam,
##             txDB = txDB,
##             upstream = 1e3,
##             downstream = 1e3,
##             nBins = 40,
##             unique = FALSE,
##             plot = TRUE,
##             verbose = TRUE )

## ----fig.width=5, fig.height=5, fig.align='center', eval=TRUE------------
png('distribution.png')
plotSizeDistribution( clusters = wavclusters, showCov = TRUE, col = "skyblue2" )
dev.off()

## ----fig.width=5, fig.height=5, fig.align='center', eval=FALSE-----------
png('statistics.png')
plotStatistics( clusters = wavclusters,
                 corMethod = "spearman",
                lower = panel.smooth )
dev.off()

## ----eval=TRUE-----------------------------------------------------------
sessionInfo()

fclusters <- filterClusters( clusters = clusters,
                             highConfSub = highConfSub,
                             coverage = coverage,
                             model = model,
                             genome = Hsapiens,
                             refBase = 'T',
                             minWidth = 12 )
annotateClusters(fclusters, txDB = NULL, genome = "hg19", tablename = "ensGene", plot = TRUE, verbose = TRUE)