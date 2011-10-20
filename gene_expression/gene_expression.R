library(DESeq)

SAM <- read.delim("/home/ele/Data/454/mapping/all_vs_full_imputed.sam", header=FALSE)
names(SAM)[c(1,3)] <- c("read", "mapped")

## partition the sam in an easy part of reads with single matches and
## a part with problematic reads with multiple matches
easy <- !SAM$read%in%SAM$read[duplicated(SAM$read)]
SAM.easy <- SAM[easy,]
SAM.multi <- SAM[!easy,]

## get counts for the easy single matches
SAM.easy <- merge(SAM.easy[,c("read", "mapped")], R[, c("read", "lib")], all.x=TRUE)


## Get additional info from multiple maps but split according to
## easy.map get the single counts for each of the multiple

##so <- by(SAM.multi, SAM.multi$read , function (x) easy.counts[x$mapped])
## so <- so[!unlist(lapply(so, is.null))]

## so.diff <- unlist(lapply(so, function (x) x/sum(x)))
## names(so.diff) <- gsub(".*\\.", "", names(so.diff))

## added <- sapply(1:length(easy.counts), function (x){
##   sum(easy.counts[x], so.diff[names(easy.counts[x])], na.rm=TRUE)})

## names(added) <- names(easy.counts)


countsTable <- as.data.frame.matrix(table(SAM.easy$mapped, SAM.easy$lib))

## Try to exclude contigs before DESeq
## countsTable2 <- countsTable[grep("Contig", row.names(countsTable)),]

conds <- factor(c("EU", "EU", "L2", "M", "TW", "TW"))

cds <- newCountDataSet(countsTable, conds)
cds <- estimateSizeFactors(cds)
cds <- estimateVarianceFunctions(cds)
res <- nbinomTest(cds, "TW", "EU")

## very absurd these diagnostic plots
## scvPlot(cds)

## residualsEcdfPlot( cds, "TW" )
## residualsEcdfPlot( cds, "EU" )


### This shows nothing as depth is not big enough to allow
### detection of significance
## plot( 
##      res$baseMean, 
##      res$log2FoldChange, 
##      log="x", pch=20, cex=.4, 
##      col = ifelse( res$padj < .1, "red", "black" ),
##      main="My data")

vsd <- getVarianceStabilizedData( cds )
dists <- dist( t( vsd ) )
idists <- as.matrix(dists)
colnames(idists) <- gsub(".easy.counts", "", colnames(idists))
rownames(idists) <- gsub(".easy.counts", "", colnames(idists))
heatmap (idists , symm=TRUE, margins = c (7,7))

sig <- res[res$padj < .1 & !is.na(res$padj),"id"]
countsTable[row.names(countsTable)%in%sig,]
## apply(countsTable[row.names(countsTable)%in%annot.sig,], 2 , sum)

