## read.pop.var <- function (path) {
##   VAR <- read.delim(path, header=FALSE, as.is=TRUE)
##   names(VAR) <- c("contig", "base", "from", "to", "nfrom", "nto", "perc", "sfrom", "sto", "qfrom", "qto", "pval")
##   ## Not a iupac ambiguity in reference
##   VAR <- subset(VAR, VAR$from%in%c("A", "C", "G", "T"))
##   VAR <- VAR[c(1:6)]
##   return(VAR)
## }

## files <- list.files("/home/ele/Data/454/mapping/mapping_each_lib/",
##                     ".varsnp", full.name=TRUE)

## VARpop.l <- list()
## for (f in files){
##   na <- gsub("final.(\\d+|\\w+).*", "\\1", basename(f))
##   VARpop.l[[na]] <- read.pop.var(f)
##   names(VARpop.l[[na]])[5:6] <- paste(na, names(VARpop.l[[na]])[5:6], sep="-")
## }

## VP <- merge(VARqp[,1:6], VARpop.l[[1]],  by=c("contig", "base", "from", "to"), all.x=TRUE)
## for (i in 2:length(VARpop.l)){
##   VP <- merge(VP, VARpop.l[[i]],  by=c("contig", "base", "from", "to"), all.x=TRUE)
## }

## lib.l <- (apply(VP, 1, function (x) {
##   names(x[7:length(x)])[!is.na(x[7:length(x)])]}))

## names(lib.l) <- apply(VP, 1, function (x) paste (x[1:4], collapse="."))

## ## to see what is shared between each lib
## ## table(unlist(lapply(lib.l, paste, collapse="-")))

## pop.l <- lapply(lib.l, function (x) gsub( "KS4F-nto|UW07F-nto|L2R3-nto", "EURO", x))
## pop.l <- lapply(pop.l, function (x) gsub( "M175-nto|179F-nto|10F-nto", "ASIA", x))

## ## to see what is shared between each lib
## table(unlist(lapply(pop.l, paste, collapse="-")))

## e.l <- lapply(pop.l, function (x) {
##   c(EU=length(grep("EURO", x)),
##     TW=length(grep("ASIA", x)))
## })

## summary.factor(unlist(lapply(e.l, paste, collapse="EU.TW")))

## GATK!!! readgroup to bamfile with picard 

## or alternatively (easier to parse)
## /home/ele/tools/samtools/samtools mpileup -uf
## ../fullest_assembly_imputed.fasta *.sorted.bam |
## /home/ele/tools/samtools/bcftools/bcftools view -gcv -

VCF <- read.delim("/home/ele/Data/454/mapping/mapping_each_lib/all_imp.vca",
                  skip=26, as.is=TRUE)

## does a e.g. tA point tow wrongly imputed sequence???
## wrong.imputed <- VCF[grepl("(a|c|g|t){1}(A|C|G|T){1}", VCF$ALT), ]
## wrong.imputed.contigs <- wrong.imputed$X.CHROM

## get only snps
VCF <- VCF[nchar(VCF$REF)<2 & nchar(VCF$ALT)<2, ]
VCF <- VCF[grepl("a|c|g|t", VCF$REF, ignore.case=TRUE) &
           grepl("a|c|g|t", VCF$ALT, ignore.case=TRUE), ]


super.geno <- function (VCF.obj, snp.qual, gt.qual, sum=TRUE){
  VCFQ <- VCF[as.numeric(VCF$QUAL)>snp.qual, ]
  VCFQ <- VCFQ[, c(1, 2, 4, 5, 10:15)]
  names(VCFQ) <- gsub("final.(\\d*|\\w*).sff.trimmed.sff.fasta.sorted.bam",
                      "\\1", names(VCFQ))

  ## as the library-format we have GT:PL:GQ

  ## GT genotype, encoded as alleles values separated , e.g. The allele
  ## values are 0 for the reference allele (what is in the reference
  ## sequence), 1 for the first allele listed in ALT, 2 for the second
  ## allele list in ALT and so on. For diploid calls examples could be
  ## 0/1 or 1|0 etc. For haploid calls, e.g. on Y, male X,
  ## mitochondrion, only one allele value should be given. All samples
  ## must have GT call information; if a call cannot be made for a
  ## sample at a must be specified for each missing allele in the GT
  ## field (for example ./. for a diploid). The meanings of the
  ## separators are: / : genotype unphased | : genotype phased

  ## GL : genotype likelihoods
  ## comprised of comma separated floating point log10-scaled
  ## likelihoods for all possible genotypes given the set of alleles
  ## defined in the REF and ALT fields. In presence of the GT field the
  ## same ploidy is expected and the canonical order is used; without GT
  ## field, diploidy is assumed. If A is the allele in REF and
  ## B,C,... are the alleles as ordered in ALT, the ordering of
  ## genotypes for the likelihoods is given by: F(j/k) = (k*(k+1)/2)+j.
  ## In other words, for biallelic sites the ordering is: AA,AB,BB; for
  ## triallelic sites the ordering is: AA,AB,BB,AC,BC,CC, etc.  For
  ## example: GT:GL 0/1:-323.03,-99.29,-802.53 (Numeric)

  ## PL : the phred-scaled genotype likelihoods rounded to the closest
  ## integer (and otherwise defined precisely as the GL field)
  ## (Integers)

  ## GQ : conditional genotype quality, encoded as a phred quality
  ## -10log_10p(genotype call is wrong, conditioned on the site's being
  ## variant) (numeric)

  get.likely.gt <- function (vcf.col, gt.qual=10){
    gt.l <- strsplit(vcf.col, ":|,")
    gt.l <- lapply(gt.l, function(x) x[c(1,5)])
    g <- lapply(gt.l, function (x){
      if (as.numeric(as.character(x[[2]]))>gt.qual) { return(x[1]) }
      else{ return(NA) }
    })
    unlist(g)
  }

  gt.l <- apply(VCFQ[,5:10], 2, get.likely.gt)
  gt.df <- as.data.frame( gt.l)
  gt.df$"L2R3" <- NULL
  
  ## summary(gt.df)

  V <- cbind(VCFQ[ , c(1:4)], gt.df)
  names(V)[1:4] <- c( "contig", "base", "from", "to")

  VAR.pop <- merge(VARqp, V, by=c( "contig", "base", "from", "to"))

  summary(VAR.pop[, names(VAR.pop)%in%names(V)])

  relative.h <- function (x){
    length(x[x%in%c("0/1")])/length(x[x%in%c("0/0", "1/1")])
  }

  if (sum) {return (apply(VAR.pop[, c("10F", "179F", "KS4F", "M175", "UW07F")],
                          2, relative.h))}
  else {return(VAR.pop)}
}

super.geno(VCF, 0, 30)
super.geno(VCF, 0, 40)
super.geno(VCF, 0, 50)
super.geno(VCF, 0, 60)

super.geno(VCF, 10, 30)
super.geno(VCF, 10, 40)
super.geno(VCF, 10, 50)
super.geno(VCF, 10, 60)

super.geno(VCF, 30, 30)
super.geno(VCF, 30, 40)
super.geno(VCF, 30, 50)

het.table <- as.data.frame(super.geno(VCF, 10, 10, TRUE))
VAR.pop <- super.geno(VCF, 10, 10, FALSE)

sep.gt <- function (gt.col){
  l <- strsplit(as.character(gt.col), "/")
  s <- lapply(l, function (x) {if (is.na(x[1])){ rep(x,2)}
              else{x}})
  unlist(s)}

Rh <- VAR.pop[, c("10F", "179F", "KS4F", "M175", "UW07F")]
### remove loci with only one allele
Rh <- Rh[ apply(Rh, 1, function (x) ("0/1"%in%x & "0/0"%in%x) | ("0/1"%in%x & "1/1"%in%x) ), ]

number.hetero.info.snps <- nrow(Rh)

Rh.data <- as.data.frame(rbind(sep.gt(Rh$"10F"),
                                sep.gt(Rh$"179F"),
                                sep.gt(Rh$"KS4F"),
                                sep.gt(Rh$"M175"),
                                sep.gt(Rh$"UW07F")))

library(Rhh)
ir.test <- hh(Rh.data, 100, "ir")
hl.test <- hh(Rh.data, 100, "hl")
hh.test <- hh(Rh.data, 100, "sh")

## mean(hl.test)
## quantile(hl.test, probs=c(0.025, 0.975))

## Internal relatedness, a multilocus heterozygosity measure developed
## by Amos et al. (2001): Amos W, Worthington Wilmer J, Fullard K et
## al (2001) The influence of parental relatedness on reproductive
## success. Proc R Soc Lond B 268:2021-2027

het.table <- cbind(het.table, ir(Rh.data))

## homozygosity by loci: Aparicio et al. (2007)?(2006) : Aparicio JM,
##  Ortego J, Cordero PJ (2006) What should we weigh to estimate
##  heterozygosity, alleles or loci? Mol Ecol 15:4659-4665

het.table <- cbind(het.table, hl(Rh.data))

## standardized heterozygosity: Coltman et al. (1999): Coltman DW,
## Pilkington JG, Smith JA et al (1999) Parasite-mediated selection
## against inbred Soay sheep in a free-living island
## population. Evolution 53:1259-1267

het.table <- cbind(het.table, sh(Rh.data))


f
o
o
o
o


## READ: pairwise synonymous diversity und Wattersons theta

## all.pop.VARfreq <- c(VARpop$perc.common, VARpop$perc.EU, VARpop$perc)
## all.pop.VARfreq <- as.data.frame(cbind(freq=all.pop.VARfreq, pop=rep(c("Common", "Europe", "Taiwan"), each=nrow(VARpop))))
## all.pop.VARfreq$freq <- as.numeric(as.character(all.pop.VARfreq$freq))

## png("popfreq.png")
## histogram(~freq |pop, data = all.pop.VARfreq,
##           breaks=100,
##           xlim=c(0,100),
##           panel=function (...) {
##             panel.grid(v=9, h=0)
##             panel.histogram(...)},
##           layout=c(1,3))
## dev.off()


## ## Treating every SNP as independent ... surely wrong
## venn.pop.snp <- as.data.frame(cbind(!is.na(VARpop$pop.TW), !is.na(VARpop$pop.EU)))
## names(venn.pop.snp) <- c("Taiwan", "Europe")
## vp <- venn.pop.snp[venn.pop.snp$Taiwan | venn.pop.snp$Europe,]
## ## vennDiagram(vp)

## venn.pop.snp <- cbind(contig=as.character(VARpop$contig), venn.pop.snp)

## ## exclude all contigs with not at least one snp in each population
## both <- c(by(venn.pop.snp, venn.pop.snp$contig, function (x) any(x$Europe)&any(x$Taiwan)))
## venn.pop.snp <- merge(venn.pop.snp, as.data.frame(both), by.x="contig", by.y="row.names")

## vp.1 <- venn.pop.snp[venn.pop.snp$both,2:3]
## ## vennDiagram(vp.1)

## ## Count EU/TW SNPs per contig
## TW.p.c.s <- by(VARpop, VARpop$contig, function (x) nrow(x[!is.na(x$pop.TW),]))
## EU.p.c.s <- by(VARpop, VARpop$contig, function (x) nrow(x[!is.na(x$pop.EU),]))
## Co.p.c.s <- by(VARpop, VARpop$contig, function (x) nrow(x[!is.na(x$pop.EU) & !is.na(x$pop.TW),]))

## ## Try normalizing by the coverage at snp positions
## TW.p.c.r <- by(VARpop, VARpop$contig, function (x) mean(x$nfrom+x$nto, na.rm=TRUE))
## EU.p.c.r <- by(VARpop, VARpop$contig, function (x) mean(x$nfrom.EU+x$nto.EU, na.rm=TRUE))


## snp.per.contig <- as.data.frame(cbind(TW=TW.p.c.s, EU=EU.p.c.s, Common=Co.p.c.s, TWr=TW.p.c.r, EUr=EU.p.c.r))
## snp.per.contig <- snp.per.contig[!apply(snp.per.contig, 1, function (x) any(is.na(x))),]

## ratio.e.pop.s <- sum(snp.per.contig$EU-snp.per.contig$Common)/sum(snp.per.contig$TW-snp.per.contig$Common)

## snp.per.screened <- snp.per.contig[!apply(snp.per.contig, 1, function (x) any(is.na(x))),]

## own.ratio <- sum((snp.per.screened$EU-snp.per.screened$Common)/snp.per.screened$EUr, na.rm=T)/sum((snp.per.screened$TW-snp.per.screened$Common)/snp.per.screened$TWr, na.rm=T)
