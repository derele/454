source("snp_analysis.R")

produce.population.mapping <- function (VARobj){
EUVAR <- read.var ("/home/ele/Data/454/mapping/mapping_each_population/europe_vs_fullest_imputed.varsnp")
EUVAR <- cbind(EUVAR, pop.EU=TRUE)

TWVAR <- read.var ("/home/ele/Data/454/mapping/mapping_each_population/taiwan_vs_fullest_imputed.varsnp")
TWVAR <- cbind(TWVAR, pop.TW=TRUE)

VARpop <- merge(VARobj, EUVAR, by=c("contig", "base", "from", "to"), all.x=TRUE, suffixes = c(".common", ".EU"))

VARpop <- merge(VARpop, TWVAR, by=c("contig", "base", "from", "to"), all.x=TRUE, suffixes =c (".EU", ".TW"))
return(VARpop)
}

VARpop <- produce.population.mapping(VAR)

## Lesen: pairwise synonymous diversity und Wattersons theta

all.pop.VARfreq <- c(VARpop$perc.common, VARpop$perc.EU, VARpop$perc)
all.pop.VARfreq <- as.data.frame(cbind(freq=all.pop.VARfreq, pop=rep(c("Common", "Europe", "Taiwan"), each=nrow(VARpop))))
all.pop.VARfreq$freq <- as.numeric(as.character(all.pop.VARfreq$freq))

png("popfreq.png")
histogram(~freq |pop, data = all.pop.VARfreq,
          breaks=100,
          xlim=c(0,100),
          panel=function (...) {
            panel.grid(v=9, h=0)
            panel.histogram(...)},
          layout=c(1,3))
dev.off()


## Treating every SNP as independent ... surely wrong
venn.pop.snp <- as.data.frame(cbind(!is.na(VARpop$pop.TW), !is.na(VARpop$pop.EU)))
names(venn.pop.snp) <- c("Taiwan", "Europe")
vp <- venn.pop.snp[venn.pop.snp$Taiwan | venn.pop.snp$Europe,]
## vennDiagram(vp)

venn.pop.snp <- cbind(contig=as.character(VARpop$contig), venn.pop.snp)

## exclude all contigs with not at least one snp in each population
both <- c(by(venn.pop.snp, venn.pop.snp$contig, function (x) any(x$Europe)&any(x$Taiwan)))
venn.pop.snp <- merge(venn.pop.snp, as.data.frame(both), by.x="contig", by.y="row.names")

vp.1 <- venn.pop.snp[venn.pop.snp$both,2:3]
## vennDiagram(vp.1)

## Count EU/TW SNPs per contig
TW.p.c.s <- by(VARpop, VARpop$contig, function (x) nrow(x[!is.na(x$pop.TW),]))
EU.p.c.s <- by(VARpop, VARpop$contig, function (x) nrow(x[!is.na(x$pop.EU),]))
Co.p.c.s <- by(VARpop, VARpop$contig, function (x) nrow(x[!is.na(x$pop.EU) & !is.na(x$pop.TW),]))

## Try normalizing by the coverage at snp positions
TW.p.c.r <- by(VARpop, VARpop$contig, function (x) mean(x$nfrom+x$nto, na.rm=TRUE))
EU.p.c.r <- by(VARpop, VARpop$contig, function (x) mean(x$nfrom.EU+x$nto.EU, na.rm=TRUE))


snp.per.contig <- as.data.frame(cbind(TW=TW.p.c.s, EU=EU.p.c.s, Common=Co.p.c.s, TWr=TW.p.c.r, EUr=EU.p.c.r))
snp.per.contig <- snp.per.contig[!apply(snp.per.contig, 1, function (x) any(is.na(x))),]

ratio.e.pop.s <- sum(snp.per.contig$EU-snp.per.contig$Common)/sum(snp.per.contig$TW-snp.per.contig$Common)

snp.per.screened <- snp.per.contig[!apply(snp.per.contig, 1, function (x) any(is.na(x))),]

own.ratio <- sum((snp.per.screened$EU-snp.per.screened$Common)/snp.per.screened$EUr, na.rm=T)/sum((snp.per.screened$TW-snp.per.screened$Common)/snp.per.screened$TWr, na.rm=T)
