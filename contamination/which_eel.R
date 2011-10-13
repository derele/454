CO <- read.delim("/home/ele/Data/454/which_eel/fish_reads.cobl")

## Get the library-ids for the reads
F10 <- readLines("/home/ele/Data/454/rRNA_screened/10F.sff.trimmed.acc")
F179 <- readLines("/home/ele/Data/454/rRNA_screened/179F.sff.trimmed.acc")
FKS4 <- readLines("/home/ele/Data/454/rRNA_screened/KS4F.sff.trimmed.acc")
FUW07 <- readLines("/home/ele/Data/454/rRNA_screened/UW07F.sff.trimmed.acc")
L2R3 <- readLines("/home/ele/Data/454/rRNA_screened/L2R3.sff.trimmed.acc")

C <- CO[CO$dbhit!="No_hit",]
C$bitscore <- as.numeric(C$bitscore)

require(MASS)

kde2dplot <- function(d,                # a 2d density computed by kde2D
                      ncol=50,          # the number of colors to use
                      zlim=c(0,max(z)), # limits in z coordinates
                      nlevels=20,
                      ...)       # see option nlevels in contour
		      {
z   <- d$z
nrz <- nrow(z)
ncz <- ncol(z)

couleurs  <- tail(topo.colors(trunc(1.4 * ncol)),ncol)

image(d,col=couleurs, ...)
contour(d,add=T,nlevels=nlevels)
box()
}


d.euro.blast <- kde2d(C[C$dbhit=="european_eel.fasta", "gc"],
                     log(C[C$dbhit=="european_eel.fasta","bitscore"]))

d.tw.blast <- kde2d(C[C$dbhit=="Li_68_ContigsSinglets.fasta", "gc"],
                     log(C[C$dbhit=="Li_68_ContigsSinglets.fasta", "bitscore"]))


par(mfrow=c(2,1))
kde2dplot(d.euro.blast)
kde2dplot(d.tw.blast)


summary(C[C$seq_name%in%F10, "dbhit"])
summary(C[C$seq_name%in%F179, "dbhit"])

## This looks like FKS4 could be taiwanese as well
summary(C[C$seq_name%in%FKS4, "dbhit"])

summary(C[C$seq_name%in%FUW07, "dbhit"])
summary(C[C$seq_name%in%L2R3, "dbhit"])

## Now with higher bit-scores:

summary(C[C$seq_name%in%F10 & C$bitscore>100, "dbhit"])
summary(C[C$seq_name%in%F179  & C$bitscore>100, "dbhit"])

## This looks like FKS4 could be taiwanese as well
summary(C[C$seq_name%in%FKS4  & C$bitscore>100, "dbhit"])

summary(C[C$seq_name%in%FUW07  & C$bitscore>100, "dbhit"])
summary(C[C$seq_name%in%L2R3  & C$bitscore>100, "dbhit"])

## Now reading rather the full cobl-table:

A <- read.delim("/home/ele/Data/454/which_eel/all_crassus.cobl")
Af <- A[grepl("rRNA", A$dbhit),]
A.rRNA <- Af[(grepl("Anguilla", Af$annot)), ]
A.rRNA <- A.rRNA[as.numeric(as.character(A.rRNA$evalue))<1e-100,]

summary(A.rRNA[A.rRNA$seq_name%in%F10, "annot"], maxsum = 5)
summary(A.rRNA[A.rRNA$seq_name%in%F179, "annot"], maxsum = 5)

## This looking like FKS4 Taiwanese as well
summary(A.rRNA[A.rRNA$seq_name%in%FKS4, "annot"], maxsum = 5)

summary(A.rRNA[A.rRNA$seq_name%in%FUW07, "annot"], maxsum = 5)
summary(A.rRNA[A.rRNA$seq_name%in%L2R3, "annot"], maxsum = 5)



KS4_tmp <- read.delim("/home/ele/Data/454/trimmed/tmp/KS4.cobl")

summary(KS4_tmp$dbhit)

KS4_tmp<- KS4_tmp[grepl("rRNA", KS4_tmp$dbhit),]
KS4_tmp <- KS4_tmp[(grepl("Anguilla", KS4_tmp$annot)), ]

summary(KS4_tmp$annot, maxsum = 15)



## Not possible to decide on 18S, as having same bit scores 
## F <- read.delim("/home/ele/Data/454/which_eel/fish_reads_18S.cobl")
## F <- F[!is.na(as.numeric(as.character(F$evalue))),]

## summary(F[F$seq_name%in%F10, "annot"], maxsum = 5)
## summary(F[F$seq_name%in%F179, "annot"], maxsum = 5)

## ## This looking like FKS4 Taiwanese as well
## summary(F[F$seq_name%in%FKS4, "annot"], maxsum = 5)

## summary(F[F$seq_name%in%FUW07, "annot"], maxsum = 5)
## summary(F[F$seq_name%in%L2R3, "annot"], maxsum = 5)

