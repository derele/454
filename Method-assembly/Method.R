###################################################
### chunk number 1: load
###################################################
#line 20 "/home/ele/thesis/454/Method-assembly/Method.Rnw"
library(lattice)
library(ggplot2)
library(MASS)
library(xtable)
library(limma)
library(Design)


###################################################
### chunk number 2: functions_read
###################################################
#line 30 "/home/ele/thesis/454/Method-assembly/Method.Rnw"


read.sequences <- function (file){
  content <- readLines(file)
  header.lines <- grep( "^>", content)
  start.lines <- header.lines+1
  end.lines <- c(header.lines[-1]-1, length(content))
  sq <- sapply(1:length(start.lines), function (i) {
    list(content[start.lines[i]:end.lines[i]])})
  names(sq) <- substr(content[header.lines],2, nchar(content[header.lines]))
  sq <- unlist(lapply(sq, paste, collapse=""))
  return(sq[nchar(sq)>0])
}

mean.qual <- function(q){
  lapply(q, function(x) {w <- unlist(strsplit(x, " "))
                                    mean(as.numeric(w[nchar(w)>0]))})
}

raw.fasta <- read.sequences("/home/ele/Data/454/assemblies/mira/Acrassus_in.454.fasta")

raw.qual <- read.sequences("/home/ele/Data/454/assemblies/mira/Acrassus_in.454.fasta.qual")

qualmeans <- mean.qual(raw.qual)

ids <- substr(names(raw.fasta),1, 14)
seqlength <- nchar(raw.fasta)

RAW <- as.data.frame(cbind(ids=as.character(ids),
                           qualmeans=as.numeric(qualmeans),
                           seqlength=as.numeric(seqlength)))

nFLX <- length(grep("^F.*",perl=T,names(raw.fasta),value=T))
nTIT <- length(grep("^G.*",perl=T,names(raw.fasta), value=T))

## raw.fasta contains all trimmed bases in lower-case
trimmed.fasta <- gsub("[a-z]", "", raw.fasta)
nBases <- sum(nchar(trimmed.fasta))



###################################################
### chunk number 3: mirace
###################################################
#line 86 "/home/ele/thesis/454/Method-assembly/Method.Rnw"

read.ace <- function (path){
  command <- paste("perl -ne '$c =  $1 if /^CO (\\S+)/;",
                   "print \"$c $1\\n\"if /^AF (\\w+)(\\|*|\\.*)/'",
                   path)                    # last bit to get rid of
                                            # newblers .info and my
                                            # own |readcount
#  closeAllConnections()
  ace <- read.table(pipe(command), sep=" ")
  names(ace) <- c("contig", "read")
  return(ace)
}

MIR <- read.ace("/home/ele/Data/454/assemblies/mira/Acrassus_assembly/Acrassus_d_results/Acrassus_out.ace")


###################################################
### chunk number 4: newace
###################################################
#line 106 "/home/ele/thesis/454/Method-assembly/Method.Rnw"
  
N25 <- read.ace("/home/ele/Data/454/assemblies/newbler/P_2011_06_01_09_13_23_runAssembly/454Isotigs.ace")
## The newbler Singletons without the addition of reads
## in contigs not present in the fasta
wrongNew <- merge(N25, RAW, by.x="read", by.y="ids", all=TRUE)
wrongNewblerSing <- unique(wrongNew[is.na(wrongNew$contig),]$read)

CNewbler <- read.sequences("/home/ele/Data/454/assemblies/newbler/P_2011_06_01_09_13_23_runAssembly/454Isotigs.fna")
Newblercontig.names <- substr(names(CNewbler), 1, 11)
N25 <- N25[N25$contig%in%Newblercontig.names,]

## Merge the data-frame to get our goal data-strukture
MIR.merged <- merge(MIR, RAW, by.x="read", by.y="ids", all=TRUE)
MERGED <- merge(N25, MIR.merged, by="read", all=TRUE)

MiraSing <- unique(MERGED[is.na(MERGED$contig.y),]$read)
NewblerSing <- unique(MERGED[is.na(MERGED$contig.x),]$read)

contigs <- read.sequences("/home/ele/Data/454/assemblies/second_order/first_order.fasta.cap.contigs")
contigs <- c(contigs, read.sequences("/home/ele/Data/454/assemblies/second_order/first_order.fasta.cap.singlets"))


###################################################
### chunk number 5: newdist
###################################################
#line 152 "/home/ele/thesis/454/Method-assembly/Method.Rnw"
## Plot the distribution of read-splitting occuring in Newbler ###
## Numbers of reads per Newbler-contig
## are visble by the amount a read-name is repeated

n25d <- as.data.frame(table(N25$read))
names(n25d)[2] <- "new.split"

MERGED <- merge(MERGED, n25d, by.x="read",
                by.y="Var1", all=TRUE)

n25.dist <- ggplot(MERGED, aes(x=new.split, y=..count..)) +
  geom_bar(binwidth=1, color="white") +
  scale_y_log10()



###################################################
### chunk number 6: pc
###################################################
#line 195 "/home/ele/thesis/454/Method-assembly/Method.Rnw"
pc <- read.delim("/home/ele/Data/454/assemblies/pc/contig_stats.txt")[, -c(24:30)]
pc$Filename <- gsub("\\/.*", "", pc$Filename)
names(pc)[1:9] <- c("assembly", "Max length", "Number of contigs", "Number of Bases", "N50", "Number of congtigs in N50", "mean GC", "non ATGC bases", "Mean length") 

P <- as.data.frame(t(pc[,c(2:9)]))
names(P) <- pc$assembly
xtable(P)


###################################################
### chunk number 7: 
###################################################
#line 214 "/home/ele/thesis/454/Method-assembly/Method.Rnw"

Bm454 <- read.delim("/home/ele/Data/454/assemblies/blast_eval/454Isotigs_vs_Bm.covred", header=FALSE)
Bm454$assembler <- "newbler"
Bm454$species <- "Bm"

Ce454 <- read.delim("/home/ele/Data/454/assemblies/blast_eval/454Isotigs_vs_Ce.covred", header=FALSE)
Ce454$assembler <- "newbler"
Ce454$species <- "Ce"

Bmmir <- read.delim("/home/ele/Data/454/assemblies/blast_eval/Acrassus_vs_Bm.covred", header=FALSE)
Bmmir$assembler <- "mira"
Bmmir$species <- "Bm"

Cemir <- read.delim("/home/ele/Data/454/assemblies/blast_eval/Acrassus_vs_Ce.covred", header=FALSE)
Cemir$assembler <- "mira"
Cemir$species <- "Ce"

CeSndO <- read.delim("/home/ele/Data/454/assemblies/blast_eval/SndO_vs_Ce.covred", header=FALSE)
CeSndO$assembler <- "SndO"
CeSndO$species <- "Ce"

BmSndO <- read.delim("/home/ele/Data/454/assemblies/blast_eval/SndO_vs_Bm.covred", header=FALSE)
BmSndO$assembler <- "SndO"
BmSndO$species <- "Bm"

BmMN <- read.delim("/home/ele/Data/454/assemblies/blast_eval/MN_vs_Bm.covred", header=FALSE)
BmMN$assembler <- "SndO.MN"
BmMN$species <- "Bm"

CeMN <- read.delim("/home/ele/Data/454/assemblies/blast_eval/MN_vs_Ce.covred", header=FALSE)
CeMN$assembler <- "SndO.MN"
CeMN$species <- "Ce"


BE <- rbind(Bm454, Ce454, Bmmir, Cemir, BmSndO, CeSndO, BmMN, CeMN)[,-c(1,2)]
names(BE)[1:3] <- c("uncovered", "covered", "CovRed")

suj.summary <- function (x){
  w <- sum(x[,1])/sum(c(x[,1], x[,2]))*100
  y <- nrow(x[x[,1]>0,])/nrow(x)*100
  z <- sum(x[x[,1]>0,3])
  c(percbaseCov=w, percprotCov=y,  sumCovRed=z)
}

suj.sum <- do.call("rbind", by(BE, as.factor(BE$species):as.factor(BE$assembler), suj.summary))
xtable(suj.sum)


###################################################
### chunk number 8: sing
###################################################
#line 268 "/home/ele/thesis/454/Method-assembly/Method.Rnw"
frass <- unique(MERGED[is.na(MERGED$contig.x) & is.na(MERGED$contig.y),]$read)
frass.fasta <- trimmed.fasta[names(trimmed.fasta)%in%frass]


###################################################
### chunk number 9: capace
###################################################
#line 273 "/home/ele/thesis/454/Method-assembly/Method.Rnw"
SnO <- read.ace("/home/ele/Data/454/assemblies/second_order/first_order.fasta.cap.ace")
MERGED <- merge(MERGED, SnO, by.x="contig.x", by.y="read", all.x=TRUE)
names(MERGED)[ncol(MERGED)] <- "SndO.NewblerContig"

MERGED <- merge(MERGED, SnO, by.x="contig.y", by.y="read", all.x=TRUE)
names(MERGED)[ncol(MERGED)] <- "SndO.MiraContig"


###################################################
### chunk number 10: categorize
###################################################
#line 307 "/home/ele/thesis/454/Method-assembly/Method.Rnw"

MN <- c(as.character(MERGED[MERGED$SndO.MiraContig%in%MERGED$SndO.NewblerContig,]$SndO.MiraContig),
        as.character(MERGED[MERGED$SndO.NewblerContig%in%MERGED$SndO.MiraContig,]$SndO.NewblerContig))
MN <- unique(MN[!is.na(MN)])

M_n <- MERGED[!MERGED$SndO.MiraContig%in%MERGED$SndO.NewblerContig,]$SndO.MiraContig
M_n <- unique(M_n[!is.na(M_n)])

N_n <- MERGED[!MERGED$SndO.NewblerContig%in%MERGED$SndO.MiraContig,]$SndO.NewblerContig
N_n <- unique(N_n[!is.na(N_n)])

N_1<- by(MERGED, MERGED$contig.x,
         function (x) all(is.na(x$SndO.NewblerContig)))
N_1 <- names(N_1[N_1 & !is.na(N_1)])

M_1<- by(MERGED, MERGED$contig.y,
         function (x) all(is.na(x$SndO.MiraContig)))
M_1 <- names(M_1[M_1 & !is.na(M_1)])

M_n.contigs <- unique(MERGED[MERGED$SndO.MiraContig%in%M_n,]$contig.y)

N_n.contigs <- unique(MERGED[MERGED$SndO.NewblerContig%in%N_n,]$contig.x)

MN.newbler <- unique(MERGED[MERGED$SndO.NewblerContig%in%MN, ]$contig.x)

MN.mira <- unique(MERGED[MERGED$SndO.MiraContig%in%MN, ]$contig.y)

## need unique here to get no newbler duplicates
MN.reads.mira <- unique(MERGED[MERGED$SndO.MiraContig%in%MN ,]$read)

## no unique here to get also the dplicates
MN.reads.newbler <- MERGED[MERGED$SndO.NewblerContig%in%MN ,]$read

## also here remove reads duplicated by newbler for  mira reads
MN.reads.inc <- MERGED[(MERGED$read%in%MN.reads.mira & !duplicated(MERGED$read)) |
                   (MERGED$read%in%MN.reads.newbler & !duplicated(MERGED$read)) ,]$read

MN.reads.both <- MERGED[MERGED$read%in%MN.reads.mira  & !duplicated(MERGED$read) &
                   MERGED$read%in%MN.reads.newbler ,]$read

N_n.reads <- unique(MERGED[MERGED$SndO.NewblerContig%in%N_n,]$read)

M_n.reads <- unique(MERGED[MERGED$SndO.MiraContig%in%M_n,]$read)

M_1.reads <- unique(MERGED[MERGED$contig.y%in%M_1,]$read)

N_1.reads <- unique(MERGED[MERGED$contig.x%in%N_1,]$read)


###################################################
### chunk number 11: tabulate
###################################################
#line 358 "/home/ele/thesis/454/Method-assembly/Method.Rnw"
SndO.con <- cbind("M_1"=" ",
                   "M_n"=length(M_n),
                   "MN"=length(MN),
                   "N_n"=length(N_n),
                   "N_1"=" ")

FstO.con <- cbind(length(M_1),
                  length(M_n.contigs),
                  paste("mira=", length(MN.mira), "/",
                        "newbler=", length(MN.newbler), sep=""),
                  length(N_n.contigs),
                  length(N_1))

reads <- cbind(length(M_1.reads),
                 length(M_n.reads),
                 length(MN.reads.inc),
                 unique(length(MN.reads.both)),
                 length(N_n.reads),
                 length(N_1.reads))

reads <- c(reads[c(1,2)],
           paste("one=", reads[3], "/", "both=", reads[4], sep=""),
           reads[c(5,6)])


MNtable <- rbind(SndO.con, FstO.con, reads)
 
row.names(MNtable) <- c("Snd.o.con", "Fst.o.con", "reads")

xtable(MNtable, caption="\\small{\\textbf{Number of reads, first-order contigs (Fst.o.con) and second-order contigs (Snd.o.con) for different categories of contigs (M\\_1 and N\\_1 = first-order contigs not assembled in second-order assembly, from mira and newbler respectively; M\\_n and N\\_n =  assembled in second-order contigs only with contigs from the same first-order assembly; MN = assembled in second-order contigs with first order contigs from both first order assemblies}}", label="tab:categ")



###################################################
### chunk number 12: convenn
###################################################
#line 400 "/home/ele/thesis/454/Method-assembly/Method.Rnw"

vennData <- cbind(Mira_in_MN=MERGED$read%in%MN.reads.mira,
                  duplicated=duplicated(MERGED$read)&
                  (!MERGED$read%in%MN.reads.mira|MERGED$read%in%MN.reads.newbler),
                  Newbler_in_MN=MERGED$read%in%MN.reads.newbler)

#vennDiagram(vennData, cex=0.9)



###################################################
### chunk number 13: contigdist
###################################################
#line 425 "/home/ele/thesis/454/Method-assembly/Method.Rnw"
mira.contigs.per.SndO <- by(MERGED, MERGED[,"SndO.MiraContig"], function (x) length(unique(x[!is.na(x$contig.y), ]$contig.y)))
mira.contigs.per.SndO[is.na(mira.contigs.per.SndO)] <- 0
mira.contigs.per.SndO <- as.data.frame(unlist(list(mira.contigs.per.SndO)))

newbler.contigs.per.SndO <- by(MERGED, MERGED[,"SndO.NewblerContig"], function (x) length(unique(x[!is.na(x$contig.x),]$contig.x)))
newbler.contigs.per.SndO[is.na(newbler.contigs.per.SndO)] <- 0
newbler.contigs.per.SndO <- as.data.frame(unlist(list(newbler.contigs.per.SndO)))

contigs.per.SndO <- merge(newbler.contigs.per.SndO, mira.contigs.per.SndO,by="row.names")
names(contigs.per.SndO) <- c("contig", "Newbler_contigs", "Mira_contigs")



###################################################
### chunk number 14: readdist
###################################################
#line 439 "/home/ele/thesis/454/Method-assembly/Method.Rnw"
mira.reads.per.SndO <- by(MERGED, MERGED[,"SndO.MiraContig"], function (x) length(unique(x[!is.na(x$contig.y), ]$read)))
mira.reads.per.SndO[is.na(mira.reads.per.SndO)] <- 0
mira.reads.per.SndO <- as.data.frame(unlist(list(mira.reads.per.SndO)))

newbler.reads.per.SndO <- by(MERGED, MERGED[,"SndO.NewblerContig"], function (x) length(unique(x[!is.na(x$contig.x),]$read)))
newbler.reads.per.SndO[is.na(newbler.reads.per.SndO)] <- 0
newbler.reads.per.SndO <- as.data.frame(unlist(list(newbler.reads.per.SndO)))

reads.per.SndO <- merge(newbler.reads.per.SndO, mira.reads.per.SndO, by="row.names")
names(reads.per.SndO) <- c("contig", "reads_through_Newbler", "reads_through_Mira")

contig.df <- merge(as.data.frame(contigs), reads.per.SndO, by.x="row.names", by.y="contig")
names(contig.df)[1:2] <- c("contig", "seq")
contig.df <- merge(contig.df, contigs.per.SndO)



###################################################
### chunk number 15: readconplot
###################################################
#line 466 "/home/ele/thesis/454/Method-assembly/Method.Rnw"

c.p.SndO.plot <- ggplot(contig.df,
                        aes(x=Newbler_contigs, y=Mira_contigs)) +
        geom_jitter(alpha=0.5) +
        scale_y_continuous(limits=c(0, 20))


r.p.SndO.plot <- ggplot(contig.df,
                        aes(x=reads_through_Newbler, y=reads_through_Mira)) +
       geom_point() +
       scale_x_log10() +
       scale_y_log10()



###################################################
### chunk number 16: 
###################################################
#line 523 "/home/ele/thesis/454/Method-assembly/Method.Rnw"
##  number of reads being split in first order assembly
num.new.split <- do.call("rbind",
                     as.list(by(MERGED, MERGED$SndO.NewblerContig,
                                function (x) nrow(x[x$new.split>1, ]))))
## The NA's, where there is no newbler contig in the SnO.contig can be replace by 0
num.new.split[is.na(num.new.split)] <- 0
contig.df <- merge(contig.df, num.new.split, by.x="contig", by.y="row.names")
names(contig.df)[ncol(contig.df)] <- "num.new.split"

## per SndOcontig maximal number of N25 contigs read is split into
max.new.split <- do.call("rbind",
                     as.list(by(MERGED, MERGED$SndO.NewblerContig,
                                function (x) max(x$new.split))))
## The NA's, where there is no newbler contig in the SnO.contig can be replace by 0
max.new.split[is.na(max.new.split)] <- 0
contig.df <- merge(contig.df, max.new.split, by.x="contig", by.y="row.names")
names(contig.df)[ncol(contig.df)] <- "max.new.split"


## the number of reads merged back into one second order contig
num.SndO.merge <- do.call("rbind",
                          as.list(by(MERGED, MERGED$SndO.NewblerContig,
                                     function (x) {
                                       nrow(x[x$SndO.MiraContig!=x$SndO.NewblerContig
                                              & !is.na(x$SndO.MiraContig)
                                              & !is.na(x$SndO.NewblerContig),])
                                     })))
num.SndO.merge[is.na(num.SndO.merge)] <- 0
contig.df <- merge(contig.df, num.SndO.merge, by.x="contig", by.y="row.names")
names(contig.df)[ncol(contig.df)] <- "num.SndO.merge"




###################################################
### chunk number 17: cluster
###################################################
#line 562 "/home/ele/thesis/454/Method-assembly/Method.Rnw"
SndO <- cbind(as.character(MERGED$SndO.NewblerContig), as.character(MERGED$SndO.MiraContig))
SndO <- SndO[SndO[,1]!=SndO[,2] & !is.na(SndO[,1]) & !is.na(SndO[,2]),]
a <- apply(SndO, 1 , list)

add.up.clusters <- function (clustlist){
  clustlist <- unique(clustlist)
  clus <- list()
  for(k in 1:length(clustlist)){
    e <- sapply(1:length(clustlist),function (i) {
      any(unlist(clustlist[k]) %in% unlist(clustlist[i]))
    })
    clus[[k]] <- unique(unlist(clustlist[e]))
  }
  return(unique(clus))
}

## make it work recursively to fully collapse clusters
recursive.cluster <- function(clusterlist){
  res <- list(NA, clusterlist)
  while(length(res[[length(res)]])!=length(res[[length(res)-1]])){
    res[[length(res)+1]] <- add.up.clusters(res[[length(res)]])
  }
  return(res[[length(res)]])
}

SndO.clusters <- recursive.cluster(a)
write(unlist(lapply(SndO.clusters, paste, collapse=" ")), file="SndO.clusters")


SndO.cluster.size <- as.data.frame(do.call("rbind",
                                           lapply(SndO.clusters,
                                                  function (x)
                                                  cbind(contig=x, cluster.size=length(x)))))
SndO.cluster.size$cluster.size <- as.numeric(as.character(SndO.cluster.size$cluster.size))

contig.df <- merge(contig.df, SndO.cluster.size,
                   by.x="contig", by.y="contig",
                   all=TRUE)

contig.df$cluster.size <- ifelse(is.na(contig.df$cluster.size), 1, contig.df$cluster.size)



###################################################
### chunk number 18: blasts
###################################################
#line 607 "/home/ele/thesis/454/Method-assembly/Method.Rnw"

CEblast <- read.delim("/home/ele/Data/454/assemblies/blast_eval/SndO_vs_Ce.blt", header=F)
BMblast <- read.delim("/home/ele/Data/454/assemblies/blast_eval/SndO_vs_Bm.blt", header=F)
## get best hits
CEblast <- CEblast[unique(CEblast$V1),]
BMblast <- BMblast[unique(BMblast$V1),]

BM <- merge(as.data.frame(MN), BMblast, by.x="MN", by.y="row.names", all.x=TRUE)

## isoforms?
Ce <- as.data.frame(cbind(row.names(CEblast), grepl("\\.\\d+\\w", CEblast[,2])))
names(Ce) <- c("contig", "spliced")
Ce <- merge(as.data.frame(MN), Ce, by.x="MN", by.y="contig", all.x=TRUE)
names(Ce) <- c("contig", "Ce.spliced")



###################################################
### chunk number 19: 
###################################################
#line 625 "/home/ele/thesis/454/Method-assembly/Method.Rnw"

## hit.model <- glm(!is.na(Ce$Ce.spliced) ~ Ce$max.split + Ce$num.SndO.merge +  Ce$max.SndO.split + Ce$cluster.size, family=binomial(link="logit" ))
## hit.model2 <- step(hit.model)

## splice.model <- glm(as.logical(Ce$Ce.spliced) ~ Ce$max.split + Ce$num.SndO.merge +  Ce$max.SndO.split + Ce$cluster.size, family=binomial(link="logit" ))
## hit.model2 <- step(hit.model)





###################################################
### chunk number 20: 
###################################################
#line 638 "/home/ele/thesis/454/Method-assembly/Method.Rnw"
## save what is needed downstream
save(file="Method.Rdata", list=c("contig.df", "MN", "M_n", "N_n", "M_1", "N_1", "frass", "read.sequences"))
save.image("Method_image.Rdata")
write(MN, "/home/ele/thesis/454/Method-assembly/MN.acc")
write(paste(">", names(frass.fasta), "\n" , frass.fasta, sep=""),
      "/home/ele/thesis/454/Method-assembly/frass.fasta")


