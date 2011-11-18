###################################################
### chunk number 1: load.libs
###################################################
#line 476 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
library(VennDiagram)
library(xtable)
library(reshape)
library(ggplot2)

library(GSEABase)
library(GOstats)
library(GO.db)

library(limma)
library(XML)
library(coin)
library(multcomp)

library(DESeq)
library(ShortRead)
library(biomaRt)

library(BioIDMapper)

library(Rhh)

source("/home/ele/thesis/454/common_R_functions.R")


###################################################
### chunk number 2: trim
###################################################
#line 517 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
trimdat <- readLines("/home/ele/Data/454/trimmed/trimming.stats")
trimfile <-  grep("Trimming sff file *", trimdat, value=TRUE)
trimfile <- sapply(strsplit(trimfile, " |,"), function (x) basename(x[4]))
lib <- gsub(".sff'", "", trimfile)

raw_reads<- grep("Sequences analyzed:.*", trimdat, value=TRUE)
raw_reads <- as.numeric(gsub("Sequences analyzed: *", "", raw_reads))
short <- grep("by 'short':.*", trimdat, value=TRUE)
short <- as.numeric(gsub(" *by 'short': *", "", short))
lowq <- grep("by 'low_qual':.*", trimdat, value=TRUE)
lowq <- as.numeric(gsub(" *by 'low_qual': *", "", lowq))
dust <- grep("by 'dust':.*", trimdat, value=TRUE)
dust <- as.numeric(gsub(" *by 'dust': *", "", dust))
shortq <- grep("by 'shortq':.*", trimdat, value=TRUE)
shortq <- as.numeric(gsub(" *by 'shortq': *", "", shortq))
valid <- raw_reads-(short+lowq+dust+shortq)

TRIM <- data.frame(row.names=lib, raw_reads, short, lowq, dust, shortq, valid)
TRIM["total",] <- colSums(TRIM)

## subsitute the publication style library abrevations
rownames(TRIM) <- gsub("179F", "T1", rownames(TRIM))
rownames(TRIM) <- gsub("10F", "T2", rownames(TRIM))
rownames(TRIM) <- gsub("KS4F", "E1", rownames(TRIM))
rownames(TRIM) <- gsub("UW07F", "E2", rownames(TRIM))
rownames(TRIM) <- gsub("L2R3", "L2", rownames(TRIM))
rownames(TRIM) <- gsub("M175", "M", rownames(TRIM))



###################################################
### chunk number 3: read.libs
###################################################
#line 548 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

RE <- list()
for (l in (c("10F", "179F", "M175", "L2R3", "KS4F", "UW07F"))){
   command <- paste("/home/ele/tools/Newbler2.6/bin/sffinfo -s /home/ele/Data/454/trimmed/",
                    l, ".sff.trimmed.sff;", sep="")
   connect <- pipe(command)
   RE[l] <- list(read.sequences(connect))
   closeAllConnections()
}
names <- gsub(" .*$", "", unlist(lapply(RE, names)))
lib <- melt(RE)

## subsitute the publication style library abrevations
lib$L1 <- gsub("179F", "T1", lib$L1)
lib$L1 <- gsub("10F", "T2", lib$L1)
lib$L1 <- gsub("KS4F", "E1", lib$L1)
lib$L1 <- gsub("UW07F", "E2", lib$L1)
lib$L1 <- gsub("L2R3", "L2", lib$L1)
lib$L1 <- gsub("M175", "M", lib$L1)

RE <- cbind(read=names, lib)
names(RE)[2:3] <- c("seq", "lib")
RE$seq <- as.character(RE$seq)
RE$length <- nchar(RE$seq)

eel.all <- read.sequences(pipe("/home/ele/tools/Newbler2.6/bin/sffinfo -s /home/ele/Data/454/Li68_assembly/sff/Li68.sff"))

eel.clean <- read.sequences(pipe("/home/ele/tools/Newbler2.6/bin/sffinfo -s /home/ele/Data/454/Li68_assembly/sff/Li68.sff.trimmed.sff"))

eel.rRNA.contigs <- read.sequences("/home/ele/Data/454/Li68_assembly/Aj_liver_rRNA.fasta")

eel.contigs <- read.sequences("/home/ele/Data/454/Li68_assembly/Aj_liver.fasta")



###################################################
### chunk number 4: load.cobl
###################################################
#line 584 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

eelmRNA <- read.blast.best("/home/ele/Data/454/pre_assembly_screening/trimmed_all_vs_eelmRNA.blt")
eelmRNA$dbhit <- "eelmRNA"

eelrRNA <- read.blast.best("/home/ele/Data/454/pre_assembly_screening/trimmed_all_vs_eelrRNA.blt")
eelrRNA$dbhit <- "eelrRNA"

AcrRNA <- read.blast.best("/home/ele/Data/454/pre_assembly_screening/trimmed_all_vs_AcrRNA.blt")
AcrRNA$dbhit <- "AcrRNA"
AcrRNA[grepl("Cerco|Flag", AcrRNA$V2), "dbhit"] <- "Cercozoa"

B <- rbind(eelmRNA, eelrRNA, AcrRNA)

B <- merge(B[, c(1, 3, 7, 8, 12, 13)], RE[,c(1,3,4)], 
           by.x="V1", by.y="read")

names(B)[1:5] <- c("read", "iden", "s.hit", "e.hit", "bits")

Br <- reduce.blast.by.length(B, len=0.8, iden=95)
R <- merge(Br, RE, all.y=TRUE)

## no dbhit means still valid
R$dbhit[is.na(R$dbhit)] <- "valid"

R$lib <- as.factor(R$lib)
R$dbhit <- as.factor(R$dbhit)

rRNA.lib.reads <- ggplot(R, aes(lib, ..count.., fill=dbhit)) +
  geom_bar() +
  coord_flip() +
  scale_y_continuous("Number of Reads") +
  scale_x_discrete("Library")

rRNA.lib.bases <- ggplot(R, aes(lib, ..count.., fill=dbhit, weight=length)) +
  geom_bar() +
  coord_flip() +
  scale_y_continuous("Number of Bases") +
  scale_x_discrete("Library")




###################################################
### chunk number 5: plot.rRNA
###################################################
#line 627 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
png("../figures/rRNA_plots.png", width=1000, height=1000, res=144)

# Set up the page
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)

# Make each plot, in the correct location
print(rRNA.lib.reads, vp = vplayout(1, 1 ))
print(rRNA.lib.bases, vp = vplayout(2, 1 ))
dev.off() 



###################################################
### chunk number 6: count.screened
###################################################
#line 643 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

pre.screen.tab <- as.data.frame.matrix(table(R$lib, R$dbhit))
pre.screen.span <- rbind( by(R, list(R$lib, R$dbhit), function (x) sum(nchar(x$seq))))

pre.screen.tab <- cbind(pre.screen.tab, valid.span=pre.screen.span[,"valid"])

nscreened <- nrow(R[R$dbhit%in%c("eelrRNA", "AcrRNA", "eelmRNA"), ])
nscreened.r <- nrow(R[R$dbhit%in%c("eelrRNA", "AcrRNA"),])
nscreened.m <- nrow(R[R$dbhit%in%"eelmRNA",])

AcrRNA_vs_nt <- read.blast.best("/home/ele/Data/454/pre_assembly_screening/AcrRNA_vs_nt.blt")
AcrRNA_cerco <- AcrRNA_vs_nt[grepl("Cerco|Flag", AcrRNA_vs_nt$V1), ]

cerco <- cbind( sequence.identifier=as.character(AcrRNA_cerco$V2), 
                sequence.identitiy=AcrRNA_cerco$V3,
                hsp.length=AcrRNA_cerco$V8-AcrRNA_cerco$V7)






###################################################
### chunk number 7: load.meth
###################################################
#line 688 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
load("/home/ele/thesis/454/Method-assembly/Method.Rdata")

nMN <- nrow(contig.df[contig.df$category%in%"MN",])

nbad.c <- nrow(contig.df[contig.df$category%in%
                         c("M_n","N_n","N_1","M_1"),])

nSing <- nrow(contig.df[contig.df$category%in%0,])


###################################################
### chunk number 8: mean.cov
###################################################
#line 720 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

cov.mn <- tapply(con.pile$coverage,
                 con.pile$contig%in%contig.df[contig.df$category=="MN", "contig"],
                 mean)
                            
cov.mn.u <- tapply(con.pile.uniq$uniq_coverage,
                   con.pile.uniq$contig%in%contig.df[contig.df$category=="MN", "contig"],
                   mean)

cov.nu.total <- mean(con.pile$coverage)
cov.u.total <- mean(con.pile.uniq$uniq_coverage)



###################################################
### chunk number 9: host.screen
###################################################
#line 735 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

## new style uncompromising screening 
eelmRNA <- read.blast.best("/home/ele/Data/454/post_assembly_screening/fullest_vs_eelmRNA.blt")
eelmRNA$dbhit <- "eelmRNA"

eelrRNA <- read.blast.best("/home/ele/Data/454/post_assembly_screening/fullest_vs_eelrRNA.blt")
eelrRNA$dbhit <- "eelrRNA"

AcrRNA <- read.blast.best("/home/ele/Data/454/post_assembly_screening/fullest_vs_AcrRNA.blt")
AcrRNA$dbhit <- "AcrRNA"

nempep4.nuc <- read.blast.best("/home/ele/Data/454/post_assembly_screening/fullest_vs_nembase_nuc.blt")
nempep4.nuc$dbhit <- "nempep4"

## looking at best across dbs using nempep to rescue 
E <- rbind(eelmRNA, eelrRNA, AcrRNA, nempep4.nuc)
E <- E[, c(1, 3, 7, 8, 12, 13)]
names(E) <- c("contig", "iden", "s.hit", "e.hit", "bits", "contamination")

E <- E[order(E$contig,E$bits, decreasing=TRUE), ]
E <- E[!duplicated(E$contig),]

## get the length of the contigs
E <- merge(E,
           contig.df[, names(contig.df)%in%c("contig", "length")])
## do the reduction based on length
reduce.blast.by.length <- function (bt, len, iden){
bt[((bt$length*len) < (bt$e.hit-bt$s.hit) & bt$iden>iden), ]
}
ER <- reduce.blast.by.length(E, len=0.5, iden=70)

contig.df <- merge(contig.df, ER[ , c("contig", "contamination")], all.x=TRUE)
contig.df$contamination[is.na(contig.df$contamination)] <- "valid_no_hit"
contig.df$contamination[contig.df$contamination=="nempep4"] <- "valid_nempep"

## this is rather annotation
## annotation with Bm orthologous sequence
BMblast <- read.blast.best("/home/ele/Data/454/annotation/blast/full_vs_bm.blt")
names(BMblast)[names(BMblast)%in%c("V1")] <- "contig"
names(BMblast)[names(BMblast)%in%c("V2")] <- "Bm.hit"
names(BMblast)[names(BMblast)%in%c("V12")] <- "Bm.bit"
contig.df <- merge(contig.df, BMblast[, c("contig", "Bm.hit", "Bm.bit")], all.x=TRUE)


## anntoation with Ce orthologous sequence
CEblast <- read.blast.best("/home/ele/Data/454/annotation/blast/full_vs_ce.blt")
names(CEblast)[names(CEblast)%in%c("V1")] <- "contig" 
names(CEblast)[names(CEblast)%in%c("V2")] <- "Ce.hit"
names(CEblast)[names(CEblast)%in%c("V12")] <- "Ce.bit"

contig.df <- merge(contig.df, CEblast[, c("contig", "Ce.hit", "Ce.bit")], all.x=TRUE)

## Patches ## patches ## patches

bm.rRNA.contigs <- contig.df[grepl("A8Q5G4|A8PJZ7|A8NLF8|A8PEQ0",
                                   contig.df$Bm.hit), "contig"]

contig.df[contig.df$contig%in%bm.rRNA.contigs, "contamination"] <- "AcrRNA"

## a valid Ce protein
contig.df[contig.df$contig%in%"Acrassus_rep_c255", "contamination"] <- "valid_nempep"

## nt ###################################################################
nt <- read.blast.best("/home/ele/Data/454/post_assembly_screening/fullest_vs_nt.blt")
gi.nt.blast <- sub("gi\\|(\\d+)\\|.*", "\\1", nt$V2)
nt <- cbind(nt, gi.nt.blast)

gis.nt <- unique(gi.nt.blast)
write(as.character(gis.nt), "/home/ele/Data/454/post_assembly_screening/tmp-nt.acc")
taxa.nt <- read.delim( pipe ("/home/ele/thesis/454/tax4gi_nt.pl /home/ele/Data/454/post_assembly_screening/tmp-nt.acc"), header=FALSE, sep=",", as.is=TRUE)
names(taxa.nt) <- c("gi.nt.blast", "taxid.nt", "family.nt", "phylum.nt", "kingdom.nt")
nt <- merge(nt, taxa.nt, all.x=TRUE)

contig.df <- merge(contig.df, nt[, c(2, 14:17)], by.x="contig", by.y="V1", all.x=TRUE)
contig.df[is.na(contig.df$taxid.nt) | contig.df$taxid.nt==0,
          (length(contig.df)-2):length(contig.df)] <- "No_hit"


## nr ##########################################################################
nr <- read.blast.best("/home/ele/Data/454/post_assembly_screening/fullest_vs_nr.blt")
gi.nr.blast <- sub("gi\\|(\\d+)\\|.*", "\\1", nr$V2)
nr <- cbind(nr, gi.nr.blast)

gis.nr <- unique(gi.nr.blast)
write(as.character(gis.nr), "/home/ele/Data/454/post_assembly_screening/tmp-nr.acc")

taxa.nr <- read.delim( pipe ("/home/ele/thesis/454/tax4gi_nr.pl /home/ele/Data/454/post_assembly_screening/tmp-nr.acc"), header=FALSE, sep=",", as.is=TRUE)
## 1 entry changed from my blast to taxonomy version
## -> no information available, moved to "no-hit category"
taxa.nr[taxa.nr$taxid.nr==0, c(3:5)] <- NA
names(taxa.nr) <- c("gi.nr.blast", "taxid.nr", "family.nr", "phylum.nr", "kingdom.nr")
nr <- merge(nr, taxa.nr, all.x=TRUE)

contig.df <- merge(contig.df, nr[, c(2, 14:17)], by.x="contig", by.y="V1", all.x=TRUE)
contig.df[is.na(contig.df$taxid.nr) | contig.df$taxid.nr==0,
          (length(contig.df)-2):length(contig.df)] <- "No_hit"

##################################################################################

phyl.plot.nt <- ggplot(subset(contig.df, contig.df$kingdom.nt!="No_hit"), aes(x=phylum.nt)) +
  facet_wrap(~ kingdom.nt, scales="free") +
  geom_bar() + 
  opts(axis.text.x=theme_text(angle=-90, hjust=0))

phyl.plot.nr <- ggplot(subset(contig.df, contig.df$kingdom.nr!="No_hit"), aes(x=phylum.nr)) +
  facet_wrap(~ kingdom.nr, scales="free") +
  geom_bar() + 
  opts(axis.text.x=theme_text(angle=-90, hjust=0))

gc.cont.phylum <- ggplot(subset(contig.df, contig.df$kingdom.nr%in%"Metazoa" &
                                (contig.df$phylum.nr%in%"Nematoda"|
                                 contig.df$phylum.nr%in%"Chordata")),
                   aes(x=gc, color=contamination)) + 
  geom_density() +
  facet_wrap( ~ phylum.nr)

png("../figures/phylum_plots.png", width=2000, height=1000, res=144)

# Set up the page
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)

# Make each plot, in the correct location
print(phyl.plot.nt, vp = vplayout(1, 1 ))
print(phyl.plot.nr, vp = vplayout(1, 2 ))
dev.off() 




###################################################
### chunk number 10: post.con
###################################################
#line 868 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
post.table <- rbind(table(contig.df$contamination),
                    rbind(tapply(contig.df$uniq_coverage, contig.df$contamination, mean)))

row.names(post.table) <- c("number",
                           "mean coverage")

nvalid <- nrow(contig.df[grepl("valid", contig.df$contamination) ,])

## Examples of selecting from the data
select_AC <- grepl("valid", contig.df$contamination) &
                            contig.df$kingdom.nt%in%c("Metazoa", "No_hit") &
                            contig.df$kingdom.nr%in%c("Metazoa", "No_hit") &
                            !contig.df$phylum.nr%in%c("Chordata") &
                            !contig.df$phylum.nt%in%c("Chordata") 


select_AC_MN <- select_AC & contig.df$category%in%"MN"

contig.df$Ac <- select_AC
contig.df$AcMN <- select_AC_MN

n_su_Ac <- nrow(contig.df[contig.df$Ac, ])
n_suMN_Ac <- nrow(contig.df[contig.df$AcMN, ])

contig.df$seq.origin <- contig.df$contamination
contig.df[!contig.df$phylum.nr%in%c("Nematoda", "No_hit"),  "seq.origin"] <- contig.df[!contig.df$phylum.nr%in%c("Nematoda", "No_hit"), "phylum.nr"]

other.king <- summary.factor(contig.df[!contig.df$Ac &
                                       !contig.df$seq.origin%in%c("Nematoda", "No_hit",
                                                                  "eelmRNA", "Chordata", "AcrRNA"),
                             "kingdom.nr"])

other.phyl <- summary.factor(contig.df[!contig.df$Ac &
                                       !contig.df$seq.origin%in%c("Nematoda", "No_hit",
                                                                  "eelmRNA", "Chordata", "AcrRNA"),
                             "phylum.nr"])



###################################################
### chunk number 11: conservation
###################################################
#line 936 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
nr.all <- read.delim("/home/ele/Data/454/post_assembly_screening/fullest_vs_nr.blt", header=FALSE)

gi.nr.all.blast <- sub("gi\\|(\\d+)\\|.*", "\\1", nr.all$V2)
nr.all <- cbind(nr.all, gi.nr.all.blast)

gis.nr.all <- unique(gi.nr.all.blast)
write(as.character(gis.nr.all), "/home/ele/Data/454/post_assembly_screening/tmp-nr.all.acc")

taxa.nr.all <- read.delim( pipe ("/home/ele/thesis/454/tax4gi_nr.pl /home/ele/Data/454/post_assembly_screening/tmp-nr.all.acc"), header=FALSE, sep=",", as.is=TRUE)
## 1 entry changed from my blast to taxonomy version
## -> no information available, moved to "no-hit category"

names(taxa.nr.all) <- c("gi.nr.all.blast", "taxid.nr.all", "family.nr.all", "phylum.nr.all", "kingdom.nr.all")
nr.all <- merge(nr.all, taxa.nr.all, all.x=TRUE)


names(nr.all)[2] <- "contig"

nr.all <- nr.all[order(nr.all$contig, nr.all$V12, decreasing=TRUE),]

is.conserved <- function(taxa.blast, level, value, bit.threshold) {
  B <- subset(taxa.blast, as.numeric(V12)>as.numeric(bit.threshold) )
  if(nrow(B) < 1){return("No_hit")}
  if(all(B[, level]%in%value)){return(FALSE)}
  if(any(!B[,level]%in%value)){return(TRUE)}
  else{stop("Error in consevation function: Unrecognized case")}
}

## Clade
spirurina <- c("Anisakidae", "Ascarididae", "Heterocheilidae", "Quimperiidae", "Raphidascarididae", "Toxocaridae", "Cosmocercidae", "Kathlaniidae", "Ascaridiidae", "Aspidoderidae", "Heterakidae", "Cucullanidae", "Acuariidae", "Camallanidae", "Diplotriaenidae", "Anguillicolidae", "Daniconematidae", "Dracunculidae", "Philometridae", "Skyrjabillanidae", "Homungellidae", "Mesidionematidae", "Onchocercidae", "Setariidae", "Ungellidae", "unclassified Filarioidea", "Gnathostomatidae", "Habronematidae", "Hedruridae", "Tetrameridae", "Physalopteridae", "Skrjabillanidae", "Cystidicolidae", "Gongylonematidae", "Rhabdochonidae", "Spirocercidae", "unclassified Spiruroidea", "Thelaziidae")

conserved.clade.50 <- as.data.frame(t(rbind(by(nr.all, nr.all$contig,
                                               is.conserved, "family.nr.all", spirurina,  50))))
names(conserved.clade.50) <- "conserved.clade.50"

conserved.clade.80 <- as.data.frame(t(rbind(by(nr.all, nr.all$contig,
                                               is.conserved, "family.nr.all", spirurina,  80))))
names(conserved.clade.80) <- "conserved.clade.80"

## Nematodes

conserved.phylum.50 <- as.data.frame(t(rbind(by(nr.all, nr.all$contig,
                                                is.conserved, "phylum.nr.all", "Nematoda",  50))))
names(conserved.phylum.50) <- "conserved.phylum.50"


conserved.phylum.80 <- as.data.frame(t(rbind(by(nr.all, nr.all$contig,
                                                is.conserved, "phylum.nr.all", "Nematoda",  80))))
names(conserved.phylum.80) <- "conserved.phylum.80"

## Animals

conserved.kingdom.50 <- as.data.frame(t(rbind(by(nr.all, nr.all$contig,
                                                 is.conserved, "kingdom.nr.all", "Metazoa",  50))))
names(conserved.kingdom.50) <- "conserved.kingdom.50"

conserved.kingdom.80 <- as.data.frame(t(rbind(by(nr.all, nr.all$contig,
                                                 is.conserved, "kingdom.nr.all", "Metazoa",  80))))
names(conserved.kingdom.80) <- "conserved.kingdom.80"


## merge it
contig.df <- merge(contig.df, as.data.frame(conserved.clade.50),
                   by.x="contig", by.y="row.names", all.x=TRUE)

contig.df <- merge(contig.df, as.data.frame(conserved.clade.80),
                   by.x="contig", by.y="row.names", all.x=TRUE)

contig.df <- merge(contig.df, as.data.frame(conserved.phylum.50),
                   by.x="contig", by.y="row.names", all.x=TRUE)

contig.df <- merge(contig.df, as.data.frame(conserved.phylum.80),
                   by.x="contig", by.y="row.names", all.x=TRUE)

contig.df <- merge(contig.df, as.data.frame(conserved.kingdom.50),
                   by.x="contig", by.y="row.names", all.x=TRUE)

contig.df <- merge(contig.df, as.data.frame(conserved.kingdom.80),
                   by.x="contig", by.y="row.names", all.x=TRUE)

## Compute novelity from conservation

contig.df$novel.50 <- ifelse(contig.df$conserved.clade.50%in%TRUE &
                             contig.df$conserved.phylum.50%in%TRUE &
                             contig.df$conserved.kingdom.50%in%TRUE,
                             "conserved", 
                             ifelse(contig.df$conserved.clade.50%in%TRUE &
                                    contig.df$conserved.phylum.50%in%TRUE &
                                    contig.df$conserved.kingdom.50%in%FALSE,
                                    "novel.in.metazoa", 
                                    ifelse(contig.df$conserved.clade.50%in%TRUE &
                                           contig.df$conserved.phylum.50%in%FALSE &
                                           contig.df$conserved.kingdom.50%in%FALSE,
                                           "novel.in.nematoda",
                                           ifelse(contig.df$conserved.clade.50%in%FALSE &
                                                  contig.df$conserved.phylum.50%in%FALSE &
                                                  contig.df$conserved.kingdom.50%in%FALSE,
                                                  "novel.in.clade3",
                                                  ifelse(contig.df$conserved.clade.50%in%c("No_hit", NA) &
                                                         contig.df$conserved.phylum.50%in%c("No_hit", NA) &
                                                         contig.df$conserved.kingdom.50%in%c("No_hit", NA) ,
                                                         "novel.in.Ac",
                                                         NA)))))

contig.df$novel.80 <- ifelse(contig.df$conserved.clade.80%in%TRUE &
                             contig.df$conserved.phylum.80%in%TRUE &
                             contig.df$conserved.kingdom.80%in%TRUE,
                             "conserved", 
                             ifelse(contig.df$conserved.clade.80%in%TRUE &
                                    contig.df$conserved.phylum.80%in%TRUE &
                                    contig.df$conserved.kingdom.80%in%FALSE,
                                    "novel.in.metazoa", 
                                    ifelse(contig.df$conserved.clade.80%in%TRUE &
                                           contig.df$conserved.phylum.80%in%FALSE &
                                           contig.df$conserved.kingdom.80%in%FALSE,
                                           "novel.in.nematoda",
                                           ifelse(contig.df$conserved.clade.80%in%FALSE &
                                                  contig.df$conserved.phylum.80%in%FALSE &
                                                  contig.df$conserved.kingdom.80%in%FALSE,
                                                  "novel.in.clade3",
                                                  ifelse(contig.df$conserved.clade.80%in%c("No_hit", NA) &
                                                         contig.df$conserved.phylum.80%in%c("No_hit", NA) &
                                                         contig.df$conserved.kingdom.80%in%c("No_hit", NA),
                                                         "novel.in.Ac",
                                                         NA)))))

## ## A venn diagramm for conservation and missingness
## clade.hit.50 <- contig.df[contig.df$conserved.clade.50%in%c(TRUE, FALSE),
##                      "contig"]
## clade.hit.50.l = match(clade.hit.50, contig.df$contig)[!is.na(match(clade.hit.50, contig.df$contig))]

## phylum.hit.50 <- contig.df[contig.df$conserved.phylum.50%in%c(TRUE, FALSE),
##                      "contig"]
## phylum.hit.50.l = match(phylum.hit.50, contig.df$contig)[!is.na(match(phylum.hit.50,
##   contig.df$contig))]

## kingdom.hit.50 <- contig.df[contig.df$conserved.kingdom.50%in%c(TRUE,FALSE),
##                      "contig"]
## kingdom.hit.50.l = match(kingdom.hit.50, contig.df$contig)[!is.na(match(kingdom.hit.50, contig.df$contig))]

## ## Bm.hit.50 <- contig.df[!is.na(contig.df$Bm.hit) &
## ##                        as.numeric(as.character(contig.df$Bm.bit))>50,
## ##                        "contig"]
## ## Bm.hit.50.l = match(Bm.hit.50, contig.df$contig)[!is.na(match(Bm.hit.50, contig.df$contig))]


## ## Ce.hit.50 <- contig.df[!is.na(contig.df$Ce.hit) &
## ##                        as.numeric(as.character(contig.df$Ce.bit))>50,
## ##                        "contig"]
## ## Ce.hit.50.l = match(Ce.hit.50, contig.df$contig)[!is.na(match(Ce.hit.50, contig.df$contig))]


## venn.diagram(list(Kingdom = kingdom.hit.50.l,
##                   Phylum = phylum.hit.50.l,
##                   Clade = clade.hit.50.l),
##              filename = "/home/ele/thesis/454/figures/conservation_venn.tiff")



###################################################
### chunk number 12: imputedP4E
###################################################
#line 1099 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

get.seq.and.orf <- function (filepath){
  IMPFAS <- read.sequences(filepath)
  contig <-  gsub("\\s.*", "", names(IMPFAS))
  mstrand <- grepl("minus", names(IMPFAS))
  method <- sapply(names(IMPFAS),
                             function (x){
                               unlist(strsplit(x, " "))[[length(unlist(strsplit(x, " ")))]]})
  
  starts <- as.numeric(as.character(lapply(gregexpr("[a-z]", IMPFAS), function (x) x[[1]])))
  st <- gregexpr("[a-z]", IMPFAS)
  sst <- as.numeric(as.character(lapply(st, function (x) {
    length(attr(x, "match.length"))})))
  end <- sst + starts
  
  names(IMPFAS) <- contig
  names(method) <- contig
  
  STST <- as.data.frame(cbind(contig,
                              mstrand,
                              method,
                              starts,
                              end,
                              imp = IMPFAS))
  return(STST)
}

imputed <- get.seq.and.orf("/home/ele/Data/454/prot4EST/fullest_assembly_imputed.fasta")
## to get only the MN subset
## imputed <- imputed[imputed$contig%in%MN,]

imputed$method <- factor(imputed$method, levels=unique(imputed$method))


imputed$starts <- as.numeric(as.character(imputed$starts))
imputed$starts[imputed$starts==-1] <- 0
imputed$end <- as.numeric(as.character(imputed$end))
imputed$end[imputed$end==-1] <- 0
imputed$imp <- as.character(imputed$imp)

imp.sum <- table(imputed$mstrand, imputed$method)

row.names(imp.sum) <- ifelse(row.names(imp.sum), "minus strand", "plus strand")
names(imp.sum) <- gsub("p4e..", "", names(imp.sum))


imputed$over.5p <-  ifelse(as.logical(imputed$mstrand),
                           imputed$end+3 < nchar(imputed$imp),
                           imputed$starts > 3)

imputed$startcod <- ifelse(as.logical(imputed$mstrand),
                           revcom(substr(imputed$imp, imputed$end-3, imputed$end-1)),
                           substr(imputed$imp, imputed$starts, imputed$starts+2))

imputed$stopcod <- ifelse(as.logical(imputed$mstrand),
                          revcom(substr(imputed$imp, imputed$starts-3, imputed$starts-1)),
                          ## lower-case letters have to be removed
                          ## where ther a not 3 bases overhang 
                          gsub("[a-z]", "",
                               substr(imputed$imp, imputed$end, imputed$end+2 )))

imputed$full.5p <-  ifelse(imputed$startcod=="atg" & imputed$over.5p, 
                           TRUE,
                           FALSE)

imputed$full.3p <-  imputed$stopcod%in%c("TAG", "TAA", "TGA")


imputed$full.length <- imputed$full.3p & imputed$full.5p

imputed[(imputed$end-imputed$starts)%%3>0, 6] <- toupper(imputed[(imputed$end-imputed$starts)%%3>0, 6])
imputed[(imputed$end-imputed$starts)%%3>0, c(2:5, 7:12)] <- NA

contig.df <- merge(contig.df, imputed, all.x=TRUE)

stopcods <- table(contig.df$stopcod[contig.df$stopcod%in%c("TAG", "TAA", "TGA")])

same.same <- summary.factor(tolower(as.character(contig.df[contig.df$Ac, "imp"]))==tolower(as.character(contig.df[contig.df$Ac, "seq"])))

imp.pile <- read.table(pipe("cut -f1,2,4 /home/ele/Data/454/mapping/all_vs_full_imputed_uq.pileup"))
names(imp.pile) <- c("contig", "base", "imputed.coverage")
per.con.imp <- data.frame(imputed.coverage= tapply(imp.pile$imputed.coverage, imp.pile$contig, mean, rm.na=T))
contig.df <- merge(contig.df, per.con.imp, by.x="contig", by.y="row.names", all.x=TRUE)

## finally also add the peptide sequence
Ac.pep <- read.sequences("/home/ele/Data/454/prot4EST/A.crassus_p4ePro.fsa")
names(Ac.pep) <- (gsub(" .*|_\\d.*|\t.*|\\s.*" , "", names(Ac.pep)))

contig.df <- merge(contig.df, as.data.frame(Ac.pep), by.x="contig", by.y="row.names", all.x=TRUE)



###################################################
### chunk number 13: snp
###################################################
#line 1193 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

read.var <- function (path) {
  VAR <- read.delim(path, header=FALSE, as.is=TRUE)
  names(VAR) <- c("contig", "base", "from", "to", "nfrom", "nto", "perc", "sfrom", "sto", "qfrom", "qto", "pval")
  VAR$perc <- as.numeric(gsub("%", "", VAR$perc))
  ## Not a iupac ambiguity in reference
  ## this was not the case anyways 
  VAR <- subset(VAR, VAR$from%in%c("A", "C", "G", "T"))
  ## Reference not wrong
  VAR$nfrom <- as.numeric(as.character(VAR$nfrom))
  VAR <- subset(VAR, VAR$nfrom>1)
  VAR <- VAR[1:12]
### remove snps looking like multiple allele
  combined <- paste(VAR$contig, VAR$base)
  VAR <- VAR[! combined %in% combined[duplicated(combined)],]
  return(VAR)
}



VAR <- read.var ("/home/ele/Data/454/mapping/all_vs_full_imputed_uq.varsnp")
## ## screen for fish and bacterial contamination
## VAR <- VAR[VAR$contig%in%nematode.contig.df$contig,]
## ## use only the good category
## VAR <- VAR[VAR$contig%in%MN,]

## reduce to get only the sure Ac 
VAR <- VAR[VAR$contig%in%contig.df[contig.df$Ac, "contig"], ]

VAR <- merge(VAR, imputed)

factors <- c("contig", "from", "to", "method")
numerics <- c("base", "nfrom", "nto", "perc", "sfrom", "sto", "qfrom", "qto", "pval", "starts", "end")

VAR[, factors] <- lapply(VAR[, factors], as.factor)
VAR[, numerics] <- lapply(VAR[, numerics], function(x) as.numeric(as.character(x)))
VAR[, "mstrand"] <- as.logical(VAR[, "mstrand"])


VAR$inORF <- as.factor(VAR$base>VAR$starts & VAR$base<VAR$end)

## Exclude all snps in contigs without prediction: 
VAR <- VAR[!is.na(VAR$inORF),]

frame <- sapply(1:nrow(VAR), function (i){
  if (VAR[i,"mstrand"]){
    (VAR[i,"end"]-VAR[i,"base"])%%3
  }
  else {
    (VAR[i,"base"]-VAR[i,"starts"])%%3
  }
})
  
VAR$inFRAME <- as.factor(frame+1)
VAR[!as.logical(VAR$inORF), "inFRAME"] <- NA

transversion.transition <- function (VARobj){
  get.trans <- function (x) {
    trans <- summary(x$from:x$to)
    trans <- trans[trans!=0]
  }
  transf <- get.trans(VARobj)
  get.vers <- function (x){
    transitions <- c("A:G", "G:A", "C:T", "T:C")
    sition <- sum(x[names(x)%in%transitions])
    version <- sum(x[!names(x)%in%transitions])
    res <- cbind(transitions=sition, transversions=version, ratioTS.TV=sition/version)
    return(as.data.frame(res))
  }
  get.vers(transf)
}

find.homopolymers <- function (VARobj, width=5, size=4){
  su <- vector(length=nrow(VARobj))
  for (b in 1:nrow(VARobj)){
    su[b] <- substr(VARobj[b, "imp"], VARobj[b, "base"]-width, VARobj[b, "base"]+width)
  }
  regex <- paste("\\.*T{", size, ",}|A{", size, ",}|G{", size, ",}|C{", size, ",}\\.*", sep="")
  homopol <- sapply(su, function (x) grepl(regex, x , ignore.case=TRUE))
  names(homopol) <- su
  return(homopol)
}

## Look through the ratios of transversion vs. transition to find an
## optimal exclusion space for width of the window and size of the
## homopolymer

TS.TV <- lapply (3:9,  function(width) {
  lapply (3:6,  function (size) {
    transversion.transition(VAR[!find.homopolymers(VAR, width=width, size=size), ])
  })
})


T <- melt(TS.TV, id=c("transitions", "transversions", "ratioTS.TV"))
T <- rename(T, c(L2="size", L1="width"))
T$width <- T$width+2
T$size <- T$size+2

## How the ratio of transv-transs changes of parameterspace in homopolymer exclusion

trans.vers.parameter.plot <- ggplot(T, aes(ratioTS.TV,
                                           transversions+transitions,
                                           shape=as.factor(size),
                                           color=as.factor((width*2)+1))) +
  geom_point(size=4) +
  geom_smooth(aes(group = size<4), se=FALSE, method="lm") +
  scale_x_continuous("ratio transistions/transversions") +
  scale_y_continuous("total numter of snps found") +
  scale_color_discrete("width of screening\nwindow around snp") +
  scale_shape_discrete("homopolymer length\nthreshold for exclusion") #+
##  opts(title="Effect of parameters for homopolymer\nexclusion on the number of snps found and transversion/transition ratio")

## ## BASE ONTOLOGY
## #
## # Code from Mark Blaxter, modified by John Davey,
## # Translated from Perl to R by Emanuel Heitlinger:

## # A phase 1 any change is nonsynonymous
## # B phase 2 any change is nonsynonymous
## # C phase 3 any change is nonsynonymous
## # D phase 1 change to CT is nonsynonymous
## # E phase 2 change to CT is nonsynonymous
## # F phase 3 change to CT is nonsynonymous
## # G phase 1 change to AG is nonsynonymous
## # H phase 2 change to AG is nonsynonymous
## # I phase 3 change to AG is nonsense
## # K phase 1 change to GT is nonsynonymous
## # L phase 2 change to A is nonsense, to anything else is nonsynonymous
## # J phase 3 change to G is nonsynonymous
## # M phase 3 change to G is nonsense, to A is nonsynonymous
## # N phase 3 any change synonymous
## # O phase 1 change to T nonsense, others nonsynonymous
## # P phase 3 change to AG is nonsynonymous
## # Q phase 1 change to T nonsense, to G nonsynonymous
## # R phase 2 change to AG nonsense, others nonsynonymous
## # S phase 3 change to A nonsense, others nonsynonymous
## # T phase 3 change to A nonsense, G nonsynonymous

## # W all changes are unknown # EH added 08/23/2011

## #        a           g           c           t
## #
## # a     aaa K OBF   aga R QBF   aca T ABN   ata I ABJ
## #       aag K OBF   agg R KBF   acg T ABN   atg M ABC
## #       aac N ABP   agc S ABP   acc T ABN   atc I ABJ
## #       aat N ABP   agt S ABP   act T ABN   att I ABJ
## #
## # g     gaa E OBF   gga G OBN   gca A ABN   gta V ABN
## #       gag E OBF   ggg G ABN   gcg A ABN   gtg V ABN
## #       gac D ABP   ggc G ABN   gcc A ABN   gtc V ABN
## #       gat D ABP   ggt G ABN   gct A ABN   gtt V ABN
## #
## # c     caa Q OBF   cga R QBN   cca P ABN   cta L GBN
## #       cag Q OBF   cgg R KBN   ccg P ABN   ctg L GBN
## #       cac H ABP   cgc R ABN   ccc P ABN   ctc L ABN
## #       cat H ABP   cgt R ABN   cct P ABN   ctt L ABN
## #
## # t     taa * AEF   tga * AEC   tca S ARN   tta L GRF
## #       tag * ABF   tgg W ALS   tcg S ALN   ttg L GLF
## #       tac Y ABI   tgc C ABT   tcc S ABN   ttc F ABP
## #       tat Y ABI   tgt C ABT   tct S ABN   ttt F ABP


base.ontology.encode <- function(x){  
  ## # wrap in a control function for iupac and other "bad"
  ## # bases
  if (nchar(gsub("[bdefhijklmnopqrsuvwxyz]", "", x, ignore.case=TRUE)) != 3){
    return (paste(rep("W", times = nchar(x)), collapse=""))
  }
  else {
    ## # Set up Base Ontology vector
    base.ontology.encode.string = c(
      "aaa" = "OBF",
      "aag" = "OBF",
      "aac" = "ABP",
      "aat" = "ABP",
      "aga" = "QBF",
      "agg" = "KBF",
      "agc" = "ABP",
      "agt" = "ABP",
      "aca" = "ABN",
      "acg" = "ABN",
      "acc" = "ABN",
      "act" = "ABN",
      "ata" = "ABJ",
      "atg" = "ABC",
      "atc" = "ABJ",
      "att" = "ABJ",
      ##
      "gaa" = "OBF",
      "gag" = "OBF",
      "gac" = "ABP",
      "gat" = "ABP",
      "gga" = "OBN",
      "ggg" = "ABN",
      "ggc" = "ABN",
      "ggt" = "ABN",
      "gca" = "ABN",
      "gcg" = "ABN",
      "gcc" = "ABN",
      "gct" = "ABN",
      "gta" = "ABN",
      "gtg" = "ABN",
      "gtc" = "ABN",
      "gtt" = "ABN",
      ##
      "caa" = "OBF",
      "cag" = "OBF",
      "cac" = "ABP",
      "cat" = "ABP",
      "cga" = "QBN",
      "cgg" = "KBN",
      "cgc" = "ABN",
      "cgt" = "ABN",
      "cca" = "ABN",
      "ccg" = "ABN",
      "ccc" = "ABN",
      "cct" = "ABN",
      "cta" = "GBN",
      "ctg" = "GBN",
      "ctc" = "ABN",
      "ctt" = "ABN",
    ##
      "taa" = "AEF",
      "tag" = "ABF",
      "tac" = "ABI",
      "tat" = "ABI",
      "tga" = "AEC",
      "tgg" = "ALS",
      "tgc" = "ABT",
      "tgt" = "ABT",
      "tca" = "ARN",
      "tcg" = "ALN",
      "tcc" = "ABN",
      "tct" = "ABN",
      "tta" = "GRF",
      "ttg" = "GLF",
      "ttc" = "ABP",
      "ttt" = "ABP"
      );
    return(base.ontology.encode.string[x])
  }
}
  
base.ontology.decode = list(
  "A" = c("A" = "Nonsynonymous",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "B" = c("A" = "Nonsynonymous",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "C" = c( "A" = "Nonsynonymous",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "D" = c( "A" = "Synonymous", "C" = "Nonsynonymous",
    "G" = "Synonymous", "T" = "Nonsynonymous" ),
  "E" = c( "A" = "Synonymous", "C" = "Nonsynonymous",
    "G" = "Synonymous", "T" = "Nonsynonymous" ),
  "F" = c( "A" = "Synonymous", "C" = "Nonsynonymous",
    "G" = "Synonymous", "T" = "Nonsynonymous" ),
  "G" = c( "A" = "Nonsynonymous",  "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Synonymous" ),
  "H" = c( "A" = "Nonsynonymous",  "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Synonymous" ),
  "I" = c( "A" = "Nonsense",  "C" = "Synonymous",
    "G" = "Nonsense",  "T" = "Synonymous" ),
  "J" = c( "A" = "Synonymous", "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Synonymous" ),
  "K" = c( "A" = "Synonymous", "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "L" = c( "A" = "Nonsense",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "M" = c( "A" = "Nonsynonymous",  "C" = "Synonymous",
    "G" = "Nonsense",  "T" = "Synonymous" ),
  "N" = c( "A" = "Synonymous", "C" = "Synonymous",
    "G" = "Synonymous", "T" = "Synonymous" ),
  "O" = c( "A" = "Nonsynonymous",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsense" ),
  "P" = c( "A" = "Nonsynonymous",  "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Synonymous" ),
  "Q" = c( "A" = "Synonymous", "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsense" ),
  "R" = c( "A" = "Nonsense",  "C" = "Nonsynonymous",
    "G" = "Nonsense",  "T" = "Nonsynonymous" ),
  "S" = c( "A" = "Nonsense",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "T" = c( "A" = "Nonsense",  "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Synonymous" ),
  "X" = c( "A" = "Nonsense",  "C" = "Nonsense",
    "G" = "Nonsense",  "T" = "Nonsense" ),
  "W" = c("A" = NA,  "C" = NA,
    "G" = NA, "T" = NA , "R" = NA,
    "Y" = NA, "S" = NA, "W" = NA,
    "K" = NA, "M" = NA, "B" = NA,
    "D" = NA, "H" = NA, "V" = NA,
    "N" = NA, "X" = NA)
  );

contig.df$coding <- ifelse(!is.na(contig.df$mstrand) &
                           as.logical(contig.df$mstrand),  
                           ## revcom to get the sequence on the plus strand always
                           ## sometimes, there are Ns after the end of orfs!!
                           revcom (gsub("[ACGTRYSWKMBDHVNX]", "", contig.df$imp)),
                           gsub("[ACGTRYSWKMBDHVNX]", "", contig.df$imp)
                           )

## get the coding region coverd >8x in snp calling
imp.cov.pile <- imp.pile[imp.pile$imputed.coverage>7, ]
b <- by(imp.cov.pile, imp.cov.pile$contig,
        function (x) c(x[1,"base"], x[nrow(x), "base"]))
cov.8.borders <- as.data.frame(do.call("rbind", b))
names(cov.8.borders) <- c("cov.8.start", "cov.8.stop")

cov.df <- merge(contig.df, cov.8.borders, by.x="contig", by.y="row.names")

seq.cov.8 <- by(cov.df, cov.df$contig, function(x){
  substring(x$imp, x$cov.8.start, x$cov.8.stop)})
Sc <- as.data.frame(cbind(seq.cov.8))

cov.df <- merge(cov.df, Sc, by.x="contig", by.y="row.names")

cov.df$coding.cov8 <- ifelse(!is.na(cov.df$mstrand) &
                                as.logical(cov.df$mstrand),  
                                ## revcom to get the sequence on the plus strand always
                                ## sometimes, there are Ns after the end of orfs!!
                                revcom (gsub("[ACGTRYSWKMBDHVNX]", "", cov.df$seq.cov.8)),
                                gsub("[ACGTRYSWKMBDHVNX]", "", cov.df$seq.cov.8)
                                )



add.ontology.sites <- function (contig.df.obj, cod.col.name){
  coding <- contig.df.obj[, cod.col.name]
  codons <- lapply(coding, function (x){
    substring(x,seq(1,nchar(x), by=3), seq(3,nchar(x), by=3))})
  ontology <- lapply(codons, function (x) { sapply (x, base.ontology.encode)})
  ontology <- as.character(sapply(ontology, paste, collapse=""))
   ## reverse back to be able to use normal snp cooridingates
  ontology <- ifelse(as.logical(contig.df.obj$mstrand),
                     strReverse(ontology),
                     ontology)
   s.sites <- sapply(1:length(ontology), function (i) {
     split.ont <- unlist(strsplit(ontology[[i]], ""))
     split.cod <- unlist(strsplit(coding[[i]], ""))
     decoded <- lapply(split.ont, function (a){
       base.ontology.decode[[a]]})
     reduced.decoded <- lapply(1:length(decoded), function (x) {
       subset(decoded[[x]], names(decoded[[x]])!=toupper(split.cod[[x]]))})
     s <- lapply(reduced.decoded, function (w) {
       length(w[grepl("Syn", w)])/length(w[!is.na(w)])})
     n <- lapply(reduced.decoded, function (w) {
       length(w[grepl("Non", w)])/length(w[!is.na(w)])})
     nsyn.sites <- sum(unlist(n), na.rm=TRUE)
     syn.sites <- sum(unlist(s), na.rm=TRUE)
     cbind(nsyn.sites, syn.sites)
   })
  s.sites <- t(s.sites)
  c.df <- data.frame(cbind(ontology, s.sites))
  names(c.df)[c(ncol(c.df)-1, ncol(c.df))] <- c("nsyn.sites", "syn.sites")
  c.df$nsyn.sites <- as.numeric(as.character(c.df$nsyn.sites))
  c.df$syn.sites <- as.numeric(as.character(c.df$syn.sites))
  row.names(c.df) <- contig.df.obj$contig
  return(c.df)
}

coding.df <- add.ontology.sites(contig.df[nchar(contig.df$coding)>2, ], "coding")
coding.cov.df <- add.ontology.sites(cov.df[nchar(cov.df$coding.cov8)>2, ], "coding.cov8")
names(coding.cov.df)[2:3] <-  c("cov8.nsyn.sites", "cov8.syn.sites")

contig.df <- merge(contig.df, coding.df, by.x="contig", by.y="row.names", all.x=TRUE )
contig.df <- merge(contig.df, coding.cov.df[,2:3], by.x="contig", by.y="row.names", all.x=TRUE )

VAR <- merge(VAR, contig.df[,c("contig", "ontology", "nsyn.sites", "syn.sites", "cov8.nsyn.sites", "cov8.syn.sites")] )

add.effect.by.ontology <- function (VARobj){
  orfbase <- (VARobj$base-VARobj$starts+1)
  ontology <- as.character(VAR$ontology)
  onto.codes <- sapply(1:nrow(VARobj), function (i) {
    ifelse(as.logical(VARobj[i,]$inORF) &
           orfbase[i]<nchar(ontology[i]) ,
           unlist(strsplit(ontology[i], ""))[[orfbase[[i]]]], NA)})

  effect <- sapply(1:length(onto.codes), function (i) {
    base.ontology.decode[[onto.codes[[i]]]][[VARobj[i,]$to]]})
  effect[unlist(lapply(effect, is.null))] <- "outside ORF"
  effect <- unlist(effect)
  VARobj <- cbind(VARobj, effect)
  return(VARobj)
}

VAR <- add.effect.by.ontology(VAR)

get.dn.ds <- function(VARobj){
  (nrow(VARobj[grepl("Non*", VARobj$effect),])/
   sum(VARobj[!duplicated(VARobj$contig), "cov8.nsyn.sites"], na.rm=T))/
  (nrow(VARobj[VARobj$effect=="Synonymous",])/
   sum(VARobj[!duplicated(VARobj$contig), "cov8.syn.sites"], na.rm=T))}
  
## Look through the overall dn/ds to find an
## optimal exclusion space for width of the window and size of the
## homopolymer

O.DN.DS <- sapply (3:9,  function(width) {
  sapply (3:6,  function (size) {
    get.dn.ds(VAR[!find.homopolymers(VAR, width=width, size=size), ])
  })
})


D <- melt(O.DN.DS)
names(D)[1:2] <- c("size", "width")
D$width <- D$width+2
D$size <- D$size+2
T <- merge(D, T)
T <- rename(T, replace=c(value="dn.ds"))

dn.ds.parameter.plot <- ggplot(T, aes(dn.ds,
                                      transversions+transitions,
                                      shape=as.factor(size),
                                      color=as.factor(width*2+1))) +
                               geom_point(size=4) +
                               geom_smooth(aes(group = size<4), se=FALSE, method="lm") +
                               scale_x_reverse("ratio non-synonymous/synonymous snps") +
                               scale_y_continuous("total numter of snps found") +
                               scale_color_discrete("width of screening\nwindow around snp") +
                               scale_shape_discrete("homopolymer length\nthreshold for exclusion")
                               ## opts(title="Effect of parameters for
                               ## homopolymer\nexclusion on the number
                               ## of snps found and\nratio of non-synonymous/synonymous snps")

#### Decided on the SNP sequence quality-screening
VARq <- (VAR[!find.homopolymers(VAR),])

#### Now consider mapping quality
VARq.plot <- VARq
VARq.plot$perc[VARq.plot$perc>50] <- 51

position.perc.plot <- ggplot(VARq.plot[!is.na(VARq.plot$inFRAME),], aes(x = perc, fill = inFRAME)) +
  stat_bin(binwidth=1, color="black", breaks=seq(0, 52, by=2),
           position=position_dodge(width=(1.4))) +
  scale_x_continuous("percent of minority allele", 
                     breaks= seq(0, 52, by=2),
                     labels = c(seq(0, 48, by=2), ">50", "")) +
  scale_y_continuous("number of snps found") +
  scale_fill_discrete("position in frame")

effect.perc.plot <- ggplot(VARq.plot[!is.na(VARq.plot$effect),], aes(x = perc, fill = effect)) +
  stat_bin(binwidth=1, color="black", breaks=seq(0, 52, by=2),
           position=position_dodge(width=c(1.7))) +
  scale_x_continuous("percent of minority allele",
                     breaks= seq(0, 52, by=2),
                     labels = c(seq(0, 48, by=2), ">50", "")) +
  scale_y_continuous("number of snps found") +
  scale_fill_discrete("effect of snp")

VARq.plot$fromto <- VARq.plot$nfrom + VARq.plot$nto
VARq.plot$fromto[VARq.plot$fromto>100] <- 101

position.cov.plot <- ggplot(VARq.plot[!is.na(VARq.plot$inFRAME),], aes(x = fromto, fill = inFRAME)) +
  stat_bin(binwidth=1, color="black", breaks=seq(0, 105, by=5),
           position=position_dodge(width=(2.8))) +
  scale_x_continuous("coverage at snp positon",
                     breaks = seq(0, 105, by=5),
                     labels = c(seq(0, 95, by=5), ">100", "")) +
  scale_y_continuous("number of snps found") +
  scale_fill_discrete("position in frame")

effect.cov.plot <- ggplot(VARq.plot[!is.na(VARq.plot$effect),], aes(x = fromto, fill = effect)) +
  stat_bin(binwidth=1, color="black", breaks=seq(0, 105, by=5),
           position=position_dodge(width=(2.8))) +
  scale_x_continuous("coverage at snp positon",
                     breaks = seq(0, 105, by=5),
                     labels = c(seq(0, 95, by=5), ">100", "")) +
  scale_y_continuous("number of snps found") +
  scale_fill_discrete("effect of snp")

#### discard snps with <7% minority allel
VARqp <- VARq[VARq$perc>7, ]

## Without screening the percentage

TSVs <- do.call("rbind",
                     by(VARq,
                        VARq$contig,
                        function (x) transversion.transition(x)))

contig.df.wo <- merge(contig.df, TSVs, by.x="contig", by.y="row.names", all.x=TRUE)

N.by.S <- by(VARq, VARq$contig, get.dn.ds)
contig.df.wo <- merge(contig.df.wo, t(as.data.frame.list(N.by.S)), by.x="contig", by.y="row.names", all.x=TRUE)
names(contig.df.wo)[ncol(contig.df.wo)] <- "dn.ds"
contig.df.wo$dn.ds[is.infinite(contig.df.wo$dn.ds)] <- NA

sy.per <- as.data.frame(rbind(table(VARq$contig, VARq$effect)))

contig.df.wo <- merge(contig.df.wo, sy.per, by.x="contig", by.y="row.names", all.x=TRUE)

dn.ds.coverage.plot.wo <- ggplot(contig.df.wo, aes(x = dn.ds, y= imputed.coverage)) +
  geom_point() +
  scale_y_log10() +
  stat_smooth(method="lm")

dn.ds.total.snps.plot.wo <- ggplot(contig.df.wo, aes(y=transversions+transitions, x=dn.ds)) +
  geom_point() +
  scale_y_continuous("total number of snps in contig") +
  stat_smooth(method="lm")

## With screening the percentage

TSVs <- do.call("rbind",
                     by(VARqp,
                        VARqp$contig,
                        function (x) transversion.transition(x)))

contig.df <- merge(contig.df, TSVs, by.x="contig", by.y="row.names", all.x=TRUE)


N.by.S <- by(VARqp, VARqp$contig, get.dn.ds)
contig.df <- merge(contig.df, t(as.data.frame.list(N.by.S)), by.x="contig", by.y="row.names", all.x=TRUE)
names(contig.df)[ncol(contig.df)] <- "dn.ds"
contig.df$dn.ds[is.infinite(contig.df$dn.ds)] <- NA

sy.per <- as.data.frame(rbind(table(VARqp$contig, VARqp$effect)))

contig.df <- merge(contig.df, sy.per, by.x="contig", by.y="row.names", all.x=TRUE)

dn.ds.coverage.plot <- ggplot(contig.df, aes(x = dn.ds, y= imputed.coverage)) +
  geom_point() +
  scale_y_log10() +
  stat_smooth(method="lm")

dn.ds.total.snps.plot <- ggplot(contig.df, aes(y=transversions+transitions, x=dn.ds)) +
  geom_point() +
  scale_y_continuous("total number of snps in contig") +
  stat_smooth(method="lm")

### The contigs looking like TS.TV~30 are the infinite values!

dens.dn.ds <- ggplot(contig.df[!is.na(contig.df$dn.ds),],
                     aes(dn.ds, ..count..)) + geom_bar(binwidth=0.1)


## no need to  beta calc from http://www.biomedcentral.com/1471-2164/11/310
## and the original idea http://www.biomedcentral.com/1471-2164/9/312



###################################################
### chunk number 14: 
###################################################
#line 1739 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

png("../figures/snp_ex_parameter_plots.png", width=2000, height=1000, res=144)

# Set up the page
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)

# Make each plot, in the correct location
print(trans.vers.parameter.plot, vp = vplayout(1, 1 ))
print(dn.ds.parameter.plot, vp = vplayout(1, 2 ))

dev.off() 


png("../figures/snp_pos_eff_plots.png", width=5000, height=3000, res=288)

# Set up the page
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)

# Make each plot, in the correct location
print(position.perc.plot, vp = vplayout(1, 1 ))
print(position.cov.plot, vp =vplayout(1, 2 ))
print(effect.perc.plot, vp = vplayout(2, 1 ))
print(effect.cov.plot, vp = vplayout(2, 2 ))
dev.off() 

ggsave("../figures/dens_dn_ds.png", dens.dn.ds)


png("../figures/dn_ds_scales.png", width=2000, height=1000, res=144)

# Set up the page
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)

# Make each plot, in the correct location
print(dn.ds.total.snps.plot.wo, vp = vplayout(1, 1 ))
print(dn.ds.coverage.plot.wo, vp = vplayout(1, 2 ))

print(dn.ds.total.snps.plot, vp = vplayout(2, 1 ))
print(dn.ds.coverage.plot, vp = vplayout(2, 2 ))

dev.off() 


###################################################
### chunk number 15: snps.again
###################################################
#line 1793 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"


nCov8bases <- nrow(imp.pile[imp.pile$contig%in%contig.df[contig.df$Ac, "contig"] &
                            imp.pile$imputed.coverage>7, ])

nrawOrf <- tapply(VAR$inORF, VAR$inORF, length)
nrawFrame <- tapply(VAR$inFRAME, VAR$inFRAME, length)

raw.ts.tv <- transversion.transition(VAR)
tsv.raw.orf <- round(do.call(rbind, by(VAR, VAR$inORF, transversion.transition)), 2)
raw.dn.ds <- round(get.dn.ds(VAR), 2)


s.base <- sum(contig.df[contig.df$Ac, "cov8.syn.sites"], na.rm=TRUE)
n.base <- sum(contig.df[contig.df$Ac, "cov8.nsyn.sites"], na.rm=TRUE)

sSNP <- sum(contig.df[contig.df$Ac, "Synonymous"], na.rm=TRUE)
nSNP <- sum(contig.df[contig.df$Ac, "Nonsynonymous"], na.rm=TRUE)

per.base <- round(nrow(VARqp)/(nCov8bases/1000), 2)

s.per.s.base <- round(sSNP/(s.base/1000), 2)
n.per.n.base <- round(nSNP/(n.base/1000), 2)



###################################################
### chunk number 16: snp.by.population
###################################################
#line 1819 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

## from within /home/ele/Data/454/mapping/mapping_each_lib:  
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
  
  V <- cbind(VCFQ[ , c(1:4)], gt.df)
  names(V)[1:4] <- c( "contig", "base", "from", "to")

  ## merging with the general VARscan estimate to remove bad 
  VAR.pop <- merge(VARqp, V, by=c( "contig", "base", "from", "to"))
    
  relative.h <- function (x){
    length(x[x%in%c("0/1")])/length(x[x%in%c("0/0", "1/1")])
  }

  if (sum) {return (apply(VAR.pop[, c("10F", "179F", "M175", "KS4F", "UW07F")],
                          2, relative.h))}
  else {return(VAR.pop)}
}


## consistently over heterozygoszy-calling-thresholds european samples
## are less heterozygot

## super.geno(VCF, 0, 30)
## super.geno(VCF, 0, 40)
## super.geno(VCF, 0, 50)
## super.geno(VCF, 0, 60)

## super.geno(VCF, 10, 30)
## super.geno(VCF, 10, 40)
## super.geno(VCF, 10, 50)
## super.geno(VCF, 10, 60)

## super.geno(VCF, 30, 30)
## super.geno(VCF, 30, 40)
## super.geno(VCF, 30, 50)
  
het.table <- as.data.frame(super.geno(VCF, 10, 10, TRUE))
VAR.pop <- super.geno(VCF, 10, 10, FALSE)

sep.gt <- function (gt.col){
  l <- strsplit(as.character(gt.col), "/")
  s <- lapply(l, function (x) {if (is.na(x[1])){ rep(x,2)}
              else{x}})
  unlist(s)}

Rh <- VAR.pop[, c("10F", "179F", "KS4F", "M175", "UW07F")]
### remove loci with only one allele

## like this internal relationship does not work
## Rh <- Rh[ apply(Rh, 1, function (x) ("0/1"%in%x ) | ("0/0"%in%x & "1/1"%in%x) ), ]

## like this internal relationship does  work
Rh <- Rh[ apply(Rh, 1, function (x) ("0/1"%in%x & "0/0"%in%x  ) |
                ("0/0"%in%x & "1/1"%in%x)|
                ("0/1"%in%x & "1/1"%in%x)), ]

number.hetero.info.snps <- nrow(Rh)

Rh.data <- as.data.frame(rbind(sep.gt(Rh$"10F"),
                               sep.gt(Rh$"179F"),
                               sep.gt(Rh$"M175"),
                               sep.gt(Rh$"KS4F"),
                               sep.gt(Rh$"UW07F")))

## shows that it worked:
ir.test <- hh(Rh.data, 100, "ir")
hl.test <- hh(Rh.data, 100, "hl")
hh.test <- hh(Rh.data, 100, "sh")

## mean(ir.test)
## quantile(ir.test, probs=c(0.025, 0.975))

## mean(hl.test)
## quantile(hl.test, probs=c(0.025, 0.975))

## mean(hh.test)
## quantile(hh.test, probs=c(0.025, 0.975))

## Internal relatedness, a multilocus heterozygosity measure developed
## by Amos et al. (2001): Amos W, Worthington Wilmer J, Fullard K et
## al (2001) The influence of parental relatedness on reproductive
## success. Proc R Soc Lond B 268:2021-2027
het.table <- cbind(het.table, ir(Rh.data))

## homozygosity by loci: Aparicio et al. (2007)?(2006) : Aparicio JM,
## Ortego J, Cordero PJ (2006) What should we weigh to estimate
## heterozygosity, alleles or loci? Mol Ecol 15:4659-4665
het.table <- cbind(het.table, hl(Rh.data))

## standardized heterozygosity: Coltman et al. (1999): Coltman DW,
## Pilkington JG, Smith JA et al (1999) Parasite-mediated selection
## against inbred Soay sheep in a free-living island
## population. Evolution 53:1259-1267
het.table <- cbind(het.table, sh(Rh.data))
names(het.table) <- c("rel.het", "int.rel", "ho.loci", "std.het")



###################################################
### chunk number 17: annotation
###################################################
#line 1988 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
GO.annot <- read.delim("/home/ele/Data/454/annotation/annot8r/output_Ac/GO.csv", sep=",", header=FALSE)
names(GO.annot) <- c("pept_id", "go_term", "pcf", "descr", "slim", "besthit", "bestscore", "bestev", "hitnum", "maxhits", "fraction") 

## reset the shorted names used in annotation to the proper contig-names
GO.annot$pept_id <- gsub("Ac_", "Acrassus_", GO.annot$pept_id)

GO.annot$go_term <- as.character(gsub(" ", "", GO.annot$go_term))
GO.annot$pcf <- gsub(" C", "CC", GO.annot$pcf)
GO.annot$pcf <- gsub(" F", "MF", GO.annot$pcf)
GO.annot$pcf <- gsub(" P", "BP", GO.annot$pcf)

EC.annot <- read.delim("/home/ele/Data/454/annotation/annot8r/output_Ac/EC.csv", sep=",", header=FALSE)
names(EC.annot) <- c("pept_id", "ec_id", "descr", "besthit", "bestscore", "bestev", "hitnum", "maxhits", "fraction") 
EC.annot$pept_id <- gsub("Ac_", "Acrassus_", EC.annot$pept_id)

KEGG.annot <- read.delim("/home/ele/Data/454/annotation/annot8r/output_Ac/KEGG.csv", sep=",", header=FALSE)
names(KEGG.annot) <- c("pept_id", "ko_id", "path", "descr", "besthit", "bestscore", "bestev", "hitnum", "maxhits", "fraction") 
KEGG.annot$pept_id <- gsub("Ac_", "Acrassus_", KEGG.annot$pept_id)

IPR.annot <- read.delim("/home/ele/Data/454/annotation/iprscan/MN.pep.fasta.iprscan", header=FALSE)
names(IPR.annot)[1] <- "pept_id"

## annotation with Bm orthologous sequence
all.bm.hit <- as.character(contig.df[contig.df$Ac, "Bm.hit"])
all.bm.hit <- all.bm.hit[!is.na(all.bm.hit)]
all.bm.hit <- gsub("UniRef100_", "", all.bm.hit)

uniProt <- useMart("unimart", dataset="uniprot")
all.bm.names <- getBM(attributes =c("name", "protein_name"),
                          filter="accession",values=all.bm.hit, mart=uniProt)

all.bm.names$name <- gsub("_BRUMA", "", all.bm.names$name)
all.bm.names$Bm.hit <- paste("UniRef100_", all.bm.names$name, sep="")
names(all.bm.names)[names(all.bm.names)=="protein_name"] <- "Bm.annot"
contig.df <- merge(all.bm.names[, 2:3], contig.df, by="Bm.hit", all.y=TRUE, sort=FALSE)

## anntoation with Ce orthologous sequence
all.ce.hit <- as.character(contig.df[contig.df$Ac, "Ce.hit"])
all.ce.hit <- all.ce.hit[!is.na(all.ce.hit)]

WORMB <- useMart("WS220" , dataset="wormbase_gene")
all.ce.rnai <- getBM(attributes = c("transcript", "brief_ident",
                        "rnai_phenotype_phenotype_label"),
                          filter="transcript",values=all.ce.hit, mart=WORMB)

rnai <- as.data.frame(c(by(all.ce.rnai, all.ce.rnai$transcript,
                           function (x) any(grepl("lethal", x$rnai_phenotype_phenotype_label)))))
names(rnai) <- "Ce.rnai"

ce.names <- all.ce.rnai[, c("transcript","brief_ident")]
names(ce.names) <- c("Ce.hit", "Ce.annot")
ce.names <- ce.names[!duplicated(ce.names$Ce.hit), ]
ce.names <- merge(ce.names, rnai, by.x="Ce.hit", by.y="row.names", all.x=TRUE)

contig.df <- merge(contig.df, ce.names, by="Ce.hit", all.x=TRUE, sort=FALSE)


## anntoation with nempep descriptions
nempep.nuc <- read.sequences("/home/ele/db/blastdb/nembase4/nembase4_nuc.fasta")
nempep <- data.frame(Cluster.ID=substring(names(nempep.nuc), 1, 8))
nempep$nempep.annot <- gsub("\\w{3}\\d{5}_1 Cluster: ", "",  names(nempep.nuc))
nempep <- nempep[!duplicated(nempep$Cluster.ID),]

BL.nempep <- read.blast.best("/home/ele/Data/454/post_assembly_screening/fullest_vs_nembase_pro.blt")
BL.nempep$Pept.ID <- substring(BL.nempep$V2, 1, 8)
## notice difference between cluster identifiers and peptide identifiers
BL.nempep$Cluster.ID <- gsub("(\\w)(\\w)P", "\\1\\2C", BL.nempep$Pept.ID)
names(BL.nempep)[names(BL.nempep)=="V12"] <- "nempep.bit"
names(BL.nempep)[names(BL.nempep)=="V1"] <- "contig"

nemp <- merge(BL.nempep[, c(1, 12, 14)], nempep, by="Cluster.ID")
names(nemp)[names(nemp)=="Cluster.ID"] <- "nempep.hit"
contig.df <- merge(contig.df, nemp, by="contig", all.x=TRUE)

n.nempep <- length(contig.df[contig.df$Ac &
                             !grepl("no annotation", contig.df$nempep.annot) &
                             !is.na(contig.df$nempep.annot), "nempep.annot"])

## annotagion with uniprot-name of nr-orthologous sequence
gis.nr <- unique(as.character(nr$gi.nr.blast))

## write(gis.nr, "/home/ele/Data/454/annotation/all_gis.acc")
gi2uniprot <- read.delim("/home/ele/Data/454/annotation/gi2uniprot.tab")
names(gi2uniprot) <- c("gi", "uniprot")
gi2uniprot$accession <- gsub("_.*", "", gi2uniprot$uniprot)

gi.uni <- getBM(attributes = c("accession", "protein_name"),
                filter="accession",values=gi2uniprot$accession, mart=uniProt)

gi.uni.acc <- merge(gi2uniprot, gi.uni)

nr <- merge(nr, gi.uni.acc[,c(2, 4)], by.x="gi.nr.blast", by.y="gi", all.x=TRUE)
nr <- nr[!duplicated(nr$V1),]

names(nr)[names(nr)%in%c("V1")] <- "contig"
names(nr)[names(nr)%in%c("V12")] <-  "nr.bit"
names(nr)[names(nr)%in%c("protein_name")] <- "nr.uniprot.annot"

contig.df <- merge(contig.df,
                   nr[, c("gi.nr.blast", "contig", "nr.bit", "nr.uniprot.annot")],
                   all.x=TRUE)

n.uniprot <- length(contig.df[contig.df$Ac & !is.na(contig.df$nr.uniprot.annot), "nr.uniprot.annot"])



###################################################
### chunk number 18: ann.venn
###################################################
#line 2095 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
GOannotations <- unique(GO.annot$pept_id[GO.annot$pept_id%in%contig.df[contig.df$Ac,"contig"]])
nGO <- length(GOannotations)

ECannotations <- unique(EC.annot$pept_id[EC.annot$pept_id%in%contig.df[contig.df$Ac,"contig"]])
nEC <- length(ECannotations)

KEGGannotations <- unique(KEGG.annot$pept_id[KEGG.annot$pept_id%in%contig.df[contig.df$Ac,"contig"]])
nKEGG <- length(KEGGannotations)


Allannotations <- unique(c(GOannotations, ECannotations, KEGGannotations))
nallannotations <- length(Allannotations)

annotationVenn <- venn.diagram(list(GO   = match(GOannotations, Allannotations),
                                       EC   = match(ECannotations, Allannotations),
                                       KEGG = match(KEGGannotations, Allannotations)),
                                  filename = NULL)

GOannotations.MN <- unique(GO.annot$pept_id[GO.annot$pept_id%in%contig.df[contig.df$AcMN,"contig"]])
nGO.MN <- length(GOannotations.MN)

ECannotations.MN <- unique(EC.annot$pept_id[EC.annot$pept_id%in%contig.df[contig.df$AcMN,"contig"]])
nEC.MN <- length(ECannotations.MN)

KEGGannotations.MN <- unique(KEGG.annot$pept_id[KEGG.annot$pept_id%in%contig.df[contig.df$AcMN,"contig"]])
nKEGG.MN <- length(KEGGannotations.MN)

IPRannotations <- unique(as.character(IPR.annot[!IPR.annot$V4%in%"Seg",
                                                "pept_id"]))
nIPR <- length(IPRannotations)

Allannotations.MN <- unique(c(GOannotations.MN, ECannotations.MN, KEGGannotations.MN, IPRannotations))

annotationVennMN <- venn.diagram(list(GO   = match(GOannotations.MN, Allannotations.MN),
                  EC   = match(ECannotations.MN, Allannotations.MN),
                  KEGG = match(KEGGannotations.MN, Allannotations.MN),
                  IPR = match(IPRannotations, Allannotations.MN)),
                                 filename = NULL)

not.annot <- nrow(contig.df[contig.df$Ac & !contig.df$contig%in%Allannotations,])
not.annot.MN <- nrow(contig.df[contig.df$AcMN & !contig.df$contig%in%Allannotations.MN,])

png("/home/ele/thesis/454/figures/annotataionVenn.png", width=400, height=600)
vp1 <- viewport(x=0.01, y=0.5, w=0.98, h=0.51,
                just=c("left", "bottom"))
vp2 <- viewport(x=0.01, y=0, w=0.98, h=0.49,
                just=c("left", "bottom"))
push.viewport(vp1)
grid.roundrect()
grid.text("a",x=unit(0,"npc"),y=unit(0.97,"npc"), gp=gpar(fontsize=18))
grid.text(not.annot, x=unit(0.07,"npc"), y=unit(0.1,"npc"), gp=gpar(fontfamily="serif"))
grid.draw(annotationVenn)
pop.viewport()
push.viewport(vp2)
grid.roundrect()
grid.text("b",x=unit(0,"npc"),y=unit(0.97,"npc"), gp=gpar(fontsize=18))
grid.text(not.annot.MN, x=unit(0.07,"npc"), y=unit(0.1,"npc"), gp=gpar(fontfamily="serif"))
grid.draw(annotationVennMN)
pop.viewport()
dev.off()            



###################################################
### chunk number 19: sigp
###################################################
#line 2158 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

signalip <- as.data.frame(read.delim("/home/ele/Data/454/annotation/signalp/Ac_p4ePro.signalip4",
                                     skip=2, sep="", header=FALSE,
                                     strip.white=TRUE, comment.char=";",
                                     as.is=TRUE))

signalip.names <- c("contig", "Cmax", "pos", "Ymax",
                    "pos", "Smax", "pos", "Smean", "D",
                    "sigp.pred", "Dmaxcut", "Networks-used")

names(signalip) <- signalip.names

## translate the contig-names back
signalip$contig <- gsub("Ac_", "Acrassus_", signalip$contig)

sum.sigP <- function (sigp.obj) {
  ifelse(sigp.obj$sigp.pred=="Y",
         paste("Yes", gsub("SignalP",
                           "", sigp.obj$"Networks-used"),
               sep=""), "No")}

signalip$sigP <- sum.sigP(signalip)

contig.df <- merge(contig.df, signalip[, c("contig", "sigP")], all.x=TRUE)

sigP <- contig.df[contig.df$sigP%in%"Yes-noTM", "contig"]
sigPtm <- contig.df[contig.df$sigP%in%"Yes-TM", "contig"]


signalip.bm <- as.data.frame(read.delim("/home/ele/Data/454/annotation/signalp/uniref100_bm.signalp4",
                                     skip=2, sep="", header=FALSE,
                                     strip.white=TRUE, comment.char=";",
                                     as.is=TRUE))
names(signalip.bm) <- signalip.names
signalip.bm$sigP <- sum.sigP(signalip.bm)

signalip.ce <- as.data.frame(read.delim("/home/ele/Data/454/annotation/signalp/wormpep220.signalp4",
                                        skip=2, sep="", header=FALSE,
                                        strip.white=TRUE, comment.char=";",
                                        as.is=TRUE))
names(signalip.ce) <- signalip.names
signalip.ce$sigP <- sum.sigP(signalip.ce)



###################################################
### chunk number 20: signn.venn
###################################################
#line 2205 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

nsigP <- length(contig.df[contig.df$Ac &  contig.df$sigP%in%"Yes-noTM", "sigP"])
nsigPtm <- length(contig.df[contig.df$Ac &  contig.df$sigP%in%"Yes-TM", "sigP"])


nsigP.Bm <- nrow(signalip.bm[signalip.bm$sigP%in%"Yes-noTM",])
nsigPtm.Bm <- nrow(signalip.bm[signalip.bm$sigP%in%"Yes-TM",])

nsigP.Ce <- nrow(signalip.bm[signalip.ce$sigP%in%"Yes-noTM",])
nsigPtm.Ce <- nrow(signalip.bm[signalip.ce$sigP%in%"Yes-TM",])




###################################################
### chunk number 21: annot.compare
###################################################
#line 2219 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

slim.annot8r <- function (path, species){ 
  l <- list()
  l <- lapply(c("F", "C", "P"), function (x){
    d <- read.delim(paste(path, "/piedata_", x, sep=""), skip=1)
    d$species <- species
    d
  })
  l[[1]]$ontology <- "Molecular Function"
  l[[2]]$ontology <- "Cellular compartment"
  l[[3]]$ontology <- "Biological process"
  melt(l)
}

Ac.slim <- slim.annot8r("/home/ele/Data/454/annotation/annot8r/output_Ac/", "Anguillicola crassus")
Bm.slim <- slim.annot8r("/home/ele/Data/454/annotation/annot8r/output_Bm/", "Brugia malayi")
Ce.slim <- slim.annot8r("/home/ele/Data/454/annotation/annot8r/output_Ce/", "Caenorhabditis elegans")

slim <- rbind(Ac.slim, Bm.slim, Ce.slim)

GO.bm.com <- ggplot(slim, aes(x=description, fill=species, weight=value)) +
  geom_bar(position="dodge") +
  scale_x_discrete(breaks=slim$description,
                   label=gsub("^\n", "",
                     gsub("and\nnucleic\nacid\nmetabolic\nprocess",
                          "and nucleic acid\nmetabolic process",
                          gsub(" ", "\n", slim$description)))) +
  facet_wrap(~ontology, ncol=1, scales="free") +
    opts(legend.text = theme_text(face="italic"))

ggsave("/home/ele/thesis/454/figures/go_bm_com.png", GO.bm.com, width=16, height=16)

Go.rep.Ac.Bm <- cor(Ac.slim$value, Bm.slim$value, method="spearman")
Go.rep.Ac.Ce <- cor(Ac.slim$value, Ce.slim$value, method="spearman")
Go.rep.Bm.Ce <- cor(Bm.slim$value, Ce.slim$value, method="spearman")



###################################################
### chunk number 22: overrep
###################################################
#line 2257 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

## subsets for contig.df and GO-annotation containing only dn.ds contigs
dn.ds.df <- subset(contig.df, contig.df$Ac & !is.na(dn.ds))
GO.dn.ds <- subset(GO.annot, pept_id%in%dn.ds.df$contig)

## the category of >.5 dn/ds contigs
big05 <- dn.ds.df[dn.ds.df$dn.ds>0.5 &
                  dn.ds.df$contig%in%GO.dn.ds$pept_id, "contig"]

## ## From the GOstats vignette
test.over.under.GOstats <- function (annotation, set, ontology="MF") {
  goframeData <- as.data.frame(cbind(frame.go_id=
                                     as.character(gsub(" ", "",
                                                       annotation[annotation$pcf%in%ontology,
                                                               "go_term"])),
                                     frame.Evidence="IEA",
                                     frame.gene_id=annotation[annotation$pcf%in%ontology,
                                       "pept_id"]))
  goFrame <- GOFrame(goframeData, organism="Anguillicola crassus")
  goAllFrame <- GOAllFrame(goFrame)
  gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
  over.params <- GSEAGOHyperGParams(name = "GSEA based annot arameters
                                             for over-representation",
                                    geneSetCollection = gsc,
                                    geneIds = as.character(set),
                                    universeGeneIds = unique(annotation[annotation$pcf%in%ontology,
                                      "pept_id"]),
                                    ontology = ontology,
                                    pvalueCutoff = 0.05,
                                    conditional = FALSE,
                                    testDirection = "over") 
  OVER <- hyperGTest(over.params)
  Over <- summary(OVER)
  if (nrow(Over)>0) Over <- cbind(Over, direction="Over") 
  under.params <- GSEAGOHyperGParams(name = "GSEA based annot arameters
                                             for under-representation",
                                     geneSetCollection = gsc,
                                     geneIds = as.character(set),
                                     universeGeneIds = unique(annotation[annotation$pcf%in%ontology,
                                       "pept_id"]),
                                     ontology = ontology,
                                     pvalueCutoff = 0.05,
                                     conditional = FALSE,
                                     testDirection = "under")
  UNDER <- hyperGTest(under.params)
  Under <- summary(UNDER)
  if (nrow(Under)>0) Under <- cbind(Under, direction="Under") 
  rbind(Over, Under)
}

## only for sure AC this time
mf.dn.ds <- test.over.under.GOstats(GO.dn.ds, big05, "MF")
bp.dn.ds <- test.over.under.GOstats(GO.dn.ds, big05, "BP")
cc.dn.ds <- test.over.under.GOstats(GO.dn.ds, big05, "CC")

## make a list of offspring terms to allow finding of 
GOMFOFF.list <- as.list(GOMFOFFSPRING)
GOMFOFF.list <- GOMFOFF.list[!is.na(GOMFOFF.list)]

## which are the amino acid transporter contigs (2) being in dn.ds>0.5
amino.acid.transporter.go <- mf.dn.ds[grepl("acid transmembrane transporter",
                                               mf.dn.ds$Term), 1]
amino.acid.tr.off <- unlist(GOMFOFF.list[amino.acid.transporter.go])

## get the term-contigs plus their boffspring term-contigs
## consistent with mf.dn.ds results (size)
amino.acid.tr.contigs <- unique(GO.dn.ds[GO.dn.ds$go_term%in%
                                            c(amino.acid.transporter.go, amino.acid.tr.off),
                                            "pept_id"])

peptidase.go <- mf.dn.ds[grepl("peptidase", mf.dn.ds$Term), 1]
peptidase.off <- unlist(GOMFOFF.list[peptidase.go])

## consistent with mf.dn.ds results
peptidase.contigs <- unique(GO.dn.ds[GO.dn.ds$go_term%in%
                                        c(peptidase.go, peptidase.off),
                                            "pept_id"])

## consistent with mf.dn.ds results (count)
peptidase.contigs.pos <- contig.df[contig.df$Ac & contig.df$contig%in%peptidase.contigs &
                                   contig.df$dn.ds>0.5 , "contig"]


## just to test that results for quality-contigs only are consistent.
mf.dn.ds.MN <- test.over.under.GOstats(subset(GO.dn.ds,
                                                 GO.dn.ds$pept_id%in%
                                                 contig.df[contig.df$AcMN,
                                                           "contig"]),
                                          big05, "MF")

bp.dn.ds.MN <- test.over.under.GOstats(subset(GO.dn.ds,
                                                 GO.dn.ds$pept_id%in%
                                                 contig.df[contig.df$AcMN,
                                                           "contig"]),
                                          big05, "BP")

cc.dn.ds.MN <- test.over.under.GOstats(subset(GO.dn.ds,
                                                 GO.dn.ds$pept_id%in%
                                                 contig.df[contig.df$AcMN,
                                                           "contig"]),
                                          big05, "CC")

## The same for KEGG doesn't work as I have not-uniqe annotations to a single contig
## keggframeData = data.frame(frame.path_id=as.character(gsub("^ K", "", KEGG.annot$ko_id)),
##   frame.gene_id=as.character(KEGG.annot$pept_id))

## ## HS
## keggframeData = data.frame(frame$path_id, frame$gene_id)
## keggFrame = KEGGFrame(keggframeData, organism = "Homo sapiens")

## keggFrame = KEGGFrame(keggframeData, organism = "Anguillicola crassus")

## gsc <- GeneSetCollection(keggFrame, setType = KEGGCollection())

## big05 <- as.character(keggframeData$frame.gene_id[1:500])
## universe <- as.character(keggframeData$frame.gene_id)


## kparams <- GSEAKEGGHyperGParams(name = "My Custom GSEA based annot Params",
##                                 geneSetCollection = gsc,
##                                 geneIds = big05, universeGeneIds = universe,
##                                 pvalueCutoff = 0.05, testDirection = "over")


## kOver <- hyperGTest(kparams)
## head(summary(kOver))
### WTF!!

dn.ds.df$sigP <- as.factor(dn.ds.df$sigP)
## this for the easy nn-data
u.test <- wilcox.test(dn.ds ~ sigP, data=dn.ds.df, conf.int = TRUE)


## sort the conservation according to breath 
contig.df$novel.50 <- factor(contig.df$novel.50,
                             levels=c("conserved", "novel.in.metazoa",
                               "novel.in.nematoda", "novel.in.clade3",
                               "novel.in.Ac"))

contig.df$novel.80 <- factor(contig.df$novel.80,
                             levels=c("conserved", "novel.in.metazoa",
                               "novel.in.nematoda", "novel.in.clade3",
                               "novel.in.Ac"))



###################################################
### chunk number 23: nn.dn.ds
###################################################
#line 2404 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

sigP.box <- ggplot(dn.ds.df, aes(dn.ds.df$sigP, dn.ds.df$dn.ds)) + 
  geom_boxplot() + 
  scale_y_log10("dn/ds") +
  scale_x_discrete("", 
                   breaks=levels(dn.ds.df$sigP),
                   labels=c("No signal peptide",
                     "Signal peptide")) +
  opts(title = "SignalP  prediction",
       axis.title.x = theme_text(size = 15), 
       axis.text.x = theme_text(size = 15)) +
  theme_bw()

ggsave("../figures/sigp_dn_ds.png", sigP.box)


###################################################
### chunk number 24: novel.dn.ds
###################################################
#line 2421 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

novel.dn.ds.50 <- ggplot(subset(contig.df, contig.df$Ac & !is.na(novel.50) & !is.na(dn.ds)),
                         aes(novel.50, dn.ds)) +
  geom_boxplot() +
  scale_y_log10("dn/ds") +
  scale_x_discrete("") +
  opts(title="evolutionary conseravation at bitsore threshold of 80") +
  theme_bw()

novel.dn.ds.80 <- ggplot(subset(contig.df, contig.df$Ac & !is.na(novel.80) & !is.na(dn.ds)),
                                aes(novel.80, dn.ds)) +
  geom_boxplot() +
  scale_y_log10("dn/ds") +
  scale_x_discrete("") +
  opts(title="evolutionary conseravation at bitsore threshold of 80") +
  theme_bw()

# Set up the page

png("../figures/conservation_dn_ds.png", width=2000, height=1000, res=144)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)

# Make each plot, in the correct location
print(novel.dn.ds.50, vp = vplayout(1, 1 ))
print(novel.dn.ds.80, vp = vplayout(1, 2 ))

grid.text("a", x=unit(0.01,"npc"), y=unit(0.99,"npc")) 
grid.text("b", x=unit(0.51,"npc"), y=unit(0.99,"npc")) 

dev.off() 



###################################################
### chunk number 25: novel.dn.ds.test
###################################################
#line 2458 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

### Nemenyi-Damico-Wolfe-Dunn test (joint ranking)
  ### Hollander & Wolfe (1999), page 244
  ### (where Steel-Dwass results are given)
get.ndwd <- function(data.df, test.var, fac.var) {
  NDWD <- oneway_test(test.var ~ fac.var, data = data.df,
                      ytrafo = function(data) trafo(data, numeric_trafo = rank),
                      xtrafo = function(data) trafo(data, factor_trafo = function(x)
                        model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
                      teststat = "max", distribution = approximate(B = 90000))
  return(NDWD)
}

dn.ds.df$novel.50 <- as.factor(dn.ds.df$novel.50)
dn.ds.df$novel.80 <- as.factor(dn.ds.df$novel.80)

novel.50.ndwd <- get.ndwd(dn.ds.df, dn.ds.df$dn.ds, dn.ds.df$novel.50)
novel.50.kw.ph <- pvalue(novel.50.ndwd, method = "single-step")

novel.80.ndwd <- get.ndwd(dn.ds.df, dn.ds.df$dn.ds, dn.ds.df$novel.80)
novel.80.kw.ph <- pvalue(novel.80.ndwd, method = "single-step")



###################################################
### chunk number 26: novel.nn
###################################################
#line 2483 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

sigP.50 <- as.data.frame.matrix(table(contig.df[contig.df$Ac, "novel.50"],
                                      grepl("Yes*", contig.df[contig.df$Ac, "sigP"])))
names(sigP.50) <- c("No.Sigp", "Sigp")
sigP.50 <- transform(sigP.50, proportion=(Sigp/(Sigp+No.Sigp))*100)

sigP.80 <- as.data.frame.matrix(table(contig.df[contig.df$Ac, "novel.80"],
                                      grepl("Yes*", contig.df[contig.df$Ac, "sigP"])))
names(sigP.80) <- c("No.Sigp", "Sigp")
sigP.80 <- transform(sigP.80, proportion=(Sigp/(Sigp+No.Sigp))*100)

sigP.con.test.50 <- fisher.test(contig.df[contig.df$Ac, "novel.50"]%in%
                                c("novel.in.Ac", "novel.in.nematoda"),
                                grepl("Yes*", contig.df[contig.df$Ac, "sigP"]))
sigP.con.test.80 <- fisher.test(contig.df[contig.df$Ac, "novel.80"]%in%
                                c("novel.in.Ac", "novel.in.nematoda"),
                                grepl("Yes*", contig.df[contig.df$Ac, "sigP"]))


sigP.50.p <- ggplot(subset(contig.df, contig.df$Ac & !is.na(novel.50) & !is.na(sigP)),
              aes(novel.50, fill=sigP)) +
       geom_bar(position="fill") + coord_polar(theta="y") +
  scale_y_continuous("") +
  scale_x_discrete("") +
  opts(title="proportion of TUGs in SignalP category\nby evolutionary conseravation categories at bitsore threshold of 50", plot.background = theme_blank()) +
  theme_bw()

sigP.80.p <- ggplot(subset(contig.df, contig.df$Ac & !is.na(novel.80) & !is.na(sigP)),
              aes(novel.80, fill=sigP)) +
  geom_bar(position="fill") + coord_polar(theta="y") +
  scale_y_continuous("") +
  scale_x_discrete("") +
  opts(title="proportion of TUGs in SignalP category\nby evolutionary conseravation categories at bitsore threshold of 50", plot.background = theme_blank()) + 
  theme_bw()

sigP.50.p.MN <- ggplot(subset(contig.df, contig.df$AcMN & !is.na(novel.50) & !is.na(sigP)),
              aes(novel.50, fill=sigP)) +
       geom_bar(position="fill") + coord_polar(theta="y") +
  scale_y_continuous("") +
  scale_x_discrete("") +
  opts(title="proportion of highCA contigs in SignalP category\nby evolutionary conseravation categories at bitsore threshold of 50", plot.background = theme_blank()) +
  theme_bw()

sigP.80.p.MN <- ggplot(subset(contig.df, contig.df$AcMN & !is.na(novel.80) & !is.na(sigP)),
              aes(novel.80, fill=sigP)) +
  geom_bar(position="fill") + coord_polar(theta="y") +
  scale_y_continuous("") +
  scale_x_discrete("") +
  opts(title="proportion of highCA contigs in SignalP category\nby evolutionary conseravation categories at bitsore threshold of 80", plot.background = theme_blank()) + 
  theme_bw()

png("../figures/signal_novel.png", width=2000, height=2000, res=144)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)


vp1 <- viewport(width = 0.54, height = 0.54, x = 0.25, y = 0.24)
vp2 <- viewport(width = 0.54, height = 0.54, x = 0.25, y = 0.72)
vp3 <- viewport(width = 0.54, height = 0.54, x = 0.74, y = 0.24)
vp4 <- viewport(width = 0.54, height = 0.54, x = 0.74, y = 0.72) 

print(sigP.50.p, vp = vp1)
print(sigP.80.p, vp = vp2)
print(sigP.50.p.MN, vp = vp3)
print(sigP.80.p.MN, vp = vp4)

grid.text("a", x=unit(0.01,"npc"), y=unit(0.99,"npc")) 
grid.text("b", x=unit(0.51,"npc"), y=unit(0.99,"npc")) 
grid.text("c", x=unit(0.01,"npc"), y=unit(0.52,"npc")) 
grid.text("d", x=unit(0.51,"npc"), y=unit(0.52,"npc")) 

dev.off() 

### tests for lethal Ce rnai
rnai.contigs <- table(contig.df[contig.df$Ac, "Ce.rnai"])

u.rnai <- wilcox.test(dn.ds ~ Ce.rnai, data=dn.ds.df, conf.int = TRUE)



###################################################
### chunk number 27: gene.expr
###################################################
#line 2745 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

## read the 454 mapping first  
## alternative method using library(GenomicRanges) would be:
## aligns <- readBamGappedAlignments("/home/ele/Data/RNAseq/mapping/AJ_T19M_1.bam")
## ie. bam <- scanBam("/home/ele/Data/RNAseq/mapping/AJ_T19M_1.bam", param=param)

## Rsamtools
what <- c("rname")# , "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(what = what)
  
files <- list.files("/home/ele/Data/454/mapping/mapping_each_lib/", "*sorted.bam$")

counts.for.files <- function (files){
  countsList <- list()
  for (fi in files) {
    bam <- scanBam(paste("/home/ele/Data/454/mapping/mapping_each_lib/", fi, sep=""), param=param)
    counts <- table(bam[[1]]$rname)
    counts.frame <- as.data.frame(counts)
    countsList[[fi]] <- counts.frame
  }
  return(countsList)
}

countsList <- counts.for.files(files)

counts.frame.long <- melt(countsList, id.vars="Var1")
counts.frame.long$variable <- NULL

counts.frame.wide <- reshape(counts.frame.long, v.names="value",
                             idvar="Var1", direction="wide", timevar="L1")
names(counts.frame.wide) <- gsub("^value.final\\.", "", names(counts.frame.wide))
names(counts.frame.wide) <- gsub("\\.sff.trimmed.sff.fasta.sorted.bam$", "", names(counts.frame.wide))
row.names(counts.frame.wide) <- counts.frame.wide$Var1
counts.frame.wide$Var1 <- NULL
n454.map <- sum(rowSums(counts.frame.wide))

## read the solexa tag-counts
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(what = what)
bam <- scanBam("/home/ele/Data/nlaIII-tags/UW08F_tags.bam", param=param)

nTags <- length(bam[[1]]$seq)
counts <- as.data.frame(table(bam[[1]]$rname))
names(counts) <- c("contig", "solexa.tags")

nMapped <- sum(counts$solexa.tags)
nUniq <- length(unique(bam[[1]]$seq))

nUMapped <- nTags-nMapped

lst <- lapply(names(bam[[1]]), function(elt) {
  do.call(c, unname(lapply(bam, "[[", elt)))
})
names(lst) <- names(bam[[1]])
df <- do.call("DataFrame", lst)

nUniqUMapped <- length(unique(df[is.na(df$rname), "seq"]))

counts <- merge(counts, counts.frame.wide, by.x="contig", by.y="row.names")

names(counts) <- gsub("179F", "T1", names(counts))
names(counts) <- gsub("10F", "T2", names(counts))
names(counts) <- gsub("KS4F", "E1", names(counts))
names(counts) <- gsub("UW07F", "E2", names(counts))
names(counts) <- gsub("L2R3", "L2", names(counts))
names(counts) <- gsub("M175", "M", names(counts))

counts$all.reads <- apply(counts[,2:length(counts)], 1, sum)

### merge all counts into the main data
contig.df <- merge(contig.df, counts, all.x=TRUE)

## discard all non Ac matches
cst <- contig.df[contig.df$Ac, c("E1", "E2", "L2", "M", "T1", "T2")]
rownames(cst) <- contig.df[contig.df$Ac, "contig"]

cst <- cst[rowSums(cst)>48, ]
cst <- cst[!is.na(cst$E1),]

conds.eel <- factor(c("EU", "EU", "L2", "M", "TW", "TW"))
conds.mf <- factor(c("F", "F", "L2", "M", "F", "F"))
conds.ad <- factor(c("Ad", "Ad", "L2", "Ad", "Ad", "Ad"))

exp.get.3.diff <- function (countsFrame) {
  cds.eel <- newCountDataSet(countsFrame, conds.eel)
  cds.eel <- estimateSizeFactors(cds.eel)
  cds.eel <- estimateVarianceFunctions(cds.eel)
  res.eel <- nbinomTest(cds.eel, "TW", "EU")

  cds.mf <- newCountDataSet(countsFrame, conds.mf)
  cds.mf <- estimateSizeFactors(cds.mf)
  cds.mf <- estimateVarianceFunctions(cds.mf)
  res.mf <- nbinomTest(cds.mf, "M", "F")

  cds.ad <- newCountDataSet(countsFrame, conds.ad)
  cds.ad <- estimateSizeFactors(cds.ad)
  cds.ad <- estimateVarianceFunctions(cds.ad)
  res.ad <- nbinomTest(cds.ad, "Ad", "L2")

  return(list(res.eel, res.mf, res.ad))
}

## for the contig-counts
contig.diff <- exp.get.3.diff(cst)

## collapsed for Bm orthologs 
bst <- contig.df[contig.df$Ac, c("Bm.hit", "E1", "E2", "L2", "M", "T1", "T2")]
bst <- bst[!is.na(bst$Bm.hit),]
bst <- do.call("rbind", by(bst, bst$Bm.hit, function (x) colSums(x[,2:7])))
bst <- bst[rowSums(bst)>48, ]
bst <- bst[!is.na(bst[,"E1"]),]

bm.diff <- exp.get.3.diff(bst)

## Collapsed for Ce orthologs
wst <- contig.df[contig.df$Ac, c("Ce.hit", "E1", "E2", "L2", "M", "T1", "T2")]
wst <- wst[!is.na(wst$Ce.hit),]
wst <- do.call("rbind", by(wst, wst$Ce.hit, function (x) colSums(x[,2:7])))
wst <- wst[rowSums(bst)>48, ]
wst <- wst[!is.na(wst[,"E1"]),]

ce.diff <- exp.get.3.diff(wst)

mf <- merge(contig.df[,c("contig", "Bm.hit", "Ce.hit")],
             contig.diff[[2]][,c("id", "padj")],
             by.x="contig", by.y="id")

mf <- merge(mf, bm.diff[[2]][,c("id", "padj")],
            by.x="Bm.hit", by.y="id", 
            all.x=TRUE)

mf <- merge(mf, ce.diff[[2]][,c("id", "padj")],
            by.x="Ce.hit", by.y="id", 
            all.x=TRUE)


eel <- merge(contig.df[,c("contig", "Bm.hit", "Ce.hit", "Bm.annot", "Ce.annot")],
             contig.diff[[1]][,c("id", "padj")],
             by.x="contig", by.y="id")

eel <- merge(eel, bm.diff[[1]][,c("id", "padj")],
            by.x="Bm.hit", by.y="id", 
            all.x=TRUE)

eel <- merge(eel, ce.diff[[1]][,c("id", "padj")],
            by.x="Ce.hit", by.y="id", 
            all.x=TRUE)



## resSig.ad <- subset(res.ad, padj < .1)
## sig.ad.contigs <- resSig.ad$id
## sig.ad.counts <- cst[sig.ad.contigs, ]


## vsd <- getVarianceStabilizedData( cds.mf )
## dists <- dist( t( vsd ) )
## idists <- as.matrix(dists)
##heatmap (idists , symm=TRUE, margins = c (7,7))



###################################################
### chunk number 28: save
###################################################
#line 2971 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
## write the MN.peptide sequence for iprscan prediction
MN.pep.fasta <- as.character(contig.df[contig.df$AcMN, "Ac.pep"])
names(MN.pep.fasta) <- as.character(contig.df[contig.df$AcMN, "contig"])
write.sequence(MN.pep.fasta, "/home/ele/thesis/454/MN.pep.fasta")

## for a more readable version of contig.df
readable <- !names(contig.df)%in%c("seq", "imp", "ontology", "coding", "Ac.pep")
write.csv(contig.df, "../A_crassus_contigs_full.csv")
write.csv(contig.df[,readable], "../A_crassus_contigs_readable.csv")

## for summary statistics a non-Ac reduced dataset
c.df <- contig.df[contig.df$Ac, ]

save.image("/home/ele/thesis/454/paper/paper.Rdata")


###################################################
### chunk number 29: 
###################################################
#line 3285 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
## print(xtable(TRIM, display=rep("d", times=7)), floating=FALSE)
trim_mark <- transform(TRIM, lowqal = short+lowq+dust+shortq)[1:6, ]
trim_mark <- trim_mark[,c(1,7)]

trim_mark <- merge(pre.screen.tab,
                   trim_mark, by="row.names")

trim_mark$life.st <- c("adult f", "adult f", "L2 lavae",
                       "adult m", "adult f", "adult f")

trim_mark$source.p <- c("Europe R", "Europe P", "Europe R",
                       "Asia C", "Asia C", "Asia W")

trim_mark[7,] <- c("Eel", "", "", "", "", length(eel.clean), sum(nchar(eel.clean)),
                   length(eel.all), length(eel.all)-length(eel.clean), "liver", "Taiwan")

trim_mark <- trim_mark[, c(1, 10, 11,  8, 9, 2, 4, 5, 3, 6, 7 )]

names(trim_mark)[c(1, 4)] <- c("library", "raw.reads")

write.csv(trim_mark, "/home/ele/thesis/454/paper/table1.csv")

trim_mark.tab <- xtable(trim_mark)
print(trim_mark.tab, floating=FALSE, include.rownames=FALSE)


###################################################
### chunk number 30: 
###################################################
#line 3314 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
T2 <- as.data.frame.array(tapply(contig.df$contig, contig.df$category%in%"MN", length))
rownames(T2) <- c("lowCA", "highCA")
names(T2) <- "total.contigs"
T2$rRNA.contigs <- table(contig.df$category%in%"MN", contig.df$seq.origin%in%"AcrRNA")[, 2]
T2$fish.contigs <- table(contig.df$category%in%"MN",
                         contig.df$seq.origin%in%c("eelmRNA", "Chordata"))[, 2]
T2$xeno.contigs <- table(contig.df$category%in%"MN",
                         !contig.df$seq.origin%in%c("Nematoda", "No_hit",
                                                  "eelmRNA", "Chordata", "AcrRNA"))[, 1]
T2$remaining.contigs <- table(contig.df$category%in%"MN", contig.df$Ac)[, 2]
T2$remaining.span <- tapply(contig.df$seq,
                            contig.df$category%in%"MN",
                            function (x) sum(nchar(as.character(x))))
T2$non.u.cov <- cov.mn
T2$cov <- cov.mn.u

T2$p4e <- table(c.df$category%in%"MN", gsub("p4e->", "", c.df$method))
T2$full.3p <- table(c.df$category%in%"MN", c.df$full.3p)[,2]
T2$full.5p <- table(c.df$category%in%"MN", c.df$full.5p)[,2]
T2$full.l <- table(c.df$category%in%"MN", c.df$full.length)[,2]

T2$GO <-  table(c.df$category%in%"MN", c.df$contig%in%GO.annot$pept_id)[, 2]
T2$EC <-  table(c.df$category%in%"MN", c.df$contig%in%EC.annot$pept_id)[, 2]
T2$KEGG <-  table(c.df$category%in%"MN", c.df$contig%in%KEGG.annot$pept_id)[, 2]
T2$IPR <-  table(c.df$category%in%"MN", c.df$contig%in%IPR.annot$pept_id)[, 2]

T2$nem.blast <- table(c.df$category%in%"MN", c.df$phylum.nr%in%"Nematoda")[, 2]

T2$any.blast <- table(c.df$category%in%"MN", c.df$phylum.nr%in%"No_hit")[,1]

T2 <- as.data.frame(t(T2))
T2 <- transform(T2, combined=lowCA+highCA)

T2["non.u.cov", "combined"] <- cov.nu.total
T2["cov", "combined"] <- cov.u.total

dig <- t(matrix(c(rep(0,28), rep(c(0, 3 , 3, 3), 2), rep(0,48)), nrow=4))
T2.tab <- xtable(T2, digits=dig)

print(T2.tab)


###################################################
### chunk number 31: 
###################################################
#line 3378 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
nov.tab <- rbind(summary.factor(contig.df[contig.df$Ac, "novel.50"]), 
                 summary.factor(contig.df[contig.df$Ac, "novel.80"]))
nov.tab.MN <- rbind(summary.factor(contig.df[contig.df$AcMN, "novel.50"]), 
                 summary.factor(contig.df[contig.df$AcMN, "novel.80"]))
nov.tab <- rbind(nov.tab, nov.tab.MN)
rownames(nov.tab) <- c("bit.50.all", "bit.80.all", "bit.50.highCA", "bit.80.highCA" )
print(xtable(nov.tab), floating=FALSE)



###################################################
### chunk number 32: 
###################################################
#line 3401 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
mf <- mf.dn.ds[, c(2, 5:8)]
bp <- bp.dn.ds[, c(2, 5:8)]
cc <- cc.dn.ds[, c(2, 5:8)]
nhlines <- c(0, nrow(mf), nrow(mf) + nrow(bp), nrow(mf) + nrow(bp) + nrow(cc))

onto.tab <- xtable(rbind(mf, bp, cc))

align(onto.tab) <- "rrrrp{4cm}r"

print(onto.tab, tabular.environment="longtable",
      floating=FALSE,
      hline.after=nhlines, include.rownames=FALSE)


