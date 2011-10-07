###################################################
### chunk number 1: load.libs
###################################################
#line 322 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
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

source("/home/ele/thesis/454/common_R_functions.R")


###################################################
### chunk number 2: trim
###################################################
#line 345 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
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
#line 387 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

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


###################################################
### chunk number 4: load.cobl
###################################################
#line 414 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

source("/home/ele/thesis/454/common_R_functions.R")

eelmRNA <- read.blast.best("/home/ele/Data/454/pre_assembly_screening/trimmed_all_vs_eelmRNA.blt")
eelmRNA$dbhit <- "eelmRNA"

eelrRNA <- read.blast.best("/home/ele/Data/454/pre_assembly_screening/trimmed_all_vs_eelrRNA.blt")
eelrRNA$dbhit <- "eelrRNA"

AcrRNA <- read.blast.best("/home/ele/Data/454/pre_assembly_screening/trimmed_all_vs_AcrRNA.blt")
AcrRNA$dbhit <- "AcrRNA"

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
#line 458 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
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
#line 474 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

pre.screen.tab <- table(R$lib, R$dbhit)
names(pre.screen.tab) <- c("A. crassus rRNA", "host", "host rRNA", "valid (assembled)")

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
#line 515 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
load("/home/ele/thesis/454/Method-assembly/Method.Rdata")

nMN <- nrow(contig.df[contig.df$category%in%"MN",])

nbad.c <- nrow(contig.df[contig.df$category%in%
                         c("M_n","N_n","N_1","M_1"),])

nSing <- nrow(contig.df[contig.df$category%in%0,])


###################################################
### chunk number 8: mean.cov
###################################################
#line 544 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

cov <- round(mean(con.pile[, "coverage"]),2)
cov.u <- round(mean(con.pile.uniq[, "uniq_coverage"]),2)
cov.mn <- round(mean(con.pile[con.pile$contig%in%
                                contig.df[contig.df$category=="MN", "contig"],
                                "coverage"]),2)
cov.mn.u <- round(mean(con.pile.uniq[con.pile.uniq$contig%in%
                                        contig.df[contig.df$category=="MN", "contig"],
                                        "uniq_coverage"]),2)



###################################################
### chunk number 9: host.screen
###################################################
#line 570 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

source("/home/ele/thesis/454/common_R_functions.R")

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


## nrow(contig.df[contig.df$contamination1=="valid" &
##                contig.df$contamination2%in%c("valid", "nempep4") &
##                contig.df$phylum.nr%in%c("Nematoda", "No_hit"), ])


## summary(contig.df[contig.df$contamination1=="valid" &
##                   contig.df$contamination2%in%c("valid", "nempep4") &
##                   contig.df$phylum.nr%in%c("Nematoda", "No_hit"), "uniq_coverage"])



###################################################
### chunk number 10: post.con
###################################################
#line 709 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
post.table <- rbind(table(contig.df$contamination),
                    rbind(tapply(contig.df$uniq_coverage, contig.df$contamination, mean)))

row.names(post.table) <- c("number",
                           "mean coverage")

nvalid <- nrow(contig.df[grepl("valid", contig.df$contamination) ,])

## Examples of selecting from the data

select_no_eel <- grepl("valid", contig.df$contamination) &
                            !contig.df$phylum.nt%in%"Chordata" &
                            !contig.df$phylum.nr%in%"Chordata"
select_no_eel_MN <- select_no_eel  & contig.df$category%in%"MN"

n_no_eel <- nrow(contig.df[select_no_eel,])
n_MNno_eel <- nrow(contig.df[select_no_eel_MN,])

select_AC <- grepl("valid", contig.df$contamination) &
                          contig.df$phylum.nt%in%c("Nematoda", "No_hit") &
                          contig.df$phylum.nr%in%c("Nematoda", "No_hit")
select_AC_MN <- select_AC & contig.df$category%in%"MN"

contig.df$Ac <- select_AC
contig.df$AcMN <- select_AC_MN

n_su_Ac <- nrow(contig.df[contig.df$Ac, ])
n_suMN_Ac <- nrow(contig.df[contig.df$AcMN, ])



###################################################
### chunk number 11: conservation
###################################################
#line 778 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
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
  B <- taxa.blast[taxa.blast$V12>bit.threshold, ]
  B <- subset(B, !level%in%"undef")
  if(nrow(B) ==0){return("No_hit")}
  if(all(B[, level]%in%value)){return(FALSE)}
  if(any(B[, level]%in%value) &
     any(!B[,level]%in%value)){return(TRUE)}
  if(all(!B[, level]%in%value) &
     any(!B[,level]%in%value)){return("dubious")}
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
                                                  NA))))

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
                                                  NA))))

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
#line 954 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

source("/home/ele/thesis/454/common_R_functions.R")


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

same.same <- summary.factor(tolower(as.character(contig.df$imp))==tolower(as.character(contig.df$seq)))


imp.pile <- read.table(pipe("cut -f1,2,4 /home/ele/Data/454/mapping/all_vs_full_imputed_uq.pileup"))
names(imp.pile) <- c("contig", "base", "imputed.coverage")
per.con.imp <- data.frame(imputed.coverage= tapply(imp.pile$imputed.coverage, imp.pile$contig, mean, rm.na=T))
contig.df <- merge(contig.df, per.con.imp, by.x="contig", by.y="row.names", all.x=TRUE)



###################################################
### chunk number 13: snp
###################################################
#line 1075 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

source("/home/ele/thesis/454/common_R_functions.R")  
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
#line 1618 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

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
#line 1677 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
varlog <- readLines("/home/ele/Data/454/mapping/all_vs_full_imputed_uq.varlog")
nRawSnps <- unlist(strsplit(varlog[length(varlog)], " "))[[1]]
nCov8bases <- nrow(imp.pile[imp.pile$imputed.coverage>7, ])
nbases <- nrow(imp.pile)
n.inOrf <- sum(contig.df$syn.sites, contig.df$nsyn.sites, na.rm=TRUE) 
n.outOrf <- n.inOrf-nbases

nReadSnps <- nrow(VAR)
nReadOrf <- tapply(VAR$inORF, VAR$inORF, length)
nReadFrame <- tapply(VAR$inFRAME, VAR$inFRAME, length)
tsv.read <- round(transversion.transition(VAR), 2)
tsv.read.orf <- round(do.call(rbind, by(VAR, VAR$inORF, transversion.transition)), 2)
read.dn.ds <- round(get.dn.ds(VAR), 2)

nScreeSnps <- nrow(VARq)
nScreeOrf <- tapply(VARq$inORF, VARq$inORF, length)
nScreeFrame <- tapply(VARq$inFRAME, VARq$inFRAME, length)
tsv.scree <- round(transversion.transition(VARq),2)
tsv.scree.orf <- round(do.call(rbind, by(VARq, VARq$inORF, transversion.transition)),2)
scree.dn.ds <- round(get.dn.ds(VARq),2)

nPercSnps <- nrow(VARqp)
nPercOrf <- tapply(VARqp$inORF, VARqp$inORF, length)
nPercFrame <- tapply(VARqp$inFRAME, VARqp$inFRAME, length)
tsv.perc <- round(transversion.transition(VARqp), 2)
tsv.perc.orf <- round(do.call(rbind, by(VARqp, VARqp$inORF, transversion.transition)), 2)
perc.dn.ds <- round(get.dn.ds(VARqp), 2)

snp.summary <- data.frame(c("No.SNPs", nReadSnps, nScreeSnps, nPercSnps),
                          c("in.ORF",nReadOrf[["TRUE"]], nScreeOrf[["TRUE"]], nPercOrf[["TRUE"]]),
                          check.names=FALSE)

FR <- rbind(nReadFrame, nScreeFrame, nPercFrame)
snp.summary <- cbind(snp.summary, rbind(colnames(FR), FR))

TS <- rbind(tsv.read[["ratioTS.TV"]], tsv.scree[["ratioTS.TV"]], tsv.perc[["ratioTS.TV"]])
snp.summary <- cbind(snp.summary, rbind("overall", TS))

TSf <- rbind(tsv.read.orf[["TRUE", "ratioTS.TV"]],
             tsv.scree.orf[["TRUE", "ratioTS.TV"]],
             tsv.perc.orf[["TRUE", "ratioTS.TV"]])

snp.summary <- cbind(snp.summary, rbind("ins.orf", TSf))

TSo <- rbind(tsv.read.orf[["FALSE", "ratioTS.TV"]],
             tsv.scree.orf[["FALSE", "ratioTS.TV"]],
             tsv.perc.orf[["FALSE", "ratioTS.TV"]])

snp.summary <- cbind(snp.summary, rbind("outs.orf", TSo))
snp.summary <- cbind(snp.summary, rbind("dn.ds", read.dn.ds, scree.dn.ds, perc.dn.ds))

names(snp.summary) <- c(" ", " ", "pos", "in", "codon", " ", "ti/tv ", " ", " ")

rownames(snp.summary) <- c("", "raw", "h.screened", "p.screened")

read.per.b <- nReadSnps/(nbases/1000)
scree.per.b <- nScreeSnps/(n.inOrf/1000)
perc.per.b <- nPercSnps/(n.inOrf/1000)

total.syn <- sum(contig.df[, "syn.sites"], na.rm=TRUE)
total.cov8.syn <- sum(contig.df[, "cov8.syn.sites"], na.rm=TRUE)

total.nsyn <- sum(contig.df[ , "nsyn.sites"], na.rm=TRUE )
total.cov8.nsyn <- sum(contig.df[ , "cov8.nsyn.sites"], na.rm=TRUE )

total.nSNP <- sum(contig.df$Nonsynonymous, na.rm=TRUE)
total.sSNP <- sum(contig.df$Synonymous, na.rm=TRUE)

s.per.s.base <- total.sSNP/(total.syn/1000)
s.per.s.cov8.base <- total.sSNP/(total.cov8.syn/1000)

n.per.n.base <- total.nSNP/(total.nsyn/1000)
n.per.n.cov8.base <- total.nSNP/(total.cov8.nsyn/1000)



###################################################
### chunk number 16: annolibs
###################################################
#line 1853 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
library(GOstats)
library(GO.db)
## library(AnnotationDbi)
## library(GSEABase)
## library(GeneAnswers)
## library(goTools)
## For comparison with C.elegans
## require("org.Ce.eg.db")


###################################################
### chunk number 17: annotation
###################################################
#line 1864 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
GO.annot <- read.delim("/home/ele/Data/454/annotation/annot8r/output/GO.csv", sep=",", header=FALSE)
names(GO.annot) <- c("pept_id", "go_term", "pcf", "descr", "slim", "besthit", "bestscore", "bestev", "hitnum", "maxhits", "fraction") 

## reset the shorted names used in annotation to the proper contig-names
GO.annot$pept_id <- gsub("Ac_", "Acrassus_", GO.annot$pept_id)


GO.annot$go_term <- as.character(gsub(" ", "", GO.annot$go_term))
GO.annot$pcf <- gsub(" C", "CC", GO.annot$pcf)
GO.annot$pcf <- gsub(" F", "MF", GO.annot$pcf)
GO.annot$pcf <- gsub(" P", "BP", GO.annot$pcf)

EC.annot <- read.delim("/home/ele/Data/454/annotation/annot8r/output/EC.csv", sep=",", header=FALSE)
names(EC.annot) <- c("pept_id", "ec_id", "descr", "besthit", "bestscore", "bestev", "hitnum", "maxhits", "fraction") 
EC.annot$pept_id <- gsub("Ac_", "Acrassus_", EC.annot$pept_id)

KEGG.annot <- read.delim("/home/ele/Data/454/annotation/annot8r/output/KEGG.csv", sep=",", header=FALSE)
names(KEGG.annot) <- c("pept_id", "ko_id", "path", "descr", "besthit", "bestscore", "bestev", "hitnum", "maxhits", "fraction") 
KEGG.annot$pept_id <- gsub("Ac_", "Acrassus_", KEGG.annot$pept_id)



###################################################
### chunk number 18: sigp
###################################################
#line 1887 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

signalip <- as.data.frame(read.delim("/home/ele/Data/454/annotation/signalp/Ac_p4ePro.signalip",
                                     skip=1, sep="\t", header=FALSE,
                                     strip.white=TRUE, comment.char=";",
                                     as.is=TRUE))

signalnn.names <- gsub("# ", "", signalip[1, 1])
signalnn.names <- unlist(strsplit(signalnn.names, " +"))

signalhmm.names <- gsub("# ", "", signalip[1, 2])
signalhmm.names <- unlist(strsplit(signalhmm.names, " +"))


signalip <- signalip[!grepl("#", signalip[, 1]) , ]

signalnn <- as.data.frame(do.call(rbind,
                                  strsplit(as.character(signalip[, 1]), " +")))
names(signalnn) <- signalnn.names

signalhmm <- as.data.frame(do.call(rbind,
                                  strsplit(as.character(signalip[, 2]), " +")))
names(signalhmm) <- signalhmm.names

sig.sum <- cbind(signalnn[ , c(1, 14)], signalhmm$"!")
names(sig.sum) <- c("contig", "sigp.nn", "sigp.hmm")

## translate the contig-names back
sig.sum$contig <- gsub("Ac_", "Acrassus_", sig.sum$contig)

contig.df <- merge(contig.df, sig.sum, all.x=TRUE)



###################################################
### chunk number 19: annot
###################################################
#line 1921 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
GOannotations <- unique(GO.annot$pept_id)
nGO <- length(GOannotations)

ECannotations <- unique(EC.annot$pept_id)
nEC <- length(ECannotations)

KEGGannotations <- unique(KEGG.annot$pept_id)
nKEGG <- length(KEGGannotations)

## IPRannotations <- unique(IPR.annot[grepl("GO:", IPR.annot$GO_descr), "pept_id"])
## nIPR <- length(IPRannotations)
## nIPRonly <- length(IPRannotations[!IPRannotations%in%GOannotations &
##                                   !IPRannotations%in%ECannotations &
##                                   !IPRannotations%in%KEGGannotations
##                                   ])

Allannotations <- unique(c(GOannotations, ECannotations, KEGGannotations))#,  IPRannotations))
nallannotations <- length(Allannotations)

venn.diagram(list(GO   = match(GOannotations, Allannotations),
                  EC   = match(ECannotations, Allannotations),
                  KEGG = match(KEGGannotations, Allannotations)),
                  ##                  IPR  = match(IPRannotations, Allannotations)),
                  filename = "/home/ele/thesis/454/figures/annotataionVenn.tiff")


Signn <- contig.df[contig.df$sigp.nn%in%"Y", "contig" ]
Sighmm <- contig.df[contig.df$sigp.hmm%in%"S", "contig" ]
Anchmm <- contig.df[contig.df$sigp.hmm%in%"A", "contig" ]

nSignn <- length(Signn)
nSighmm <- length(Sighmm)
nAnchmm <- length(Anchmm)

venn.diagram(list(nn.Signal   = match(Signn, contig.df$contig),
                  hmm.Signal   = match(Sighmm, contig.df$contig),
                  hmm.Anchor = match(Anchmm, contig.df$contig)),
                  filename = "/home/ele/thesis/454/figures/signalVenn.tiff")



###################################################
### chunk number 20: annot.compare
###################################################
#line 1971 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

Bm.sum.annot.C <- read.delim("/home/ele/Data/454/annotation/annot8r/output_Bm/piedata_C", skip=1)
Bm.sum.annot.C$species <- "Brugia malayi"
Ac.sum.annot.C <- read.delim("/home/ele/Data/454/annotation/annot8r/output/piedata_C", skip=1)
Ac.sum.annot.C$species <- "Anguillicola crassus"
sum.annot.C <- rbind(Bm.sum.annot.C, Ac.sum.annot.C)
sum.annot.C$ontology <- "Cellular compartment"

Bm.sum.annot.F <- read.delim("/home/ele/Data/454/annotation/annot8r/output_Bm/piedata_F", skip=1)
Bm.sum.annot.F$species <- "Brugia malayi"
Ac.sum.annot.F <- read.delim("/home/ele/Data/454/annotation/annot8r/output/piedata_F", skip=1)
Ac.sum.annot.F$species <- "Anguillicola crassus"
sum.annot.F <- rbind(Bm.sum.annot.F, Ac.sum.annot.F)
sum.annot.F$ontology <- "Molecular Function"

Bm.sum.annot.P <- read.delim("/home/ele/Data/454/annotation/annot8r/output_Bm/piedata_P", skip=1)
Bm.sum.annot.P$species <- "Brugia malayi"
Ac.sum.annot.P <- read.delim("/home/ele/Data/454/annotation/annot8r/output/piedata_P", skip=1)
Ac.sum.annot.P$species <- "Anguillicola crassus"
sum.annot.P <- rbind(Bm.sum.annot.P, Ac.sum.annot.P)
sum.annot.P$ontology <- "Biological process"

sum.annot <- rbind(sum.annot.C, sum.annot.F, sum.annot.P)



GO.bm.com <- ggplot(sum.annot, aes(x=description, fill=species, weight=occurences)) +
  geom_bar(position="dodge") +
  scale_x_discrete(breaks=sum.annot$description,
                   label=gsub("and\nnucleic\nacid\nmetabolic\nprocess",
                     "and nucleic acid\nmetabolic process", gsub(" ", "\n", sum.annot$description))) +
  facet_wrap(~ontology, ncol=1, scales="free") +
    opts(legend.text = theme_text(face="italic"))

ggsave("/home/ele/thesis/454/figures/go_bm_com.png", GO.bm.com)


###################################################
### chunk number 21: overrep
###################################################
#line 2023 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

## subsets for contig.df and GO-annotation containing only dn.ds contigs
dn.ds.df <- subset(contig.df, !is.na(dn.ds))
GO.dn.ds <- subset(GO.annot, pept_id%in%dn.ds.df$contig)

## the category of >.5 dn/ds contigs
big05 <- dn.ds.df[dn.ds.df$dn.ds>0.5 &
                  dn.ds.df$contig%in%GO.dn.ds$pept_id, "contig"]

## whole data-frames for these categories
big05.df <- contig.df[contig.df$dn.ds>0.5 & !is.na(contig.df$dn.ds), ]
EC.big05.df <- merge(EC.annot, big05.df, by.x="pept_id", by.y="contig")

## now the same for the contigs only absolutely surely  Ac
dn.ds.df.AC <- dn.ds.df[dn.ds.df$Ac, ]

GO.dn.ds.AC <- subset(GO.annot, pept_id%in%dn.ds.df.AC$contig)

big05.AC <- dn.ds.df.AC[dn.ds.df.AC$dn.ds>0.5 &
                        dn.ds.df.AC$contig%in%GO.dn.ds.AC$pept_id, "contig"]


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

## using just what is not na in dn ds as universe
mf.dn.ds <- test.over.under.GOstats(GO.dn.ds, big05, "MF")
bp.dn.ds <- test.over.under.GOstats(GO.dn.ds, big05, "BP")
cc.dn.ds <- test.over.under.GOstats(GO.dn.ds, big05, "CC")

## using just what is not na in dn ds as universe BUT
## only for sure AC this time
mf.dn.ds.AC <- test.over.under.GOstats(GO.dn.ds.AC, big05.AC, "MF")
bp.dn.ds.AC <- test.over.under.GOstats(GO.dn.ds.AC, big05.AC, "BP")
cc.dn.ds.AC <- test.over.under.GOstats(GO.dn.ds.AC, big05.AC, "CC")

## make a list of offspring terms to allow finding of 
GOMFOFF.list <- as.list(GOMFOFFSPRING)
GOMFOFF.list <- GOMFOFF.list[!is.na(GOMFOFF.list)]

## which are the amino acid transporter contigs (2) being in dn.ds>0.5
amino.acid.transporter.go <- mf.dn.ds.AC[grepl("acid transmembrane transporter",
                                               mf.dn.ds.AC$Term), 1]
amino.acid.tr.off <- unlist(GOMFOFF.list[amino.acid.transporter.go])

## get the term-contigs plus their boffspring term-contigs
## consistent with mf.dn.ds.AC results (size)
amino.acid.tr.contigs <- unique(GO.dn.ds.AC[GO.dn.ds.AC$go_term%in%
                                            c(amino.acid.transporter.go, amino.acid.tr.off),
                                            "pept_id"])

peptidase.go <- mf.dn.ds.AC[grepl("peptidase", mf.dn.ds.AC$Term), 1]
peptidase.off <- unlist(GOMFOFF.list[peptidase.go])

## consistent with mf.dn.ds.AC results
peptidase.contigs <- unique(GO.dn.ds.AC[GO.dn.ds.AC$go_term%in%
                                        c(peptidase.go, peptidase.off),
                                            "pept_id"])

## consistent with mf.dn.ds.AC results (count)
peptidase.contigs.pos <- contig.df[contig.df$contig%in%peptidase.contigs &
                                   contig.df$dn.ds>0.5 , "contig"]


## just to test that results for quality-contigs only are consistent.
mf.dn.ds.AC.MN <- test.over.under.GOstats(subset(GO.dn.ds.AC,
                                                 GO.dn.ds.AC$pept_id%in%
                                                 contig.df[contig.df$category%in%"MN",
                                                           "contig"]),
                                          big05.AC, "MF")

bp.dn.ds.AC.MN <- test.over.under.GOstats(subset(GO.dn.ds.AC,
                                                 GO.dn.ds.AC$pept_id%in%
                                                 contig.df[contig.df$category%in%"MN",
                                                           "contig"]),
                                          big05.AC, "BP")

cc.dn.ds.AC.MN <- test.over.under.GOstats(subset(GO.dn.ds.AC,
                                                 GO.dn.ds.AC$pept_id%in%
                                                 contig.df[contig.df$category%in%"MN",
                                                           "contig"]),
                                          big05.AC, "CC")

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


## this for the easy nn-data
u.test <- wilcox.test(dn.ds ~ sigp.nn, data=dn.ds.df)

u.test.AC <- wilcox.test(dn.ds ~ sigp.nn, data=dn.ds.df.AC)

## ## kruskal-Wallis tests, but they say nothing about kontrasts
## ## in group comparisons
## kw.test <- kruskal.test(dn.ds ~ sigp.hmm, data=dn.ds.df)

## ### Kruskal-Wallis test, approximate exact p-value
## kw.test2 <- kruskal_test(dn.ds ~ sigp.hmm, data = dn.ds.df, 
##                    distribution = approximate(B = 9999))
     
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

ndwd <- get.ndwd(dn.ds.df, dn.ds.df$dn.ds, dn.ds.df$sigp.hmm)
kw.ph <- pvalue(ndwd, method = "single-step")


ndwd.AC <- get.ndwd(dn.ds.df.AC, dn.ds.df.AC$dn.ds, dn.ds.df.AC$sigp.hmm)
kw.ph.AC <- pvalue(ndwd.AC, method = "single-step")

t.test.hmm <- wilcox.test(dn.ds.df$dn.ds ~ as.logical(dn.ds.df$sigp.hmm%in%"S"))
t.test.AC.hmm <- wilcox.test(dn.ds.df.AC$dn.ds ~ as.logical(dn.ds.df.AC$sigp.hmm%in%"S"))







###################################################
### chunk number 22: nn.dn.ds
###################################################
#line 2215 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

nn.box <- ggplot(dn.ds.df, aes(dn.ds.df$sigp.nn, dn.ds.df$dn.ds)) + 
  geom_boxplot() + 
  scale_y_log10("dn/ds") +
  scale_x_discrete("", 
                   breaks=levels(dn.ds.df$sigp.nn),
                   labels=c("No signal peptide",
                     "Signal peptide")) +
  opts(title = "SignalP neural networks prediction",
       axis.title.x = theme_text(size = 15), 
       axis.text.x = theme_text(size = 15))

  
hmm.box <- ggplot(dn.ds.df, aes(dn.ds.df$sigp.hmm, dn.ds.df$dn.ds)) + 
  geom_boxplot() + 
  scale_y_log10("dn/ds") +
  scale_x_discrete("", 
                   breaks=levels(dn.ds.df$sigp.hmm),
                   labels=c("Signal anchor",
                     "No signal peptide",
                     "Signal peptide")) +
  opts(title = "SignalP hmm prediction",
       axis.title.x = theme_text(size = 15), 
       axis.text.x = theme_text(size = 15))
dev.off()

png("../figures/sigp_dn_ds.png", width=2000, height=1000, res=144)

# Set up the page
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)

# Make each plot, in the correct location
print(nn.box, vp = vplayout(1, 1 ))
print(hmm.box, vp = vplayout(1, 2 ))

dev.off() 


###################################################
### chunk number 23: novel.dn.ds
###################################################
#line 2301 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
novel.dn.ds.50 <- ggplot(subset(contig.df, !is.na(novel.50) & !is.na(dn.ds)),
                         aes(novel.50, dn.ds)) +
  geom_boxplot() +
  scale_y_log10("dn/ds") +
  scale_x_discrete("evolutionary conseravation at bitscore threshold of 50")

novel.dn.ds.80 <- ggplot(subset(contig.df, !is.na(novel.80) & !is.na(dn.ds)),
                                aes(novel.80, dn.ds)) +
  geom_boxplot() +
  scale_y_log10("dn/ds") +
  scale_x_discrete("evolutionary conseravation at bitsore threshold of 80")

# Set up the page

png("../figures/conservation_dn_ds.png", width=2000, height=1000, res=144)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)

# Make each plot, in the correct location
print(novel.dn.ds.50, vp = vplayout(1, 1 ))
print(novel.dn.ds.80, vp = vplayout(1, 2 ))

dev.off() 



###################################################
### chunk number 24: novel.dn.ds.test
###################################################
#line 2330 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
dn.ds.df$novel.50 <- as.factor(dn.ds.df$novel.50)
dn.ds.df$novel.80 <- as.factor(dn.ds.df$novel.80)

novel.50.ndwd <- get.ndwd(dn.ds.df, dn.ds.df$dn.ds, dn.ds.df$novel.50)
novel.50.kw.ph <- pvalue(novel.50.ndwd, method = "single-step")

novel.80.ndwd <- get.ndwd(dn.ds.df, dn.ds.df$dn.ds, dn.ds.df$novel.80)
novel.80.kw.ph <- pvalue(novel.80.ndwd, method = "single-step")


u.novel.50 <- wilcox.test(dn.ds.df$dn.ds~
                          as.factor(dn.ds.df$novel.50%in%"novel.in.clade3"))
u.novel.80 <- wilcox.test(dn.ds.df$dn.ds~
                          as.factor(dn.ds.df$novel.80%in%"novel.in.clade3"))



###################################################
### chunk number 25: novel.nn
###################################################
#line 2359 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

nn.50 <- table(contig.df$novel.50, contig.df$sigp.nn)
nn.80 <- table(contig.df$novel.80, contig.df$sigp.nn)

hmm.50 <- table(contig.df$novel.50, contig.df$sigp.hmm)
hmm.80 <- table(contig.df$novel.80, contig.df$sigp.hmm)


nn.50.p <- ggplot(subset(contig.df, !is.na(novel.50) & !is.na(sigp.nn)),
              aes(novel.50, fill=sigp.nn)) +
       geom_bar(position="fill") + coord_polar(theta="y") +
  scale_y_continuous("evolutionary conseravation categories at bitsore threshold of 50") +
  scale_x_discrete("") +
  opts(title="proportion of sequences in SignalP-nn category")


nn.80.p <- ggplot(subset(contig.df, !is.na(novel.80) & !is.na(sigp.nn)),
              aes(novel.80, fill=sigp.nn)) +
  geom_bar(position="fill") + coord_polar(theta="y") +
  scale_y_continuous("evolutionary conseravation categories at bitsore threshold of 80") +
  scale_x_discrete("") +
  opts(title="proportion of sequences in SignalP-nn category")

hmm.50.p <- ggplot(subset(contig.df, !is.na(novel.50) & !is.na(sigp.hmm)),
                   aes(novel.50, fill=sigp.hmm)) +
  geom_bar(position="fill") + coord_polar(theta="y") +
  scale_y_continuous("evolutionary conseravation categories at bitsore threshold of 50") +
  scale_x_discrete("") +
  opts(title="proportion of sequences in SignalP-hmm category")

hmm.80.p <- ggplot(subset(contig.df, !is.na(novel.80) & !is.na(sigp.hmm)),
                   aes(novel.80, fill=sigp.hmm)) +
  geom_bar(position="fill") + coord_polar(theta="y") +
  scale_y_continuous("evolutionary conseravation categories at bitsore threshold of 80") +
  scale_x_discrete("") +
  opts(title="proportion of sequences in SignalP-hmm category")

png("../figures/signal_novel.png", width=2000, height=2000, res=144)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)

# Make each plot, in the correct location
print(nn.50.p, vp = vplayout(1, 1 ))
print(nn.80.p, vp = vplayout(1, 2 ))
print(hmm.50.p, vp = vplayout(2, 1 ))
print(hmm.80.p, vp = vplayout(2, 2 ))

dev.off() 


###################################################
### chunk number 26: gene.expr
###################################################
#line 2422 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

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

n454.map <- length(unique(SAM$read))
n454.Umap <- length(unique(SAM.easy$read))

countsTable <- as.data.frame.matrix(table(SAM.easy$mapped, SAM.easy$lib))

## just plain counting the twice hitting
## countsTable.all <- as.data.frame.matrix(table(SAM$mapped, SAM.all$lib))

## Try to exclude contigs before DESeq
## countsTable2 <- countsTable[grep("Contig", row.names(countsTable)),]

conds.eel <- factor(c("EU", "EU", "L2", "M", "TW", "TW"))

good.contigs <- contig.df[contig.df$category%in%"MN", "contig"]

cds.eel <- newCountDataSet(countsTable[rownames(countsTable)%in%good.contigs,], conds.eel)
                           
cds.eel <- estimateSizeFactors(cds.eel)
cds.eel <- estimateVarianceFunctions(cds.eel)
res.eel <- nbinomTest(cds.eel, "TW", "EU")

resSig.eel <- subset(res.eel, padj < .2)
sig.eel.contigs <- resSig.eel$id
sig.eel.counts <- countsTable[sig.eel.contigs, ]

### This shows nothing as depth is not big enough to allow
### detection of significance
## plot( 
##      res.eel$baseMean, 
##      res.eel$log2FoldChange, 
##      log="x", pch=20, cex=.4, 
##      col = ifelse( res.eel$padj < .1, "red", "black" ),
##      main="My data")

conds.mf <- factor(c("F", "F", "L2", "M", "F", "F"))
cds.mf <- newCountDataSet(countsTable[rownames(countsTable)%in%good.contigs,], conds.mf)
cds.mf <- estimateSizeFactors(cds.mf)
cds.mf <- estimateVarianceFunctions(cds.mf)
res.mf <- nbinomTest(cds.mf, "M", "F")

resSig.mf <- subset(res.mf, padj < .1)
sig.mf.contigs <- resSig.mf$id
sig.mf.counts <- countsTable[sig.mf.contigs, ]

male.nem.hits <- nrow(subset(contig.df, contig.df$contig%in%sig.mf.contigs &
                             phylum.nr%in%"Nematoda"))

conds.ad <- factor(c("Ad", "Ad", "L2", "Ad", "Ad", "Ad"))
cds.ad <- newCountDataSet(countsTable, conds.ad)
cds.ad <- estimateSizeFactors(cds.ad)
cds.ad <- estimateVarianceFunctions(cds.ad)
res.ad <- nbinomTest(cds.ad, "Ad", "L2")

resSig.ad <- subset(res.ad, padj < .1)
sig.ad.contigs <- resSig.ad$id
sig.ad.counts <- countsTable[sig.ad.contigs, ]

l2.met.hits <- nrow(subset(contig.df, contig.df$contig%in%sig.ad.contigs & kingdom.nr%in%"Metazoa"))
l2.nem.hits <- nrow(subset(contig.df, contig.df$contig%in%sig.ad.contigs & phylum.nr%in%"Nematoda"))

vsd <- getVarianceStabilizedData( cds.mf )
dists <- dist( t( vsd ) )
idists <- as.matrix(dists)
## heatmap (idists , symm=TRUE, margins = c (7,7))



###################################################
### chunk number 27: sol.comp
###################################################
#line 2524 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

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

## tag.seq <- sapply(as.list(bam[[1]]$seq), toString)

counts <- merge(counts, countsTable, by.x="contig", by.y="row.names")
counts$all.reads <- apply(countsTable, 1, sum)

#better use
# plotmatrix
## exp_plot <- ggplot(counts, aes(y=solexa.tags+1, x=all.reads+1)) +
##   geom_point(alpha=0.3) + geom_smooth(method="lm") + scale_x_log10() +
##   scale_y_log10()

contig.df <- merge(contig.df, counts)

cor.count.tab <- round(cor(contig.df[,
                                     c("solexa.tags", "E1", "E2", "L2", "M", "T1", "T2", "all.reads")],
                           method="spearman"),4)


cor.count.tab.MN <- round(cor(contig.df[contig.df$category%in%"MN",
                                        c("solexa.tags", "E1", "E2", "L2",  "M", "T1", "T2", "all.reads")],
                              method="spearman"),4)


cor.count.tab.MN.AC <- round(cor(contig.df[contig.df$AcMN,
                                           c("solexa.tags", "E1", "E2", "L2",  "M", "T1", "T2", "all.reads")],
                                 method="spearman"),4)



###################################################
### chunk number 28: save
###################################################
#line 2639 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
## for a more readable version of contig.df
readable <- !names(contig.df)%in%c("seq", "imp", "ontology", "coding")
save.image("paper.Rdata")
write.csv(contig.df, "../A_crassus_contigs_full.csv")
write.csv(contig.df[,readable], "../A_crassus_contigs_readable.csv")



###################################################
### chunk number 29: 
###################################################
#line 3181 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
print(xtable(TRIM, display=rep("d", times=7)), floating=FALSE)


###################################################
### chunk number 30: 
###################################################
#line 3186 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
print(xtable(pre.screen.tab), floating=FALSE)


###################################################
### chunk number 31: 
###################################################
#line 3192 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
print(xtable(cerco), include.rownames=FALSE, floating=FALSE)


###################################################
### chunk number 32: 
###################################################
#line 3198 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
digs <- matrix(c(0, 2), nrow=2 , ncol=6)
print(xtable(post.table, digits=digs), floating=FALSE)


###################################################
### chunk number 33: 
###################################################
#line 3204 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
nov.tab <- rbind(summary.factor(contig.df$novel.50), 
                 summary.factor(contig.df$novel.80))[, c(1,3,4,2)]
rownames(nov.tab) <- c("bit.threshold.50", "bit.threshold.80")
print(xtable(nov.tab), floating=FALSE)



###################################################
### chunk number 34: 
###################################################
#line 3214 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
print(xtable(imp.sum),
      floating=FALSE)


###################################################
### chunk number 35: 
###################################################
#line 3222 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

print(xtable(snp.summary),
      hline.after = c(-1,1),
      floating=FALSE)



###################################################
### chunk number 36: 
###################################################
#line 3233 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"
mf <- mf.dn.ds.AC[, 5:8]
bp <- bp.dn.ds.AC[, 5:8]
cc <- cc.dn.ds.AC[, 5:8]
nhlines <- c(0, nrow(mf), nrow(mf) + nrow(bp), nrow(mf) + nrow(bp) + nrow(cc))

onto.tab <- xtable(rbind(mf, bp, cc))

print(onto.tab, tabular.environment="longtable",
      floating=FALSE,
      hline.after=nhlines, include.rownames=FALSE)


###################################################
### chunk number 37: 
###################################################
#line 3253 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

print(xtable(cor.count.tab, digits=3), floating=FALSE)



###################################################
### chunk number 38: 
###################################################
#line 3260 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

print(xtable(cor.count.tab.MN, digits=3), floating=FALSE)



###################################################
### chunk number 39: 
###################################################
#line 3268 "/home/ele/thesis/454/paper/transcriptome_paper.Rnw"

print(xtable(cor.count.tab.MN.AC, digits=3), floating=FALSE)



