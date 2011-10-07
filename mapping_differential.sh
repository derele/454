#! /bin/bash

######### Taking appart the mapping to analyse variants in
######### different libraries

######### Run from /home/ele/Data/454/mapping/mapping_each_population/

### Taiwan

##  map using ssaha2
ssaha2 -rtype 454 -output sam -best 1 -outfile taiwan_vs_fullest_imputed.sam ../fullest_assembly_imputed.fasta taiwan.fasta;

## into a bam
/home/ele/tools/samtools/samtools view -b -t ../fullest_assembly_imputed.fasta.fai -S taiwan_vs_fullest_imputed.sam -o taiwan_vs_fullest_imputed.bam;

## sort the bam
/home/ele/tools/samtools/samtools sort taiwan_vs_fullest_imputed.bam taiwan_vs_fullest_imputed.sorted;

## build pileup
/home/ele/tools/samtools/samtools pileup -f ../fullest_assembly_imputed.fasta taiwan_vs_fullest_imputed.sorted.bam > taiwan_vs_fullest_imputed.pileup;

## pileup to variant-call sensitive
cat taiwan_vs_fullest_imputed.pileup | java -jar /home/ele/tools/varscan/VarScan.v2.2.2.jar  pileup2snp --min-coverage 1 --min-reads2 1 --min-avg-qual 1 --min-var-freq 0 --p-value 1 > tmp.varsnp;

## Get only Contigs
grep "Contig" tmp.varsnp > taiwan_vs_fullest_imputed.varsnp;

## remove the file contining also hits to the singletons
rm tmp.varsnp;

cat taiwan_vs_fullest_imputed.pileup | java -jar /home/ele/tools/varscan/VarScan.v2.2.2.jar  pileup2indel  --min-coverage 1 --min-reads2 1 --min-avg-qual 1 --min-var-freq 0 --p-value 1 > tmp.varindl;

grep "Contig" tmp.varindl > taiwan_vs_fullest_imputed.varindl;

## remove the file contining also hits to the singletons
rm tmp.varindl;


## Europe

##  map using ssaha2
ssaha2 -rtype 454 -output sam -best 1 -outfile europe_vs_fullest_imputed.sam ../fullest_assembly_imputed.fasta europe.fasta;

## into a bam
/home/ele/tools/samtools/samtools view -b -t ../fullest_assembly_imputed.fasta.fai -S europe_vs_fullest_imputed.sam -o europe_vs_fullest_imputed.bam;

## sort the bam
/home/ele/tools/samtools/samtools sort europe_vs_fullest_imputed.bam europe_vs_fullest_imputed.sorted;

## build pileup
/home/ele/tools/samtools/samtools pileup -f ../fullest_assembly_imputed.fasta europe_vs_fullest_imputed.sorted.bam > europe_vs_fullest_imputed.pileup;

## pileup to variant-call sensitive
cat europe_vs_fullest_imputed.pileup | java -jar /home/ele/tools/varscan/VarScan.v2.2.2.jar  pileup2snp --min-coverage 1 --min-reads2 1 --min-avg-qual 1 --min-var-freq 0 --p-value 1 > tmp.varsnp;

## Get only Contigs
grep "Contig" tmp.varsnp > europe_vs_fullest_imputed.varsnp;

## remove the file contining also hits to the singletons
rm tmp.varsnp;

cat europe_vs_fullest_imputed.pileup | java -jar /home/ele/tools/varscan/VarScan.v2.2.2.jar  pileup2indel  --min-coverage 1 --min-reads2 1 --min-avg-qual 1 --min-var-freq 0 --p-value 1 > tmp.varindl;

grep "Contig" tmp.varindl > europe_vs_fullest_imputed.varindl;

## remove the file contining also hits to the singletons
rm tmp.varindl;