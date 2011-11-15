#!/bin/bash

### use sametools faidx first to create ineexed fasta files for both mapped against ###

#################### For coverage:                             ####################################
#################### Mapping vs. the fullest_assembly          #####################################

##  map using ssaha2
ssaha2 -kmer 13 -skip 3 -seeds 6 -score 100 -cmatch 10 -ckmer 6  -output sam -best 1 -outfile all_vs_full.sam fullest_assembly.fasta all_assembled_reads.fasta

## into a bam
/home/ele/tools/samtools/samtools view -b -t fullest_assembly.fasta.fai -o all_vs_full.bam  -S all_vs_full.sam;

## sort the bam
/home/ele/tools/samtools/samtools sort all_vs_full.bam all_vs_full.sorted;

## build pileup
/home/ele/tools/samtools/samtools pileup -f fullest_assembly.fasta all_vs_full.sorted.bam > all_vs_full.pileup;

########################## unique mapping #########################################
#rev all_vs_full.sam | uniq -u -f 12 | rev > all_vs_full_uniq.sam

## into a bam, with samtool's native unique mapping -uq
/home/ele/tools/samtools/samtools view -uq1 -b -t fullest_assembly.fasta.fai -o all_vs_full_uniq.bam  -S all_vs_full.sam;

## sort the bam
/home/ele/tools/samtools/samtools sort all_vs_full_uniq.bam all_vs_full_uniq.sorted;

## build pileup
/home/ele/tools/samtools/samtools pileup -f fullest_assembly.fasta all_vs_full_uniq.sorted.bam > all_vs_full_uniq.pileup;


#################### For SNP-calling only !!!!!!!!!!!!!!!!!!! ####################################
#################### Mapping vs. the fullest_assembly_imputed #####################################

##  map using ssaha2
ssaha2 -kmer 13 -skip 3 -seeds 6 -score 100 -cmatch 10 -ckmer 6  -output sam -best 1 -outfile all_vs_full_imputed.sam  fullest_assembly_imputed.fasta all_assembled_reads.fasta

## into a bam
/home/ele/tools/samtools/samtools view -b -t fullest_assembly_imputed.fasta.fai -S all_vs_full_imputed.sam -o all_vs_full_imputed.bam;

## sort the bam
/home/ele/tools/samtools/samtools sort all_vs_full_imputed.bam all_vs_full_imputed.sorted;

## build pileup
/home/ele/tools/samtools/samtools pileup -f fullest_assembly_imputed.fasta all_vs_full_imputed.sorted.bam > all_vs_full_imputed.pileup;

## pileup to variant-call sensitive
cat all_vs_full_imputed.pileup | java -jar /home/ele/tools/varscan/VarScan.v2.2.2.jar  pileup2snp > all_vs_full_imputed.varsnp 2>all_vs_full_imputed.varlog;


## and the same uniquely
##  map using ssaha2

## into a bam
/home/ele/tools/samtools/samtools view -uq1 -b -t fullest_assembly_imputed.fasta.fai -S all_vs_full_imputed.sam -o all_vs_full_imputed_uq.bam;

## sort the bam
/home/ele/tools/samtools/samtools sort all_vs_full_imputed_uq.bam all_vs_full_imputed_uq.sorted;

## build pileup
/home/ele/tools/samtools/samtools pileup -f fullest_assembly_imputed.fasta all_vs_full_imputed_uq.sorted.bam > all_vs_full_imputed_uq.pileup;

## pileup to variant-call sensitive
cat all_vs_full_imputed_uq.pileup | java -jar /home/ele/tools/varscan/VarScan.v2.2.2.jar  pileup2snp > all_vs_full_imputed_uq.varsnp 2>all_vs_full_imputed_uq.varlog;


## using mpileup instead
/home/ele/tools/samtools/samtools mpileup -uf fullest_assembly_imputed.fasta all_vs_full_imputed_uq.sorted.bam | /home/ele/tools/samtools/bcftools/bcftools view -gcv - > all_once_imp.vca


