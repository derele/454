#! /bin/bash

######### Taking appart the mapping to analyse variants in
######### different libraries
######### Run from /home/ele/Data/454/mapping/mapping_each_population/
SAMT=/home/ele/tools/samtools/samtools

for i in $(ls final.*.sff.trimmed.sff.fasta) 
do 
    ##  map using ssaha2
    echo $i;
    ssaha2 -kmer 13 -skip 3 -seeds 6 -score 100 -cmatch 10 -ckmer 6 -output sam -best 1 -outfile $i.sam ../fullest_assembly_imputed.fasta  $i;
    ## to bam
    $SAMT view -uq1 -b -t ../fullest_assembly_imputed.fasta.fai -S $i.sam -o $i.bam;
    ## sort it
    $SAMT sort $i.bam $i.sorted;
    $SAMT pileup -f ../fullest_assembly_imputed.fasta $i.sorted.bam | java -jar /home/ele/tools/varscan/VarScan.v2.2.2.jar  pileup2snp --min-coverage 0 --min-reads2 0 --min-avg-qual 0 --min-var-freq 0 --p-value 100 > $i.varsnp;
done


$SAMT mpileup -uf ../fullest_assembly_imputed.fasta *.sorted.bam | /home/ele/tools/samtools/bcftools/bcftools view -gcv - > all_imp.vca