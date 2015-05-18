#!/bin/bash 
. $HM/functions/polya.sh
#samtools view -bq 10 ../data/hg19chr22sample.bam | bamToBed | modify_score - count | point -  > a.bed 
filter a.bed /mnt/db/Ucsc/hg19/ > b.bed
cluster b.bed 10 > c.bed

