HMTools
=======
>This is a collection of Hyun min's bioinformatics tools (human.gim@gmail.com)
>Focusing on sevral projcects:
- alternative polyadenylation
- alternative splicing
- patterns of RNA polymerase II (CTD) 

## Prerequisite 
* python: numpy, scipy
* BEDTools > v2.0
* Samtools > v1.2
* UCSC genome data
```
# download genome files in 
 > cd <genome_dir>
 > wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*.fa.gz'
```
* Or, just use a docker image [docker/hmtools](https://registry.hub.docker.com/u/oannes/hmtools/)

## Install on linux 
```
 > git https://github.com/actor0/hmtools.git
 > cp hmtools
 > hm # this will setup paths and enviromnet variables
 > source ~/.bash_profile
 > cd .. 
 > hm # this will guide you to the tools  
```

## Run
# prepare a BATCH file
```
## BATCH.txt example
FASTA=/mnt/db/Ucsc/hg19/chromosome/
## experimet name and bam files
BAM=(
        helaCtrAso_140128 /mnt/db/bams/polyaseq/HelapA_140128_CTRL_ASO_R1/a.bam
        helaU1Aso_140128 /mnt/db/bams/polyaseq/HelapA_140128_U1_ASO_R1/a.bam
)
TARGET=$HMHOME/data/hg19_ensGene3utr.bed.gz

## comparison plan: treatment vs control
COMP=(
        helaU1Aso_140128 helaCtrAso_140128
)

## define minimum distance between peak centers
MDIST=50
## define output directory
OUT=out
```
# run it now !!
```
 > hm batch_polya BATCH.txt
```

