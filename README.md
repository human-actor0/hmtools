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
## Install via docker in Mac
 1. install docker (https://docs.docker.com/installation/mac/)
```
 > boot2docker init docker
 > boot2docker start
 > boot2docker shellinit
 > eval "$(boot2docker shellinit)"
```
 1. download hmtools image   
```
 > docker pull oannes/hmtools:v0.2.1

```
 1. start docker with a volume where your bam files and results will be saved
> imagine your working dir is /workingdir at local computer
bam files are in /workingdir/bams/control.bam  treatment.bam
chromosome files are in /workingdir/chromosome/chr1.fa.gz chr2.fa.gz ... 
```
 > docker run -it -v /workingdir/ oannes/hmtools:v0.2.1 /bin/bash --login
 # ls /workingdir/*  
```

## Run
# prepare a BATCH.txt file
```
FASTA=/workindir/chromosome/
BAM="
control /workingdir/bams/control.bam 
treatment /workingdir/bams/treatment.bam 
"
TARGET=$HMHOME/data/hg19_ensGene3utr.bed.gz

COMP="
treatment control 
"

## define minimum distance between peak centers
MDIST=20
## define output directory
OUT=/workdir/print_result_here
```
# run it now !!
```
 # hm batch_polya BATCH.txt
```

