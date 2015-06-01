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

## Run
# Imagine your files are :
```
/mydir/bams/control.bam control.bam.bai treatment.bam treatment.bam.bai
/mydir/chromosome/chr1.fa.gz chr2.fa.gz ... 
```
# Note
  * do not leave a space between "=", FASTA = .. is not working 
  * change /mydir/chromosome to  /mydir/chromosome/ 
  * make sure bai index files are in the same directory with bam files (use samtools index program to make indices ) 
  * make full paths ( do not use relative directory names )

# prepare a BATCH.txt file in the working directory /mydir/BATCH.txt 

```

FASTA=/mydir/chromosome/
BAM="
control /mydir/bams/control.bam 
treatment /mydir/bams/treatment.bam 
"
TARGET=$HMHOME/data/hg19_ensGene3utr.bed.gz

COMP="
treatment control 
"

## define minimum distance between peak centers
MDIST=20
## define output directory
OUT=/mydir/results
```

# run docker 
```
 > docker run -it -v /mydir/:/mydir/ oannes/hmtools:v0.2.1 /bin/bash --login
```
# check /mydir/ 
```
 # ls /mydir/*; 
```

# test and delete tested results
```
 # hm batch_polya /mydir/BATCH.txt test
 # rm -rf /mydir/results/ 
```

# run
```
 # hm batch_polya /mydir/BATCH.txt
```

