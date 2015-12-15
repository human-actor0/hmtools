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
* UCSC genome data (chromosome.fasta)

## Installation (linux/Mac)
```
# find a place to install  
>pwd
/vol1/home/hyunmin/
# download hmtools
>git clone https://github.com/actor0/hmtools.git
# or update hmtools if it exists
>cd /vol1/home/hyunmin/hmtools
>git pull origin master
# set HMHOME variable
# insert the following to ~/.bash_profile
HMHOME=/vol1/home/hyunmin/hmtools
export HMHOME
# deploy the setting
>source ~/.bash_profile
# check HMHOME is correct
echo $HMHOME
```
## Usage
[click here](https://github.com/actor0/hmtools/wiki/Home)
