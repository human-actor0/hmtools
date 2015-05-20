
echo \
"chr22	29706220	29706264	.	0	-"\
 | $HMHOME/src/seq.sh -s - $HMHOME/data/hg19_chr22.fa.gz 

echo \
"chr22	29706220	29706264	.	0	-"\
 | $HMHOME/src/seq.sh -s - /mnt/db/Ucsc/hg19/fasta/
