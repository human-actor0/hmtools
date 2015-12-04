. $HMHOME/src/bed.sh
. $HMHOME/src/seq.sh




memset.score5(){
	cmd=`cat $HMHOME/lib/MEMset/score5.pl`
	cmd=${cmd/INPUT1/$HMHOME/lib/MEMset/me2x5}
	cmd=${cmd/INPUT2/$HMHOME/lib/MEMset/splicemodels}
	cat $1 | perl -e "$cmd" 
}
memset.score3(){
	cmd=`cat $HMHOME/lib/MEMset/score3.pl`
	cmd=${cmd/INPUT1/$HMHOME/lib/MEMset/me2x5}
	cmd=${cmd/INPUT2/$HMHOME/lib/MEMset/splicemodels}
	cat $1 | perl -e "$cmd" 
}

memset.run(){
usage="$FUNCNAME <intron.bed> <fasta>"
if [ $# -ne 2 ];then echo "$usage"; return; fi
	local tmpd=`mymktempd`
	#local tmpd=tmpd;rm -rf $tmpd; mkdir -p $tmpd
	cat $1 | awk -v OFS="\t" '{ $4=$1";"$2";"$3";"$4";"$5";"$6;}1'> $tmpd/a
	bed.5p $tmpd/a | bed.flank - 3 5 -s | seq.read $2 - -s -uc \
		| awk -v OFS="\t" '$7 != "NULL" { print $7,$4;}' | sort -k 1,1  > $tmpd/b
	bed.3p $tmpd/a | bed.flank - 19 3 -s | seq.read $2 - -s -uc  \
		| awk -v OFS="\t" '$7 != "NULL" { print $7,$4;}' | sort -k 1,1 > $tmpd/c

	cut -f 1 $tmpd/b | sort -u | sort -k1,1 | memset.score5 - > $tmpd/bb
	cut -f 1 $tmpd/c | sort -u | sort -k1,1 | memset.score3 - > $tmpd/cc

	join -a 1 -e NA $tmpd/b $tmpd/bb | sort -k2,2 > $tmpd/bbb
	join -a 1 -e NA $tmpd/c $tmpd/cc | sort -k2,2 > $tmpd/ccc
	join -j 2 $tmpd/bbb $tmpd/ccc | tr "; " "\t"
	rm -rf $tmpd
}
memset.run.test(){
echo "chr22   17265299        17280660        ENSG00000172967 0       -
chr22   17280914        17288628        ENSG00000172967 0       -
chr2    217280752       217280979       ENSG00000138375 0       +
chr22   17288973        17302496        ENSG00000172967 0       -" \
| perl -ne 'chomp; $_=~s/\s+/\t/g;print $_,"\n";' | memset.run - $HMHOME/data/hg19_chr22.fa.gz 
}


#FA=~/hmdata/ucsc/hg19/hg19.fa.gz;
#mycat $HMHOME/data/hg19_ensGene_coding.bed.gz | bed.intron - > intron.bed





cmd=`cat $HMHOME/lib/MEMset/score3.pl`
cmd=${cmd/INPUT1/$HMHOME/lib/MEMset/me2x5}
cmd=${cmd/INPUT2/$HMHOME/lib/MEMset/splicemodels}


