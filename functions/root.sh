#!/bin/bash

check(){ ## obtained from bamtools test
	if diff $1 $2; then
		echo ok
		return 1
	else
		echo fail
		return 0
	fi
}

make_tempdir(){
	mktemp -d 2>/dev/null || mktemp -d -t 0;
}

split_by_chrom(){
	awk -v OFS="\t" -v O=$2 '{
		fout=O"/"$1;
		print $0 >> fout;
	}' $1
	echo `ls $2/*`;	## return a list of splited files
}

modify_score(){
usage="$FUNCNAME <bed6> <method>
	<method>: count phred
"
	if [ $# -ne 2 ]; then echo "$usage"; return; fi
	awk -v OFS="\t" -v ME=$2 '{
		if(ME=="count"){
			$5=1;
		}else if(ME=="phred"){
			if ( $5 == 0){
				$5 = 0.1;
			}else{
				$5 = 1- exp( - $5/10 * log(10));
			}
		}
		print $0;
	}' $1;
}
sum_score(){
	tmpd=`make_tempdir`
	for f in `split_by_chrom $1 $tmpd`;do
		sort -k2,3n $f | groupBy -g 1,2,3,6 -c 5 -o sum \
		| awk -v OFS="\t" '{ print $1,$2,$3,".",$5,$4;}'
	done
}


###########################################################
# test 
###########################################################
test(){
echo \
"chr1	95	100	 a1	1	+
chr1	200	205	 a1	1	+
chr1	95	100	 a2	2	-
chr2	200	205	 a2	2	-" > a.bed

echo "test .. modify_score a.bed count"
echo \
"chr1	95	100	a1	1	+
chr1	200	205	a1	1	+
chr1	95	100	a2	1	-
chr2	200	205	a2	1	-" > exp
modify_score a.bed count > obs
check obs exp

echo "test .. modify_score a.bed phred"
echo \
"chr1	95	100	a1	0.205672	+
chr1	200	205	a1	0.205672	+
chr1	95	100	a2	0.369043	-
chr2	200	205	a2	0.369043	-" > exp 
modify_score a.bed phred > obs
check obs exp

echo "test .. split_by_chrom a.bed out"
echo \
"out/chr1 out/chr2" > exp
mkdir -p out
split_by_chrom a.bed out > obs
check obs exp
rm obs exp 
rm -rf out

echo "test .. sum_score "
echo \
"chr1	95	100	a1	1	+
chr1	95	100	a1	1	+
chr1	95	100	a2	1	-
chr2	95	100	a2	1	-" > inp
echo \
"chr1	95	100	.	2	+
chr1	95	100	.	1	-
chr2	95	100	.	1	-" > exp

sum_score inp > obs
check obs exp
rm obs exp inp
rm a.bed
}


