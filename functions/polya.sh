#!/bin/bash  
. $HM/functions/root.sh # import utilities 
. $HM/functions/stat.sh # import test_lineartrend

SEQ=$HM/standalones/seq.sh
FILTER=$HM/standalones/pa_filter_nb.sh
FILTER_M=$HM/standalones/nb.model
HG19FA=/mnt/db/Ucsc/hg19/
CLUSTER=$HM/standalones/cluster_1d_kmeans.sh

## take 3' end point, then switch strand, and then sum scores

point(){
	awk -v OFS="\t" '{ 
		if($6 == "+"){ $3=$2+1; $6="-"; }else{ $2=$3-1; $6="+"; } print $0; 
	}' $1 | sum_score -
}

filter(){ 
	usage="$FUNCNAME <bed> <fasta>"
	if [ $# -ne 2 ]; then echo "$usage"; return; fi

	eval "$SEQ -s -l 39 -r 30 $1 $2" \
	| perl -ne 'chomp;my @a=split/\t/,$_;
		my $s=pop @a; $s=~ s/,//g;
		print join("@",@a),"\t",$s,"\n"; '\
	| eval "$FILTER predict - $FILTER_M" \
	| perl -e 'use strict; my $offset=40; my %S=();
		my $total_sum=0; my $passed_sum=0;
		my $total_pos=0; my $passed_pos=0;
		while(<>){ chomp;
			chomp;my ($bed,$seq,$pos,$score) = split/\t/,$_;
			next if $score eq "";
			my @a=split /@/,$bed;
			$total_pos ++;
			$total_sum += $a[4];
			if($#a >=5 && $a[5] eq "-"){
				my @b=split/,/,$score;
				$score = join(",",reverse @b);
			}
			if($score > 0.5){
				print join( "\t",@a),"\t$score\n";
				$passed_pos ++;
				$passed_sum += $a[4];
			}
		}
	'
}

cluster(){
	MIND=$2;
	tmpd=`make_tempdir`;
	mergeBed -i $1 -s -c 6,2,5 -o distinct,collapse,collapse -d $MIND > $tmpd/a.bed
        cat $tmpd/a.bed | awk -v OFS="\t" -v M=$MIND '
        $3-$2 < M {
                L=split($5,a,","); split($6,b,",");
                n=0; s=0; y=0;
                for(i=1; i<=L; i++){
                        s += a[i] * b[i];
                        n += b[i];
                }
                print $1,$2,$3,int(s/n),int(n),$4;
        }'  > $tmpd/b.bed

        cat $tmpd/a.bed | awk -v M=$MIND '$3-$2 >= M' \
         | eval $CLUSTER -d $MIND -c 5,6 - \
         | awk -v OFS="\t" '{
                L=split($7,a,","); split($8,b,","); split($9,c,","); split($10,d,",");
                for(i=1;i<=L;i++){
                        print $1,a[i],b[i],int(c[i]),int(d[i]),$4;
                }
        }' > $tmpd/c.bed
	cat $tmpd/c.bed $tmpd/b.bed
}

_precompare(){
	sort -k1,1 -k2,3n $1 \
	| intersectBed -a stdin -b $2 -wa -wb -s \
	| groupBy -g 1,2,3,4,5,6 -c 8,11 -o collapse,collapse \
	| intersectBed -a stdin -b $3 -wa -wb -s \
	| groupBy -g 1,2,3,4,5,6,7,8 -c 10,13 -o collapse,collapse \
	| awk -v OFS="\t" '{ print $1"@"$2"@"$3"@"$4"@"$5"@"$6,$7,$8,$9,$10;}'
}
compare(){
	_precompare $1 $2 $3 | test_lineartrend - \
	| tr "@" "\t" \
	| awk -v OFS="\t" '{ if($6 == "-"){ $(NF-1) = - $(NF-1);} print $0;}'
}


###########################################################
# test 
###########################################################
test(){
echo \
"chr1	95	100	 a1	1	+
chr1	95	100	 a1	1	+
chr1	95	100	 a2	2	-
chr2	200	205	 a2	2	-" > inp

echo "test .. point inp"
echo \
"chr1	95	96	.	2	-
chr1	99	100	.	2	+
chr2	204	205	.	2	+" > exp
point inp > obs
check obs exp
rm obs exp inp

echo "test .. filter inp hg19"
echo \
"chr7	38765151	38765152	polya	10	-
chr7	38765161	38765162	polya	10	-
chr7	38765151	38765152	polya	10	+
chr7	38765161	38765162	polya	10	+" > inp
echo \
"chr7	38765151	38765152	polya	10	-	0.999999341801724
chr7	38765161	38765162	polya	10	-	0.999999989324471" > exp
filter inp $HG19FA > obs
check obs exp
rm obs exp inp

echo "test .. compare"
echo \
"chr1	1	200	g1	0	-" > c.bed
echo \
"chr1	5	6	.	1	-
chr1	7	8	.	2	-
chr1	99	100	.	3	-" > a.bed
echo \
"chr1	5	6	.	100	-
chr1	7	8	.	20	-
chr1	99	100	.	1	-" > b.bed

echo \
"chr1	1	200	g1	0	-	5,7,99	1,2,3	5,7,99	100,20,1	0.603335	1.26652022203e-11" > exp

compare c.bed a.bed b.bed > obs
check exp obs
rm a.bed b.bed c.bed exp obs

echo "test .. cluster"
echo \
"chr1	1	2	.	10	-
chr1	2	3	.	2	-
chr1	3	4	.	10	-
chr1	5	6	.	3	-
chr1	6	7	.	10	-
chr1	7	8	.	4	-" \
| cluster - 2 >obs

echo \
"chr1	1	3	1	12	-
chr1	3	4	3	10	-
chr1	5	8	6	17	-" > exp
check exp obs
rm exp obs

}

typeset -F
