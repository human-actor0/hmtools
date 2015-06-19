#!/bin/bash
. $HMHOME/src/root.sh

sort_bed(){
	cat $1 | sortBed -i stdin
	#sort -k1,1 -k2,2n $1 
}

get_chromsize(){
	if [ ! -f $1.bai ];then
		samtools index $1;
	fi
	samtools idxstats $1 | awk -v OFS="\t" '$1 != "*" && $3 > 0 { print $1,$2;}'
}

split_bam(){
	mkdir -p $2;
	for chrom in `get_chromsize $1 | cut -f1`;do
		echo " spliting $1 to $2/$chrom.bam .. " >&2;
		samtools view -b $1 $chrom > $2/$chrom
	done
	echo `ls $2/*`;
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
	awk -v OFS="\t" '{ $1=$1","$6;print $0;}' $1  \
	| sort_bed - | groupBy -g 1,2,3 -c 5 -o sum \
	| awk -v OFS="\t" '{ split($1,a,","); print a[1],$2,$3,".",$4,a[2];}'
}

bed_count(){
usage="
usage: $FUNCNAME <target> <read> [ <zero> [<strand>]]
output: target + sum of read scores
use modify_score first this to change scoring method
"
	opt_zero=${3:-0};
	opt_strand=${4:-""};
	if [ $# -lt 2 ];then echo "$usage"; return; fi

	tmpd=`make_tempdir`;
	mycat $1 > $tmpd/a
	mycat $2 > $tmpd/b
	n=`head -n 1 $tmpd/a | awk '{print NF;}'`
	intersectBed -a $tmpd/a -b $tmpd/b -wa -wb $opt_strand \
	| awk -v n=$n -v OFS="\t" '{ 
		for(i=2; i<=n;i++){
			$1=$1"@"$(i);
		} print $1,$(n+5);
	}' | groupBy -g 1 -c 2 -o sum | tr "@" "\t"

	if [ $opt_zero -ne 0 ]; then
		intersectBed -a $tmpd/a -b $tmpd/b -wa -v $opt_strand \
		| awk -v OFS="\t" '{ print $0,0;}'
	fi
	rm -rf $tmpd
}
_test_bed_count(){
echo \
'chr	1	100	n1
chr	50	200	n2'> a
echo \
'chr	1	10	r1	1	+
chr	40	50	r2	2	+
chr	50	200	r3	3	+' > b

bed_count a b  > obs
echo \
'chr	1	100	n1	6
chr	50	200	n2	3' > exp
check obs exp 

rm a b exp obs
}
#_test_bed_count

intersectBed_sum(){
usage="
	$FUNCNAME <target> <read> [<intersectBed options>]
	<intersectBed options>: -s 
"
	OPT="";
	if [ $# -lt 2 ];then echo "$usage";return; fi
	if [ $# -gt 2 ];then OPT=${@:3}; fi
	tmpd=`make_tempdir`
	mycat $1 | cut -f1-6 > $tmpd/a
	mycat $2 | cut -f1-6 > $tmpd/b

	intersectBed -a $tmpd/a -b $tmpd/b -wa -wb $OPT \
	| awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$11;}' \
	| groupBy -g 1,2,3,4,6 -c 7 -o sum \
	| awk -v OFS="\t" '{ print $1,$2,$3,$4,$6,$5;}'  
	## zero counts
	intersectBed -a $tmpd/a -b $tmpd/b -v $OPT \
	| awk -v OFS="\t" '{ print $1,$2,$3,$4,0,$6;}'
	rm -rf $tmpd
}
_test_intersectBed_sum(){
echo \
"chr1	100	200	n1	0	+
chr1	200	300	n2	0	-" > t
echo \
"chr1	100	101	r1	1	+
chr1	199	201	r3	5	-
chr1	200	201	r2	10	+" > a
echo \
"chr1	100	200	n1	1	+
chr1	200	300	n2	5	-" > exp

intersectBed_sum t a -s > obs
echo "test .. intersectBed_sum"
check obs exp
rm t a obs exp
}

bed12_to_lastexon(){
## not tested
	awk -v OFS="\t" '{ split($11,sizes,",");split($12,starts,",");
	    if($6=="+"){ i=$10;}else{ i=1;}
	    s=$2+starts[i]; e=s+sizes[i];
	    print $1,s,e,$4,i,$6;
	}' $1
}
bed12_to_exon(){
## not tested
	awk -v OFS="\t" '{ split($11,sizes,",");split($12,starts,",");
		for(i=1; i<=$10; i++){
	    		s=$2+starts[i]; e=s+sizes[i];
			print $1,s,e,$4,$5,$6;
		}
	}' $1 
}

bed12_to_intron(){
	##  [es   ]s----e[    ee]
	awk -v OFS="\t" '$10 > 1{
		## take introns
		split($11,sizes,",");
		split($12,starts,",");
		for(i=1;i< $10;i++){
			## intron
			s = $2 + starts[i]+sizes[i];	
			e = $2 + starts[i+1];
			#ls = $2 + starts[i-1]; le = ls + sizes[i-1]; rs = $2 + starts[i]; re = rs + sizes[i];
			print $1,s,e,$4,$5,$6;
		}	
	}' $1 
}
bed_flat(){
    ## input bed features
    ## [     ]-----------[        ]
    ##    [      ]-----[    ]
    ## output exon  fragments (: open [: closed intervals
    ## [ )[  ](  ]     [)[  ](    ]  
    #awk -v OFS="\t" '{ print $1,$2,$3,$4 "," $2 "," $3,$5,$6;}' | mergeBed -nms -s -scores collapse \
	sort_bed $1 \
	| mergeBed -i stdin -d -1 -c 2,3 -o distinct,distinct \
	| perl -ne ' chomp; my @a=split /\t/,$_; 
		my %pm=();
		foreach my $e (split/,/,$a[4]){ $pm{$e}=1;}
		foreach my $e (split/,/,$a[5]){ $pm{$e}=1;}
		if(scalar keys %pm ==1){
			print $a[0],"\t",$a[1],"\t",$a[2],"\n";
			next;
		}
		my @p=sort {$a<=>$b} keys %pm; 
		for(my $i=0; $i< $#p; $i++){
		    my $pi = $p[$i]; my $pj = $p[$i+1];
		    if($pj > $pi){
			print $a[0],"\t",$pi,"\t",$pj,"\n";
		    }
		} 
	' | sort -u 
}

_test_flat(){
echo \
'chr	11	13	.	0	+
chr	12	14	.	0	+
chr	1	4	.	0	+
chr	2	3	.	0	+' | flatten_bed - >obs
echo \
'chr	1	2		0	+
chr	2	3		0	+
chr	3	4		0	+
chr	11	12		0	+
chr	12	13		0	+
chr	13	14		0	+' > exp
check exp obs
rm exp obs
}
#_test_flat

bed12_to_junction(){
	bed12_to_intron $1 \
	| awk -v OFS="\t" '{ $2=$2-1; $3=$3+1;print $0;}'
}

merge_by_gene(){
## not tested
        perl -ne 'chomp;my @a=split/\t/,$_;
                $a[0]=$a[0]."@".$a[3];  ## avoid merging different genes
                $a[4]=0; 
                print join("\t",@a),"\n";' \
        | sort_bed - \
        | mergeBed -i stdin -s -c 4,5,6 -o distinct,count,distinct \
        | awk -v OFS="\t" '{ split($1,a,"@");$1=a[1];print $0;}'
}
bed12_to_3utr(){
	#$chr,,$start,,$end,,$name,,$score,,$strand,,$thickStart,,$thickEnd,,$itemRgb,,$blockCount,;
	awk -v OFS="\t" '{
		split($11,l,",");
		split($12,s,",");
		coding=1;
		if( $7 == $8) coding=0; ## noncoding
		if($6=="+"){
			start=$2+s[$10]; end=start+l[$10];
			if(coding) start=$8;
		}else{
			start=$2; end=start+l[$10]; 
			if(coding) end=$7;
		}
		if(coding && end > start) ## remove coding end points
			print $1,start,end,$4,$5,$6;
	}' $1;
}

ucsc_to_bed12(){
	cmd='
	    chomp;$_=~s/\r//g;
	    my @aa = split /\s/,$_;
	    my ($bin,$name,$chr,$strand,$start,$end,$thickStart,$thickEnd,$blockCount,$blockStarts,$blockEnds,$id,$name2) = split /\s/, $_;
	    my $itemRgb = "255,0,0";
	    my $score = 0;

	    if(defined $name2){
		$name = $name."|".$name2;
	    }
	    print $chr,"\t",$start,"\t",$end,"\t",$name,"\t",$score,"\t",$strand,"\t",$thickStart,"\t",$thickEnd,"\t",$itemRgb,"\t",$blockCount,"\t";
	    my @ss = split /,/,$blockStarts;
	    my @ee = split /,/,$blockEnds;
	    for(my $i=0;$i<$blockCount;$i++){
		my $length = $ee[$i]-$ss[$i];
		print $length,",";
	    }
	    print "\t";
	    for(my $i=0;$i<$blockCount;$i++){
		my $relstart = $ss[$i]-$start;
		print $relstart,",";
	    }
	    print "\n";
	'
	mycat $1 | perl -ne "$cmd";

}





get3utr(){
        #cat $1 | gtf_to_bed12.sh |bed12_to_lastexon.sh | perl -ne 'chomp;my @a=split/\t/,$_;
        cat $1 | bed12_to_lastexon.sh | mergeByGene
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

echo "test .. ucsc_to_bed12";
echo \
"585	ENST00000456328	chr1	+	11868	14409	14409	14409	3	11868,12612,13220,	12227,12721,14409,	0	ENSG00000223972	none	none	-1,-1,-1,
585	ENST00000515242	chr1	+	11871	14412	14412	14412	3	11871,12612,13224,	12227,12721,14412,	0	ENSG00000223972	none	none	-1,-1,-1," \
| ucsc_to_bed12 - > obs

echo \
"chr1	11868	14409	ENST00000456328|ENSG00000223972	0	+	14409	14409	255,0,0	3	359,109,1189,	0,744,1352,
chr1	11871	14412	ENST00000515242|ENSG00000223972	0	+	14412	14412	255,0,0	3	356,109,1188,	0,741,1353," > exp
check obs exp
rm obs exp
}


