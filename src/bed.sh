#!/bin/bash
. $HMHOME/src/root.sh


get_chromsize(){
	samtools idxstats $1 | awk -v OFS="\t" '$1 != "*" && $3 > 0 { print $1,$2;}'
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
	#tmpd=`make_tempdir`
	#for f in `split_by_chrom $1 $tmpd`;do
		awk -v OFS="\t" '{ $1=$1","$6;print $0;}' $1  \
		|  sortBed -i stdin | groupBy -g 1,2,3 -c 5 -o sum \
		| awk -v OFS="\t" '{ split($1,a,","); print a[1],$2,$3,".",$4,a[2];}' \
		| sortBed -i stdin
	#done
	#rm -rf $tmpd
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
			print $1,s,e,$4,0,$6;
		}
	}' $1 | sortBed -i stdin | sort -u 
}

merge_by_gene(){
## not tested
        perl -ne 'chomp;my @a=split/\t/,$_;
                $a[0]=$a[0]."@".$a[3];  ## avoid merging different genes
                $a[4]=0; 
                print join("\t",@a),"\n";' \
        | sortBed -i stdin \
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


