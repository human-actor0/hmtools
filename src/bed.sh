#!/bin/bash
. $HMHOME/src/root.sh
. $HMHOME/src/stat.sh

bed.mhits(){

usage="
FUNCT: redistribute multi-hits
USAGE: $FUNCNAME <target> <multi_reads> [-s|-S]
"; if [ $# -lt 2 ];then echo $usage; return; fi
local MAXITER=100;
local DELTA=0.01;

	intersectBed -a $1 -b $2 -wa -wb ${3:-""} \
	| perl -e 'use strict; 

		my $MAXITER='$MAXITER';
		my $DELTA='$DELTA';

		my %E=(); # gene expression frequency
		my %A=(); # read to gene fraction
		my %R=(); # read to gene assignment
		my %G=(); # gene to read assignment
		my %L=(); # gene lengths
		my $s=0;
		while(<STDIN>){chomp;
			my ($chr,$start,$end,$gene,$score,$strand, $chr2,$start1,$end1,$read,$score1,$strand1)=split /\t/,$_;
			$L{$gene}=$end-$start;
			$E{$gene} += $score1/$L{$gene};	
			$R{$read}{$gene} = $score1;
			$G{$gene}{$read} = 1;
			$s += $E{$gene};
		}
		foreach my $k (keys %E){ $E{$k} /= $s; }

		foreach my $k (keys %E){ print $k," ",$E{$k},"\n"; }

		foreach my $iter ( 1..$N){
			## update A
			foreach my $r (keys %R){
				my $s=0;
				foreach my $g (keys %{$R{$r}}){
					$A{$r}{$g} += $E{$g};
					$s += $E{$g};
				}
				foreach my $g (keys %{$R{$r}}){
					$A{$r}{$g} /= $s;
				}
			}

			## update E
			$s=0;
			foreach my $g (keys %E){
				my $s2=0;
				foreach my $r (keys %{$G{$g}}){
					$s2 += $A{$r}{$g};	
				}
				$E{$g}=$s2/$L{$g};
				$s += $E{$g};
			}
			foreach my $k (keys %E){ $E{$k} /= $s; }

			#foreach my $k1 (keys %A){ 
			#foreach my $k2 (keys %{$A{$k1}}){ print $k1," ",$k2," ",$A{$k1}{$k2},"\n"; }}
			print "ITER--$iter\n";
			foreach my $k (keys %E){ print $k," ",$E{$k},"\n"; }
		}

		
	'
}
bed.mhits.test(){
echo \
"c	0	100	t1	0	+
c	100	200	t2	0	+
c	300	400	t3	0	+" > genes

echo \
"c	10	20	r1	1	+
c	100	120	r1	1	+
c	1000	1020	r1	1	+
c	110	130	r2	1	+
c	320	340	r2	1	+
c	50	70	r2	1	+
c	50	70	r3	1	+
c	350	370	r3	1	+ " > reads

bed.multihits genes reads
}

bed_nf(){
	head -n 1 $1 | awk '{print NF;}'
}

bed_join(){
usage="
usage: $FUNCNAME <bed> [<bed>..]
function: sum and join counts per entry 
"; if [ $# -lt 1 ];then echo "$usage">&2 ;return; fi;
	perl -e '
		use strict;
		my @files=@ARGV;
		my $i=0;
		my $nf=scalar @files;
		my %data=();
		foreach my $f (@files){
			open my $in, "<", $f;
			while(<$in>){ chomp; my @a=split/\t/,$_;
				my $k= join(";",(@a[0..3],$a[5]));
				my $v= $a[4];
				if( !defined $data{$k} ){
					my @na=( 0 )x $nf;
					$data{$k} = \@na;
				}
				$data{$k}->[$i]=$v;
			}
			close($in);
			$i++;
		}
		sub sum{
			my $res=0;
			foreach my $i (@_){ $res += $i; }
			return $res;
		}
		foreach my $k (keys %data){
			my ($c, $s, $e, $n, $st) = split /;/,$k;
			print "$c\t$s\t$e\t$n\t",sum(@{$data{$k}}),"\t$st\t", join("\t",@{$data{$k}}),"\n";
		}
	' $@ 
}
test__bed_join(){
echo "c	1	2	n1	1	+
c	2	3	n2	100	+" > tmpa
echo "c	1	2	n1	10	+
c	3	4	n3	1000	+" > tmpb
echo \
"c	3	4	n3	1000	+	0	1000
c	2	3	n2	100	+	100	0
c	1	2	n1	11	+	1	10" > exp
bed_join tmpa tmpb  > obs
check exp obs
rm -f tmpa tmpb exp obs
}

bede_join(){
	echo "$@" | perl -e '
		use strict;
		my @files=split /\s+/,<STDIN>;
		my $i=0;
		my $n=scalar @files;
		my %data=();
		my %nc = ();
		foreach my $f (@files){
			open my $in, "<", $f;
			while(<$in>){ chomp; my @a=split/\t/,$_;
				my $k= join("@",@a[0..5]);
				my $v= join("@",@a[6..$#a]);
				$data{$k}{$i}=$v;
				if(!defined $nc{$i}){
					$nc{$i} = $#a - 5;
				}
			}
			close($in);
			$i++;
		}
		foreach my $k (keys %data){
			print $k,"\t";
			for(my $i=0; $i < $n; $i++){
				if(defined $data{$k}{$i}){
					print "\t",$data{$k}{$i};
				}else{
					print "\t",join("@",( 0 ) x $nc{$i});	
				}
			}
			print "\n";
		}
	' | tr "@" "\t"
}

bed_flank(){
	local L=$2; local R=$3; local S=${4:-""};
	awk -v OFS="\t" -v s=${2:-""} -v L=$L -v R=$R -v S=$S '{ 
		if(S == "s" && $6 == "-"){
			$2=$2-R; $3=$3+L;
		}else{ $2=$2-L; $3=$3+R; } print $0;
	}' $1;
}

bed_3p(){
 	awk -v OFS="\t" '{ if($6=="-"){$3=$2+1;} $2=$3-1; print $0; }' $1
}

bed_5p(){
 	awk -v OFS="\t" '{if($6=="-"){$2=$3-1;}$3=$2+1; print $0; }' $1
}
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

bychrom(){
usage="$FUNCNAME [samtools_ops] <bam> <functions>"
	local tmpa="";
	local bam="";
	for e in $@;do
		shift; tmpa+=" $e"
		if [ ${e##*\.} = "bam" ];then bam=$e; break; fi;
	done
	for chrom in `get_chromsize $bam | cut -f 1`;do
		if [ -n "$DEBUG" ];then echo " $chrom .. ">&2; fi
		samtools view $tmpa $chrom | $@
	done;
}

split_bam(){
	mkdir -p $2;
	for chrom in `get_chromsize $1 | cut -f1`;do
		echo " spliting $1 to $2/$chrom.bam .. " >&2;
		samtools view -b $1 $chrom > $2/$chrom
	done
	echo `ls $2/*`;
}

bed_splitByChrom(){
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
			if( $5 > 0){
				$5 = 1- exp( - $5/10 * log(10));
			}
		}
		print $0;
	}' $1;
}
bed_add(){
	perl -e 'use strict; my %res=();
	foreach my $f (@ARGV){
		my $fh;
		open($fh, $f) or die "$!";
		while(<$fh>){chomp; my ($k,$v)=split/\t/,$_;
			$res{$k} += $v;
		}
		close($fh);	
	}
	foreach my $k (keys %res){ print $k,"\t",$res{$k},"\n";}
	' $@;
}
bed_sum(){
	perl -e 'use strict; my %res=();
	foreach my $f (@ARGV){
		my $fh; open($fh, $f) or die "$!";
		while(<$fh>){ chomp; my @a=split/\t/,$_;
			$res{ join("@",(@a[0..3],$a[5])) } += $a[4]; 
		}
		close($fh);
	}
	foreach my $k (keys %res){ 
		my ($c,$s,$e,$n,$st) = split /@/,$k;
		print $c,"\t",$s,"\t",$e,"\t",$n,"\t",$res{$k},"\t",$st,"\n";
	}
	' $@; 
}
test__bed_sum(){
echo \
"chr1	1	2	a	1	+
chr1	1	2	a	2	-
chr1	1	2	a	3	+" | bed_sum - > obs

echo \
"chr1	1	2	a	2	-
chr1	1	2	a	4	+" > exp
check exp obs
rm -f obs exp
}
#test__bed_sum

sum_score(){
	awk -v OFS="\t" '{ $1=$1","$6;print $0;}' $1  \
	| sort_bed - | groupBy -g 1,2,3 -c 5 -o sum \
	| awk -v OFS="\t" '{ split($1,a,","); print a[1],$2,$3,".",$4,a[2];}'
}

bed_count(){
usage="
usage: $FUNCNAME <target> <read> [ <zero> [<strand>]]
output: target + sum of read scores
use modify_score to change scoring 
"
	opt_zero=${3:-0};
	opt_strand=${4:-""};
	if [ $# -lt 2 ];then echo "$usage"; return; fi

	local tmpd=`make_tempdir`;
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
	local tmpd=`make_tempdir`
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
usage="$FUNCNAME [options] <bed6>
     input:
     [     ]
        [      ]
     output:
     [ ][  ][  ] 
"; if [ $# -lt 1 ]; then echo "$usage";exit; fi

	mergeBed -i stdin -d -1 -c 2,3 -o distinct,distinct \
	| perl -ne ' chomp; my @a=split /\t/,$_; 
		my %pm=();
		foreach my $s (split/,/,$a[3]){ $pm{$s}=1;}
		foreach my $e (split/,/,$a[4]){ $pm{$e}=1;}
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
	' 
}

test__bed_flat(){
echo \
'chr	11	13	.	1	+
chr	12	14	.	2	-
chr	12	14	.	3	-
chr	1	4	.	4	+
chr	2	3	.	5	+' > inp

cat inp | awk -v OFS="\t" '{ $1=$1"@"$6; $6="";$0=$0;}1'   \
| sort -k1,1 -k2,3n \
| bed_flat -  > obs
echo \
"chr@+	1	2
chr@+	2	3
chr@+	3	4
chr@+	11	13
chr@-	12	14" > exp
check obs exp
rm obs exp
cat inp \
| sort -k1,1 -k2,3n \
| bed_flat -  > obs
echo \
"chr	1	2
chr	2	3
chr	3	4
chr	11	12
chr	12	13
chr	13	14" > exp
check obs exp
rm obs exp inp
}
#test__bed_flat


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

bed_igx(){
usage=" usage: $FUNCNAME <group.bed> <feature.bed6+>
        bed6+: bed6 + counts (tab-delimited)
"
if [ $# -ne 2 ];then echo "$usage"; return; fi
        local tmpd=`make_tempdir`;
        mycat $1 > $tmpd/a; mycat $2 > $tmpd/b
        intersectBed -b $tmpd/a -a $tmpd/b -wa -wb -s \
       	| perl -ne 'chomp; my @a=split/\t/,$_;
                my $id=join("@",@a[0..5]);
                my $j=0; for(my $i=6; $i <= $#a; $i++){
                        if($a[$i] eq $a[0]){  $j=$i; last;}
                }
                my $group=join("@",@a[$j..$#a]);
                my $xs=join("\t",@a[6..($j-1)]);
                print $id,"\t",$group,"\t",$xs,"\n";
        '
        rm -rf $tmpd
}


bed12_to_fexons(){
## consider all exons belonging to a gene
## bed12 input file should contain gene_id at 4th column
        #mycat /main/hmtools/data/hg19_ensGene_coding.bed.gz \
        bed12_to_exon $1 \
        | awk -v OFS="\t" '{ print $1"@"$4"@"$6, $2,$3;}' \
        | sort -k1,1 -k2,3n -u | bed_flat -\
        | tr "@" "\t" |  awk -v OFS="\t" '{ print $1,$4,$5,$2,0,$3;}' \
        | sort -k1,1 -k2,3n 
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


