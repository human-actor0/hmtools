#/bin/bash  
. $HMHOME/src/root.sh # import utilities 
. $HMHOME/src/bed.sh #import utilities 
. $HMHOME/src/stat.sh # import test_lineartrend test_fisherexact

SEQ=$HMHOME/bin/bed_seq.sh
FILTER=$HMHOME/src/pa_filter_nb.sh
FILTER_M=$HMHOME/src/nb.model
HG19FA=/hmdata/ucsc/hg19/chromosome/
CLUSTER=$HMHOME/src/cluster_1d_kmeans.sh
BW=$HMHOME/bin/bedGraphToBigWig

pa_test_pcpa_vs_pa(){
usage="
$FUNCNAME <gene.bed> <3utr.bed> <trt.bed> <ctr.bed> [option]
"
	local tmpd=`mymktempd`;
	local opts=${5:-""};
	intersectBed -a $1 -b $3 -wa -wb $opts \
	| awk -v OFS=";" '{ print $1,$2,$3,$4,$5,$6"\t"$11;}' \
	| sum - > $tmpd/trt 

	intersectBed -a $1 -b $4 -wa -wb $opts \
	| awk -v OFS=";" '{ print $1,$2,$3,$4,$5,$6"\t"$11;}' \
	| sum - > $tmpd/ctr 


}

pa_de(){
usage="
function: determine Positive,Negative, Unchange, Ignore
usage: $FUNCNAME <fisherTestOut> <fdr> [type]
 [type]: default oddratio, FDR at last -2 and last columns
"
	local FDR=$2;
	local tmpd=`mymktempd`;
	cat $1 > $tmpd/a
	local Type=${3:-"fisher"};
	if [ $Type = "3utr" ];then
		cat $tmpd/a | awk -v OFS="\t" -v FDR=$FDR -v MINS=$mins '{
			s=$8+$9+$10+$11; fdr=$(NF); OR=$(NF-2);
			d="I";#low sample
			if( fdr <= FDR){ if(OR > 1){ d="P"; }else{ d="N"; } }
			else if( s >= MINS){ d="U";}
			print $7,d;
		}' | perl -e 'use strict;
			my %res=();
			while(<STDIN>){chomp; my ($id,$e) = split /\t/,$_;
				$res{$id}{$e} = 1;
			}
			foreach my $k (keys %res){
				if(defined $res{$k}{P} || defined $res{$k}{N}){
					print $k,"\t","A","\n";
				}elsif(defined $res{$k}{U}){
					print $k,"\t","U","\n";
				}else{
					print $k,"\t","I","\n";
				}
			}
		' | tr "," "\t"
	else
		cat $tmpd/a | awk -v OFS="\t" -v FDR=$FDR -v MINS=$mins '{
			s=$8+$9+$10+$11; fdr=$(NF); OR=$(NF-2);
			d="I";#low sample
			if( fdr <= FDR){ if(OR > 1){ d="P"; }else{ d="N"; } }
			else if( s >= MINS){ d="U";}
			print $1,$2,$3,$4,$5,$6,d;
		}'
	fi
	rm -rf $tmpd
}

pa_de2(){
	local tmpd=`mymktempd`;
	pa_de $1 $3 ${4:-"fisher"} | awk -v OFS=";" '{print $1,$2,$3,$4,$5,$6"\t"$7;}' | sort -k 1,1 > $tmpd/a
	pa_de $2 $3 ${4:-"fisher"} | awk -v OFS=";" '{print $1,$2,$3,$4,$5,$6"\t"$7;}' | sort -k 1,1 > $tmpd/b
	join -a 1 -a 2 -e "I" -o 0,1.2,2.2 $tmpd/a $tmpd/b | tr " ;" "\t" 
	rm -rf $tmpd
}
pa_de_merge(){
	local files=`echo $@ | quote`;
	cmd='files=c('"$files"');
		d=NULL;
		for( f in files ){
			tt=read.table(f,header=F,stringsAsFactors=F);
			colnames(tt)=c("chrom","start","end","name","score","strand",f);
			if(is.null(d)){ d=tt;
			}else{ d=merge(d,tt,by=c(1:6),all=T); }
		}
		#colnames(d)=gsub("tmpd\\/","",colnames(d));
		d[ is.na(d) ] = "I";
		write.table(d,file="stdout",col.names=T,row.names=F,quote=F,sep="\t");
	'
	#echo "$cmd";
	run_R "$cmd" 
}

pa_prep_test_relativefreq(){
usage="
        $FUNCNAME <target> <trt.bed> <ctr.bed> [-s|-S]
"
if [ $# -lt 3 ];then echo "$usage"; return; fi
        local opts=${4:-""};
        bed_join $2 $3 \
        | intersectBed -a $1 -b stdin -wa -wb $opts \
        | perl -ne 'chomp; my @a=split /\t/,$_;
                $a[10]=0; ## this allows cross comparisons
                print join(";",@a[6..11]),"\t",
                        join(",",@a[0..5]),"\t",
                        join("\t",@a[12..$#a]),"\n";' \
        | igx_to_igxnx - \
        | awk -v OFS="\t" '{ print $1"@"$2,$3,$4,$5,$6;}'
}

pa_test_relativefreq(){
usage="
	$FUNCNAME <target> <trt.bed> <ctr.bed> [-s|-S] 
"
if [ $# -lt 3 ];then echo "$usage"; return; fi
	pa_prep_test_relativefreq $@ \
	| fisher_test - | padjust - -1 \
	| tr "@" "\t" | tr ";" "\t" 
}

pa_test_pcpa_to_pa(){
usage="
function: test the ratio of PCPA to PA 
usage: $FUNCNAME <gene.bed> <3utr.bed> <trt.bed> <ctr.bed> [options]
 <3utr.bed> : filter for PCPA
 [options] : [-s|-S]
"
	local tmpd=`mymktempd`;
	local opts=${5:-""};
	#local tmpd=tmpd;mkdir -p $tmpd;	
	intersectBed -a $1 -b $3 -wa -wb $opts \
	| awk -v OFS=";" '{ print $1,$2,$3,$4,$5,$6"\t"$11;}'  | sum - \
	| sort -k1,1 > $tmpd/trt
	
	intersectBed -a $3 -b $2 -v $opts \
	| intersectBed -a $1 -b stdin -wa -wb $opts \
	| awk -v OFS=";" '{ print $1,$2,$3,$4,$5,$6"\t"$11;}'  | sum - \
	| sort -k1,1 > $tmpd/trt-3utr
		

	intersectBed -a $1 -b $4 -wa -wb $opts \
	| awk -v OFS=";" '{ print $1,$2,$3,$4,$5,$6"\t"$11;}'  | sum - \
	| sort -k1,1 > $tmpd/ctr

	intersectBed -a $4 -b $2 -v $opts\
	| intersectBed -a $1 -b stdin -wa -wb $opts \
	| awk -v OFS=";" '{ print $1,$2,$3,$4,$5,$6"\t"$11;}'  | sum - \
	| sort -k1,1 > $tmpd/ctr-3utr

	join -a 1 -a 2 -e 0 -o 0,1.2,2.2 $tmpd/trt-3utr $tmpd/trt \
	| join -a 1 -a 2 -e 0 -o 0,1.2,1.3,2.2 - $tmpd/ctr-3utr \
	| join -a 1 -a 2 -e 0 -o 0,1.2,1.3,1.4,2.2 - $tmpd/ctr 	\
	| awk -v OFS="\t" '{ $3=$3-$2; $5=$5-$4;} 1' \
	| fisher_test - | padjust - -1
	rm -rf $tmpd;	
}

pa_test_diff_fisher(){
        bed_join $1 $2 \
        | awk -v OFS=";" '{ print $1,$2,$3,$4,$5,$6"\t"1"\t"$7"\t"$8;}' \
        | igx_to_igxnx - \
        | cut -f1,3-6 \
        | fisher_test2 - \
        | padjust - -1
}


pa_test_lineartrend(){
usage="$FUNCNAME <target> <trt.bed> <ctr.bed> [-S|-s]
"
if [ $# -lt 3 ];then echo "$usage"; return; fi
	local tmpd=`mymktempd`; #local tmpd=tmpd; mkdir -p $tmpd
	local opts=${4:-""};

	mycat $1 > $tmpd/t

	intersectBed -a $tmpd/t -b $2 -wa -wb $opts \
	| awk -v OFS=";" '{ print $1,$2,$3,$4,$5,$6"\t"$8"\t"$11;}' \
	| sort -k1,1 | groupBy -g 1 -c 2,3 -o collapse,collapse > $tmpd/ta

	intersectBed -a $tmpd/t -b $3 -wa -wb $opts \
	| awk -v OFS=";" '{ print $1,$2,$3,$4,$5,$6"\t"$8"\t"$11;}' \
	| sort -k1,1 | groupBy -g 1 -c 2,3 -o collapse,collapse > $tmpd/tb

	join -j 1 -o 0,1.2,1.3,2.2,2.3 $tmpd/ta $tmpd/tb | tr " " "\t" \
	| test_lineartrend - \
	| awk -v OFS="\t" '{ if(substr($1,length($1),1)=="-"){ $6=-$6;} print $0;}' \
	| padjust - -1 

	rm -rf $tmpd;
}


f(){

	grep -v NA $tmpd/tab | test_lineartrend - \
        | awk -v OFS="\t" '{ if( substr($1,length($1),1) == "-"){ $(NF-2) = - $(NF-2);} print $0;}'

	rm -rf $tmpd
}

pa_precomp_fe(){
        tmpd=`mktemp -d`
        mycat $2 | awk -v OFS="@" '{ print $1,$2,$3,$4,$6"\t"$5;}' | sort -k1,1 > $tmpd/a
        mycat $3 | awk -v OFS="@" '{ print $1,$2,$3,$4,$6"\t"$5;}' | sort -k1,1 > $tmpd/b
        join -e 0 -o 0,1.2,2.2 -a 1 -a 2 $tmpd/a $tmpd/b | tr "@ " "\t" | awk -v OFS="\t" '{$4=$4"\t0";}1' > $tmpd/c
        mycat $1 | intersectBed -a stdin -b $tmpd/c -wa -wb -s \
        | perl -ne 'chomp;my @a=split/\t/,$_; print join(",",@a[6..11]),"\t",join(",",@a[0..5]),"\t",$a[12],"\t",$a[13],"\n";' \
        | igx_to_igxnx - \
        | awk -v OFS="\t" '{ print $1"@"$2,$3,$4,$5,$6;}'
        rm -rf $tmpd
}

targetscan(){
usage="$FUNCNAME <miR_Family_Info.txt.zip> <ucsc:targetScanS.txt>";
if [ $# -ne 2 ];then echo "$usage"; return; fi
	## take 3' end point, then switch strand, and then sum scores
	tmpd=`make_tempdir`;
	#tmpd=tmpd;rm -rf $tmpd; mkdir -p $tmpd;
	mycat $1 | tail -n+2 | cut -f 1,4 \
	| perl -ne '$_=~s/\t[\w\d\.-]+:/\t/g; $_=~s/:\d+$//; print $_;' \
	| sort -k1,1 > $tmpd/targetscan_to_mirbase.txt

	mycat $2 | cut -f2-7 \
	| perl -ne '$_=~s/\t[\w\d\.-]+:/\t/g; $_=~s/:\d+$//; print $_;' \
	| sort -k4,4 > $tmpd/targetscan.bed 

	join -1 4 -2 1 $tmpd/targetscan.bed $tmpd/targetscan_to_mirbase.txt \
	| awk -v OFS="\t" '$7 ~ /hsa-/{ print $2,$3,$4,$7,$5,$6;}' 
	rm -rf $tmpd
}


bw(){
	pa=$1; csize=$2; out=$3;
	tmpd=`make_tempdir`
	#tmpd="tmpd"; mkdir -p $tmpd;
	sort_bed $1 | cut -f1-6 > $tmpd/a
	tot=`cat $tmpd/a | perl -ne '@a=split/\t/,$_; $s+=$a[4]; END{ print $s/1000000;}'`
	echo "norm_factor=$tot" > ${out}.info

	awk -v OFS="\t" -v tot=$tot '{if($6=="+"){ print $1,$2,$3,$5/tot;}}' $tmpd/a > $tmpd/p
	bedGraphToBigWig $tmpd/p $csize ${out}_fwd.bw; 

	awk -v OFS="\t" -v tot=$tot '{if($6=="-"){ print $1,$2,$3,$5/tot;}}' $tmpd/a > $tmpd/n 
	bedGraphToBigWig $tmpd/n $csize ${out}_bwd.bw
	rm -rf $tmpd;
}

## stat.sh::sum
pa_point(){
	awk -v OFS=";" '{ st="-";
		if($6 == "-"){ st="+"; $2=$3-1; }
	 	print $1,$2,st"\t"$5; 
	}' $1 | sum - | tr ";" "\t" | awk -v OFS="\t" '{ print $1,$2,$2+1,".",$4,$3; }'
}
test__point(){
echo \
"c	1	4	a1	1	+
c	1	3	a1	10	+
c	1	4	a2	2	-" \
| pa_point - > obs	
echo \
"c	1	2	.	11	-
c	3	4	.	2	+" > exp
check obs exp
rm -f obs exp
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

pa_cluster_sb(){ # snowballing
	sort -k1,1 -k2,3n $1 \
	| mergeBed -i stdin -s -c 5,6 -o sum,distinct -d $2 \
	| awk -v OFS="\t" -v D=$2 '{ print $1,$2,$3,"SB"D,$4,$5;}' 
}
cluster1(){
	MIND=$2;
	local tmpd=`make_tempdir`;
	#local tmpd="tempd"; mkdir -p $tmpd;
	sort -k1,1 -k2,3n $1 | mergeBed -i stdin -s -c 6,2,5 -o distinct,collapse,collapse -d $MIND > $tmpd/a.bed
        cat $tmpd/a.bed | awk -v OFS="\t" -v M=$MIND '
        $3-$2 < M {
                L=split($5,a,","); split($6,b,",");
                n=0; s=0; y=0;
                for(i=1; i<=L; i++){
                        s += a[i] * b[i];
                        n += b[i];
                }
                print $1,$2,$3,int(s/n),int(n),$4;
        }' > $tmpd/b.bed

        cat $tmpd/a.bed | awk -v M=$MIND '$3-$2 >= M' \
         | eval $CLUSTER -d $MIND -c 5,6 - \
         | awk -v OFS="\t" '{
                L=split($7,a,","); split($8,b,","); split($9,c,","); split($10,d,",");
                for(i=1;i<=L;i++){
                        print $1,a[i],b[i],int(c[i]),int(d[i]),$4;
                }
        }' > $tmpd/c.bed
	cat $tmpd/c.bed $tmpd/b.bed
	rm -rf $tmpd
}

cluster_many(){
usage=" $FUNCNAME <mind> <bed> [<bed> ..]"
if [ $# -lt 2 ];then echo "$usage"; return; fi
        MIND=$1;
        local tmpd=`make_tempdir`;
        cat ${@:2} | cluster - $MIND > $tmpd/a
        for f in ${@:2} ;do
                bed_count $tmpd/a $f 1 -s > $tmpd/b
                mv $tmpd/b $tmpd/a
        done
        cat $tmpd/a
        rm -rf $tmpd
}


_countbed(){
	intersectBed -a $1 -b $2 -wa -wb -s \
	| awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$11;}' \
	| groupBy -g 1,2,3,4,6 -c 7 -o sum \
	| awk -v OFS="\t" '{ print $1,$2,$3,$4,$6,$5;}'  
	## zero counts
	intersectBed -a $1 -b $2 -v -s \
	| awk -v OFS="\t" '{ print $1,$2,$3,$4,0,$6;}'
}
_pre_fisherexact(){
	#mkdir -p tmpd; tmpd="tmpd";
	tmpd=`make_tempdir`;
	mycat $1 | sort_bed - > $tmpd/t; 
	mycat $2 | cut -f1-6 | sort_bed - > $tmpd/a; 
	mycat $3 | cut -f1-6 | sort_bed - > $tmpd/b	
	## make clusters w/ the pooled
	cat $tmpd/a $tmpd/b | sum_score - | cluster - $4 > $tmpd/c
	## recount per cluster
	_countbed $tmpd/c $tmpd/a > $tmpd/ca
	_countbed $tmpd/c $tmpd/b > $tmpd/cb
	intersectBed -a $tmpd/ca -b $tmpd/cb -wa -wb -f 1 -r -s \
	| cut -f1-6,11 \
	| intersectBed -a stdin -b $tmpd/t  -wa -wb  -s \
	| awk -v OFS="@" '{ print $1,$2,$3,$4,0,$6"\t"$8,$9,$10,$11,$12,$13"\t"$5"\t"$7;}' 
	rm -rf $tmpd
}
_precompare(){
	#mkdir -p tmpd; tmpd="tmpd";
	tmpd=`make_tempdir`;
	mycat $1 | sort_bed - > $tmpd/t; 
	mycat $2 | cut -f1-6 | sort_bed - > $tmpd/a; 
	mycat $3 | cut -f1-6 | sort_bed - > $tmpd/b	
	## make clusters w/ the pooled
	cat $tmpd/a $tmpd/b | sum_score - | cluster - $4 > $tmpd/c
	## recount per cluster
	_countbed $tmpd/c $tmpd/a > $tmpd/ca
	_countbed $tmpd/c $tmpd/b > $tmpd/cb
	intersectBed -a $tmpd/t -b $tmpd/ca -wa -wb -s \
	| groupBy -g 1,2,3,4,5,6 -c 10,11 -o collapse,collapse  \
	| intersectBed -a stdin -b $tmpd/cb -wa -wb -s \
	| groupBy -g 1,2,3,4,5,6,7,8 -c 12,13 -o collapse,collapse \
	| awk -v OFS="\t" '{ print $1"@"$2"@"$3"@"$4"@"$5"@"$6,$7,$8,$9,$10;}'
	rm -rf $tmpd
}
compare_lineartrend(){
usage="$FUNCNAME <target> <polya_ctr> <polya_trt> <mind>";
	if [ $# -ne 4 ]; then echo "$usage"; return; fi
	_precompare $1 $2 $3 $4 \
	| test_lineartrend - \
	| tr "@" "\t" \
	| awk -v OFS="\t" '{ if($6 == "-"){ $(NF-2) = - $(NF-2);} print $0;}'
}
compare_fisherexact(){
usage="$FUNCNAME <target> <polya_ctr> <polya_trt> <mind>";
	if [ $# -ne 4 ]; then echo "$usage"; return; fi
	echo \
"chr	start	end	name	score	strand	peak_start	peak_end	count1	count2	total1	total2	log2fc	pval	fdr"
	_pre_fisherexact $1 $2 $3 $4 \
	| test_fisherexact -  \
	| tail -n+2 \
	| awk -v OFS="\t" '{ split($1,a,"@"); print $2,a[2],a[3],$3,$4,$5,$6,$7,$8,$9;}' \
	| tr "@" "\t"
}

batch_bam-to-point(){
usage="$FUNCNAME <batchscript>"
	if [ $# -lt 1 ];then echo "$usage"; return; fi
	eval `cat $1 | perl -ne 'chomp;$_=~s/#.*//g; $_=~s/^\s+$//g; print $_,"\n" unless $_ eq "";'`
}


batch_polya(){ 
usage="$FUNCNAME <batchscript> [test]"
	if [ $# -lt 1 ];then echo "$usage"; return; fi
	eval `cat $1 | perl -ne 'chomp;$_=~s/#.*//g; $_=~s/^\s+$//g; print $_,"\n" unless $_ eq "";'`
	BAM=( $BAM )
	COMP=( $COMP )
	
	if [ -d $OUT ];then
		echo "$OUT exists. overriding on this .. " >&2
	fi
	mkdir -p $OUT/comp $OUT/point $OUT/bw
	## todo: check files I trust you
	get_chromsize ${BAM[1]} > $OUT/chrom.size	
	for (( i=0; i < ${#BAM[@]}; i+=2 ));do
		name=${BAM[$i]}; bam=${BAM[$i+1]};
		outd=$OUT/point/$name; mkdir -p $outd
		if [ -f $outd/a.bed ];then
			echo "$outd/a.bed exists .. pass " >&2
			continue;
		fi
		if [[ $# -gt 1 &&  $2 = "test" ]];then
			TEST="head -n 10000";
		else
			TEST="cat";	
		fi
		echo "running $bam => $outd/a.bed .. " >&2
		samtools view -bq 10 $bam | bamToBed \
		| eval $TEST \
		| modify_score - "count"\
		| point - \
		| filter - $FASTA > $outd/a.bed 
		bw $outd/a.bed $OUT/chrom.size $OUT/bw/$name
		
	done
	if [[ -f $TARGET && ${#COMP[@]} -gt 1 ]];then
		for (( i=0; i < ${#COMP[@]}; i+=2 ));do
			point1=$OUT/point/${COMP[$i]}/a.bed
			point2=$OUT/point/${COMP[$i+1]}/a.bed
			outd=$OUT/comp/${COMP[$i]}_vs_${COMP[$i+1]}; mkdir -p $outd;

			echo "comparing ${COMP[$i]} vs ${COMP[$i+1]}, MDIST=$MDIST => $outd/lineartrend.txt " >&2
			compare_lineartrend $TARGET $point2 $point1 $MDIST > $outd/lineartrend.txt
		done
	else
		echo "$TAGRET not exists">&2
	fi
	
}

pre_test_fisherexact(){
usage="
usage $FUNCNAME <target> <point_ctr> <point_trt> [ <mind> [ <strand> ] ]
	mind: minimum distance between weighted centers of clusters (default 20)
	strand: intersectBed option for target and points 
		(-s : count points in the same strand as target)
"
	if [ $# -lt 3 ]; then echo "$usage"; return; fi
	target=$1; point_ctr=$2; point_trt=$3;  mind=$4; opt_strand=$5; 
	## options: mind=20; target_strand="-s"
	## make id group count_ctr count_trt
	mind=${mind:-20}; opt_strand=${opt_strand:-};
	#tmpd="tmpd"; mkdir -p $tmpd; 
	tmpd=`make_tempdir`;
	mycat $1 > $tmpd/t;
	mycat $2 | cut -f1-6 > $tmpd/a;
	mycat $3 | cut -f1-6 > $tmpd/b;
	## make clusters w/ the pooled
	cat $tmpd/a $tmpd/b | sum_score - | cluster - $mind > $tmpd/c
	## recount per cluster
	bed_count $tmpd/c $tmpd/a 1 $opt_strand \
	| bed_count - $tmpd/b 1 $opt_strand \
	| intersectBed -a $tmpd/t -b stdin -wa -wb  $opt_strand \
	| perl -ne 'chomp; my @a=split/\t/,$_; 
		print join("@",@a[6..11]),"\t",join("@",@a[0..5]),"\t",$a[12],"\t",$a[13],"\n";' 
	rm -rf $tmpd
}

cryptic_vs_cannonical(){
	gene=$1; utr=$2; ctr_point=$3; trt_point=$4;

	bed_count $gene $ctr_point \
	| awk -v OFS="@" '{ print $1,$2,$3,$4,$6"\t"$7;}' \
	| sort -k1,1 \
	> a
	bed_count $gene $trt_point \
	| awk -v OFS="@" '{ print $1,$2,$3,$4,$6"\t"$7;}' \
	| sort -k1,1 \
	> b
	join a b | tr " " "\t" | awk -v OFS="\t" '{ print "NC",$0;}' > c
	rm a b

	mycat $utr | intersectBed -a $ctr_point -b stdin -v -s | bed_count $gene -  \
	| awk -v OFS="@" '{ print $1,$2,$3,$4,$6"\t"$7;}' \
	| sort -k1,1 \
	> a
	mycat $utr | intersectBed -a $trt_point -b stdin -v -s | bed_count $gene -  \
	| awk -v OFS="@" '{ print $1,$2,$3,$4,$6"\t"$7;}' \
	| sort -k1,1 \
	> b
	join a b | tr " " "\t" | awk -v OFS="\t" '{ print "CR",$0;}' >> c
	rm a b
	test_fisherexact c 
}

###########################################################
# test 
###########################################################
test(){
echo \
"chr1	1	100	 n1	1	+
chr1	200	300	 n2	1	+" > t

echo \
"chr1	1	2	 a1	1	+
chr1	3	4	 a1	1	+
chr1	7	8	 a1	1	+
chr1	200	201	 a2	2	+" > a

echo \
"chr1	1	2	 a1	1	+
chr1	4	5	 a1	1	+
chr1	7	8	 a1	1	+
chr1	203	204	 a2	2	+" > b

#_pre_fisherexact t a b  2 
compare_fisherexact t a b 2
rm t a b

echo "test .. point inp"
echo \
"chr1	95	96	.	2	-
chr1	99	100	.	2	+
chr2	204	205	.	2	+" > exp

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
"chr1	1	200	g1	0	-	5,99	3,3	5,99	120,1	0.597334	2.01324512616e-11" > exp

compare_lineartrend c.bed a.bed b.bed 50 > obs
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

#typeset -F
