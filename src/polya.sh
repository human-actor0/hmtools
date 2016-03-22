#/bin/bash  
. $HMHOME/src/root.sh # import utilities 
. $HMHOME/src/bed.sh  #import utilities 
. $HMHOME/src/stat.sh # import test_lineartrend test_fisherexact
. $HMHOME/src/seq.sh
. $HMHOME/src/ucsc.sh;

SEQ=$HMHOME/bin/bed_seq.sh
FILTER=$HMHOME/src/pa_filter_nb.sh
FILTER_M=$HMHOME/src/nb.model
HG19FA=/hmdata/ucsc/hg19/chromosome/
CLUSTER=$HMHOME/src/cluster_1d_kmeans.sh
BW=$HMHOME/bin/bedGraphToBigWig

pa.comp_pcpa(){
usage="
FUNCT: relative FC  
USAGE: $FUNCNAME <gene.bed> <3utr> <trt.bed> <ctr.bed> <s> <fc>
"; if [ $# -lt 6 ]; then echo "$usage"; return; fi

	echo "gene@cluster trt.count trt.other ctr.count ctr.other FC sign" | tr " " "\t"
	pa.preptest $1 $3 $4 \
	| perl -ne 'chomp; my @a=split/\t/,$_;
		my ($g,$p) = split /@/,$a[0];	
		$p=~s/;/\t/g;
		print $p,"\t",join("\|",@a),"\n";
	'| tail -n+2 | intersectBed -a stdin -b $2 -wa -v -s | cut -f7- | tr "|" "\t" \
        | awk -v OFS="\t" -v S=$5 -v FC=$6 '{
                e="U"; eps=1; s1=$2+$4; s2=$3+$5;
                fc=($2+eps)*($3+$5+eps)/($2+$4+eps)/($3+eps);
                if( $2+$3>= S && $2+$4 >=S && $3+$5 >=S ){
                        if(fc >= FC){ e="P";}
                        else if(fc <= 1/FC){ e="N";}
                }else{ e="I";}
                print $0,fc,e;
        }'
}
pa.bw(){
usage="
FUNCT : make bigwing files from bed6
USAGE : $FUNCNAME <bed> <genome> <norm> <out> 
	[norm] : 0: no, 1: RPBM (read per bin per million reads)
	<genome> : hg19 
"
if [ $# -ne 4 ];then echo "$usage"; return; fi
	local tmpd=`mymktempd`
	awk -v OFS="\t" -v O=$tmpd '{
		fout=O"/fwd"; if($6=="-"){ fout=O"/bwd";}
		print $1,$2,$3,$5 >> fout;
	}' $1; 
	local tot=`cat $tmpd/*wd | wc -l`;
	echo $tot;
	ucsc.chrom_size $2 > $tmpd/size
	for f in $tmpd/*wd;do
		local n=${f##*/};
		bed.sort $f | awk -v OFS="\t" -v tot=$tot -v norm=$3 '{
			f=0; if(tot>0){ f=1000000/tot;}
			v=$4; if ( norm ){ v=v*f;}
			print $1,$2,$3,v;		
		}' > $f.1
		$HMHOME/bin/bedGraphToBigWig $f.1 $tmpd/size ${4}_${n}.bw
	done
	rm -rf $tmpd;
}
pa.bw.test(){
echo \
"c	100" > tmp.c

echo \
"c	1	2	.	1	+
c	3	4	.	20	-
c	5	6	.	10	-
c	7	8	.	2	+" | pa.bw - tmp.c 1 tmp.o 
ls -la tmp.*.bw
rm -rf tmp.* 
}

pa.cluster_sb(){ 
usage="
USAGE: $FUNCNAME <bed> <maxd> 
	<maxd> : maximum distance between features to be merged 
"; if [ $# -lt 2 ]; then echo "$usage"; return; fi
## mergeBed version issue 
	local tmpd=`mymktempd`;
	bed.split $1 $tmpd
	for f in $tmpd/*;do
		sort -k1,1 -k2,3n $f \
		| mergeBed -i stdin -s -c 2,5 -o collapse,collapse -d $2 \
		| perl -ne 'chomp; my @a=split/\t/,$_;
			my @x=split/,/,$a[4];
			my @y=split/,/,$a[5];
			my $n=0; my $s=0; my $p=0;
			for(my $i=0; $i <= $#x; $i++){
				$n++; $s+=$y[$i]; $p+=$y[$i]*$x[$i];
			}
			my $center=int( $p/$s + 0.5);
			my $score=sprintf("%.2f",$s);
			print $a[0],"\t",$a[1],"\t",$a[2],"\t",$center,"\t",$score,"\t",$a[3],"\n";
		'
	done
	rm -rf tmpd;
}
pa.cluster_sb.test(){
echo \
"c	10	11	r1	1	+
c	12	13	r2	2	+
c	13	14	r3	3	+
c	6	7	r4	0.01	+" | pa.cluster_sb - 1 
}


pa.point(){
usage="
$FUNCNAME <bam> <q>
"; if [ $# -ne 2 ];then echo "$usage"; return; fi

	samtools view -bq $2 $1 | bamToBed -split \
	| bed.5p - | bed.ss - | awk -v OFS="@" '{ print $1,$2,$3,$6"\t"1;}' \
	| stat.sum -  | tr "@" "\t" | awk -v OFS="\t" '{ print $1,$2,$3,".",$5,$4;}'
}

pa.filter(){ 
usage="
FUNCT: Filter inter-priming artefacts 
USAGE: $FUNCNAME <bed> <fasta>
REFER: heppard, S., Lawson, N.D. & Zhu, L.J., 2013.Bioinformatics (Oxford, England), 29(20), pp.2564–2571.
"
	if [ $# -ne 2 ]; then echo "$usage"; return; fi

	bed.flank $1 39 30 -s  | seq.read $2 - -s \
	| perl -ne 'chomp; my @a=split/\t/,$_; print join(";",@a[0..($#a-1)]),"\t",$a[$#a],"\n"' \
	| eval "$FILTER predict - $FILTER_M" \
	| perl -ne 'chomp; my ($bed,$seq,$pos,$score)=split/\t/,$_;
		my @a=split /;/,$bed;
		my $len=$a[2]-$a[1];
		my @p=split /,/,$pos;
		my @s=split /,/,$score;
		for(my $i=0; $i<=$#p;$i++){
			my $this_s=sprintf("%.3f",$s[$i]); 
			my $this_p=$p[$i];
			if ($#a > 4 && $a[5] eq "-"){
				$this_p=$len-$this_p-1;
			}
			print join("\t",($a[0],$a[1]+$this_p,$a[1]+$this_p+1,$a[3],$a[4],$a[5],$this_s)),"\n";
		}
	'
}

pa.filter.test(){
fa=~/hmdata/ucsc/hg19/chromosome/chr1.fa.gz 
#echo \
#"chr1	11117510	11117550	n1	0	+	
#chr1	11117510	11117550	n1	0	-" \
echo \
"chr1	182550600	182550620	n1	0	+
chr1	182550600	182550620	n1	0	-" \
| pa.filter - $fa 
}

