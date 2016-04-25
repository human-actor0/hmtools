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
pa.cluster(){
usage(){
echo "
USAGE: $FUNCNAME [options] <bed> [<bed>..]
 [options]: 
	-d <maxd> : maximum distance between peaks (10 default)
	-h : print header (default none)
"
}
local d OPTIND; d=10; h=0;
while getopts ":d:h" opt ;do
	case $opt in
		d) d=$OPTARG;;
		h) h=1;;
		\?) usage;;
	esac
done
shift $(( OPTIND -1 ))
if [ $# -lt 1 ];then usage; return; fi
	if [ $h = 1 ];then
		echo "# d=$1 files=[ ${@:1} ]";
	fi
	local tmpd=`mymktempd`;
	local files=( ${@:1} );	
	local chroms="";
	for (( i=0; i< ${#files[@]}; i++ ));do
		mkdir $tmpd/$i
		bed.split ${files[$i]} $tmpd/$i
		chroms="$chroms `ls $tmpd/$i`"
	done
	chroms=( `echo "$chroms" | tr " " "\n" | sort -u` );
	for chrom in ${chroms[@]};do
		cat $tmpd/*/$chrom | bed.sort - \
		| mergeBed -i stdin -s -c 2,3,5 -o collapse,collapse,collapse -d $d \
		| perl -ne 'chomp; my @a=split/\t/,$_;
			my @s=split/,/,$a[4]; my @e=split/,/,$a[5];
			my @x=(); for(my $i=0;$i<=$#s; $i++){ $x[$i]= int(($s[$i] + $e[$i]  - 0.5)/2); }
			my @y=split/,/,$a[6];
			my $n=0; my $s=0; my $p=0;
			for(my $i=0; $i <= $#x; $i++){
				$n++; $s+=$y[$i]; $p+=$y[$i]*$x[$i];
			}
			my $center=int( 0.5*($a[1]+$a[2]));
			if($s > 0){ $center=int( $p/$s + 0.5); }
			my $score=sprintf("%.2f",$s);
			print $a[0],"\t",$a[1],"\t",$a[2],"\t",$center,"\t",$score,"\t",$a[3],"\n";
		' > $tmpd/$chrom
	done

	for chrom in ${chroms[@]};do
		for (( i=0; i < ${#files[@]}; i++ ));do
			if [ -f $tmpd/$i/$chrom ];then
				intersectBed -a $tmpd/$chrom -b $tmpd/$i/$chrom -wa -wb -s \
				| bed.enc - | awk -v OFS="\t" -v i=$i '{ print $4","i,$11;}' | stat.sum -  >> $tmpd/s.$chrom
			fi
		done

		cat $tmpd/s.$chrom | perl -e 'use strict; my %res=(); my $n='${#files[@]}';
			while(<STDIN>){ chomp; 
				my ($idi,$s) = split /\t/,$_;
				my ($id,$i)=split /,/,$idi;
				$res{$id}{$i} = $s;
			}	
			foreach my $id (keys %res){
				print $id;
				for(my $i=0; $i<$n; $i++){
					my $v= defined $res{$id}{$i} ? $res{$id}{$i} : 0; 
					print "\t$v";	
				}
				print "\n";
			}
		' | tr "@" "\t"
	done
	rm -rf $tmpd	
}
pa.cluster.test(){
echo \
"c	1	2	n1	1	-
c1	1	2	n1	1	-
c	1	2	p1	1	+" > tmp.a
echo \
"c	2	3	p2	10	+
c	20	30	p4	0	+
c	4	5	p3	10	+" > tmp.b
	pa.cluster -d 3 tmp.a tmp.b
	rm tmp.*
}

pa.cmp(){
usage(){ echo \
" $FUNCNAME [options] <3utr.bed> <ctr.bed> <trt.bed>
	[options]:
	 -d <int> : maximum distance between peaks
"
return
}
local d OPTIND; d=10;

while getopts ":d:" opt ;do
	case $opt in
		d) d=$OPTARG;;
		\?) usage;;
	esac
done
shift $(( OPTIND -1 ))
if [ $# -ne 3 ];then usage; fi
	pa.cluster -d $d $2 $3 | bed.group - $1 \
	| perl -e 'use strict; my %res=();
	while(<STDIN>){ chomp; my @a=split/\t/,$_;
		if( $a[$#a] ne "NA"){
			$res{$a[$#a]}{$a[3]}=[$a[6],$a[7]];
		}
	}		
	foreach my $i (keys %res){
		my @x=sort {$a<=>$b} keys %{$res{$i}};
		my @y1=map {$res{$i}{$_}->[0]} @x;
		my @y2=map {$res{$i}{$_}->[1]} @x;
		print $i,"\t",join(",",@x),"\t",join(",",@y1),"\t",join(",",@x),"\t",join(",",@y2),"\n";
	}' | stat.lineartrend - \
	| tr "@" "\t" | awk -v OFS="\t" '{ if($6=="-"){ $11=-$11;} } 1';
}
pa.cmp.test(){
echo \
"c	100	200	gene1	0	+
c	300	400	gene2	0	+" > tmp.3utr

echo \
"c	10	11	pa0	6	+
c	110	111	pa1	5	+
c	120	121	pa2	4	+
c	150	151	pa3	3	+
c	310	311	pa4	2	+
c	399	400	pa5	1	+" > tmp.ctr

echo \
"c	10	11	pa0	10	+
c	110	111	pa1	20	+
c	120	121	pa2	30	+
c	150	151	pa3	40	+
c	310	311	pa4	50	+
c	399	400	pa5	60	+" > tmp.trt
head tmp.*
pa.cmp -d 2 tmp.3utr tmp.ctr tmp.trt 
rm tmp.*
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

