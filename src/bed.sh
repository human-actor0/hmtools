#!/bin/bash
. $HMHOME/src/root.sh
. $HMHOME/src/stat.sh

bed.sort(){
 	sort -i -k 1,1 -k 2n,2 $1 
}
bed.split(){
usage="
USAGE: $FUNCNAME <bed> <outdir>
"; if [ $# -ne 2 ];then echo "$usage"; return; fi
	if [ -d $2 ];then rm -rf $2; fi; mkdir -p $2;
	awk -v OFS="\t" -v O=$2 '{
		fout=O"/"$1;
		print $0 >> fout;
	}' $1
}

bed.merge(){
## -s option performs differently in different versions
usage="$FUNCNAME <bed> [mergeBed options]
"
if [ $# -lt 1 ];then echo "$usage"; return; fi
	local opts=${@:2};
	local opts_s="";
	if [ -z ${opts##*\-s*} ];then 
		opts_s=1; 
		opts=${opts/\-s/};
	fi
	if [ -z $opts_s ];then
		cat $1 | mergeBed -i stdin $opts 
	else
		awk -v OFS="\t" '{ $1=$1","$6;print $0;}' $1 | mergeBed -i stdin $opts \
		| awk -v OFS="\t" '{ split($1,a,","); $1=a[1]; $3=$3"\t"a[2]; print $0;}' 
	fi
}
bed.merge.test(){
echo \
"c	1	2	n	0	+
c	2	3	n	0	+
c	2	4	n	0	-" \
| bed.merge - -d -1 -s > obs
echo \
"c	1	3	+	n	0
c	2	4	-	n	0" > exp
cat obs
#check obs exp

}


bed.nf(){
	head -n 1 $1 | awk '{print NF;}'
}

bed.eachchrom(){
usage="
USAGE: $FUNCNAME <bed> '<functions>'
"; if [ $# -ne 2 ];then echo "$usage"; return; fi
	local tmpd=`mymktempd`;
	bed.split $1 $tmpd;
	for f in $tmpd/*;do
		cat $f | eval "$2";	
	done
	rm -rf $tmpd
}
bed.eachchrom.test(){
echo \
"c1	11	22	n	0 	+
c1	3	4	n	0 	+
c2	1	2	n	0 	+" | bed.eachchrom - 'sort -k1,1 -k2,2n'

}


bed.count(){
usage="
usage: $FUNCNAME <target> <read>  [options]
output: target + sum of read scores
 [options]: 
	-s : count on the same strand
	-S : count on the opposite strand
	-a : report all target and put zero counts
"
	local opt_strand="";
	local opt_all="";
	for opt in $@;do
		if [ $opt = "-s" ];then opt_strand="-s";
		elif [ $opt = "-S" ];then opt_strand="-S";
		elif [ $opt = "-a" ];then opt_all="-a";
		fi
	done
	if [ $# -lt 2 ];then echo "$usage"; return; fi

	#local tmpd=tmpd;mkdir -p $tmpd;
	local tmpd=`mymktempd`;
	mycat $1 | bed.split - $tmpd/a
	mycat $2 | bed.split - $tmpd/b
	for f in $tmpd/a/*;do
		c=${f##*/};
		if [ -f $tmpd/b/$c ];then
			local nf=`bed.nf $f`; 
			intersectBed -a $f -b $tmpd/b/$c -wa -wb $opt_strand  \
			| awk -v nf=$nf -v OFS="\t" '{ 
				for(i=2; i<=nf;i++){ $1=$1"@"$(i); }
				print $1,$(nf+5);
			}' | stat.sum - | tr "@" "\t" 

			if [ ! -z $opt_all ];then
				intersectBed -a $f -b $tmpd/b/$c -v -s  \
				| awk -v OFS="\t" '{ print $0,0;}'
			fi
		elif [ ! -z $opt_all ];then
			awk -v OFS="\t" '{ print $0,0;}'
		fi
	done
	rm -rf $tmpd
}
bed.count.test(){
echo \
'chr	1	100	n1	0	+
chr	50	200	n2	0	-
chr	1000	2000	n3	0	+'> a
echo \
'chr	1	10	r1	1	+
chr	40	50	r2	2	+
chr	50	200	r3	3	+' > b

bed.count a b -s -a 
rm -rf a b 
}

bed.mcount(){
usage="
FUNCT: make a merged bed and then recount 
<<<<<<< HEAD
USAGE: $FUNCNAME <dist> <bed> [<bed>,...] 
 [options]:
	<dist> : maximum distance to be merged
"
	if [ $# -lt 2 ];then echo "$usage"; return; fi
	local dist=$1; 
	local tmpd=`mymktempd`;
	cat "${@:2}" | bed.sort - | bed.merge - -d $dist -s -c 4 -o distinct \
	| awk -v OFS="\t" '{ print $1,$2,$3,$5,0,$4;}' > $tmpd/a
	for f in ${@:2};do
		bed.count $tmpd/a $f -s -a > $tmpd/b
		mv $tmpd/b $tmpd/a
	done
	cat $tmpd/a
	rm -rf $tmpd;
}
bed.mcount.test(){
echo  \
"c	1	2	n	1	+
c	4	5	n	2	+" > tmp.a
echo  \
"c	1	2	n	11	+
c	2	3	n	22	+
c	4	6	n	33	+" > tmp.b
bed.mcount 0 1 tmp.a tmp.b > obs 
echo \
"c	1	3	.	0	+	1	33
c	4	6	.	0	+	2	33" > exp
cat obs
#check exp obs
rm tmp.a tmp.b exp obs
}


bed.n2i(){
usage="$FUNCNAME <bed>" 
if [ $# -ne 1 ];then echo "$usage";return; fi
	local tmpd=`mymktempd`;
	#local tmpd=tmpd; mkdir -p tmpd;
	cat $1 | perl -e 'use strict; 
		my $fileb="'$tmpd/b'";
		my %ref=(); 

		my $id=0; my $first=1;
		while(<STDIN>){chomp; my @a=split/\t/,$_;
			if(defined $ref{$a[3]}){
				$ref{$a[3]}{n} ++;
			}else{
				if($first==1){ $first=0;
				}else{ $id++; }
				$ref{$a[3]}{n} = 1;
				$ref{$a[3]}{id} = $id;
			}
			$a[3]=$ref{$a[3]}{id};
			print $a[3],"\t",join (";",@a),"\n";
		}
		
		open(my $fh, ">",$fileb) or die "cannot open $fileb: $!";	
		foreach my $k (keys %ref){
			print $fh $ref{$k}{id},"\t",$ref{$k}{n},"\n";	
		}
		close($fh);
	' | sort -k1,1n > $tmpd/a
	rm -rf $tmpd
}

bed.n2i.test(){
echo \
"c	1	2	r3	1	+
c	1	2	r1	1	+
c	1	2	r2	1	+
c	1	2	r3	1	+
c	1	2	r2	1	+
c	1	2	r3	1	+" | bed.n2i - 0
}

bed.mhits(){

usage="
FUNCT: redistribute multi-hits using EM algorithm
USAGE: $FUNCNAME <target> <multi_reads> [-s|-S]
"; if [ $# -lt 2 ];then echo "$usage"; return; fi
local MAXITER=100;
local DELTA=0.001;
local OPT=${3:-""};

	intersectBed -a $1 -b $2 -wa -wb $OPT \
	| perl -e 'use strict; 
		my $MAXITER='$MAXITER';
		my $DELTA='$DELTA';

		my %E=(); # relative,normalized gene expression frequency
		my %A=(); # a fraction of a read attributed to a gene 
		my %L=(); # gene lengths
		my $Lread=0; my $Ltot=0;
		my $TOT=0; # total weighted reads
		my $s=0;

		sub ss{ ## root sum of squares
			my ($h1, $h2)=@_;
			my $s=0;
			foreach my $k (keys %$h1){
			if( defined $h2->{$k} ){
				$s += $h1->{$k} * $h2->{$k};
			}}	
			return $s;
		}
		sub min{
			my ($x,$y)=@_;
			if($x > $y){ return( $y);}
			return $x;
		}
		sub max{
			my ($x,$y)=@_;
			if($x > $y){ return( $x);}
			return $y;
		}

		sub normE{
			my ($E, $L, $Lread) = @_;
			## normalize by effective gene lengths
			my $s=0;
			foreach my $g (keys %$E){
				my $len = $L->{$g} - $Lread + 1;
				if($len >0){ $E->{$g} /= $len;
				}else{ $E->{$g} = 0;}
				$s += $E->{$g};
			}
			## normalize gene expressions
			foreach my $g (keys %$E){ $E->{$g} /= $s;}
		}
		sub estE{
			my ($Ein, $A,      $Eout) = @_;
			## I do not think of separating uniq and multi reads
			foreach my $r (keys %$A){
				my $s=0; my $norm=0;
				foreach my $g (keys %{$A->{$r}}){ $s += $Ein->{$g}; }	
				if($s > 0){ $norm = 1/$s; }
				foreach my $g (keys %{$A->{$r}}){
					$Eout->{$g} += $Ein->{$g} * $norm;
				}
			}
		}
		sub update{
			my ($Ein, $A,$L,$Lread,     $Eout) = @_;

			## estimate expression by reads
			estE( $Ein, $A, $Eout);
			normE( $Eout, $L, $Lread);
		}

		my $minStep0 = 1.0;
		my $minStep = 1.0;
		my $maxStep0 = 1.0;
		my $maxStep = 1.0;
		my $mStep = 4.0;
		my $OPT_TOL=0.01;


		while(<STDIN>){chomp;
			my ($chr,$start,$end,$gid,$score,$strand, $chr2,$start1,$end1,$rid,$score1,$strand1)=split /\t/,$_;
			$Lread += $end1-$start1; $Ltot ++; 
			$L{$gid}=$end-$start;
			$A{$rid}{$gid} = $score1;
			$E{$gid} = 0;
		}

		## estimated read length	
		$Lread /= $Ltot;

		## init with 1/n weighted read counts 
		foreach my $r (keys %A){
			my $n=scalar keys %{$A{$r}};
			foreach my $g (keys %{$A{$r}}){ $E{$g} += 1/$n; }
		}
		normE(\%E,\%L,$Lread);
		my %E1=(); 
		my %E2=(); 
		my %Eprime=();
		my $Eres = "";
		my %r=(); my %r2=(); my %v=();
		foreach my $iter ( 1..$MAXITER){
			print "ITER--$iter: \n";
			update(\%E,\%A,\%L,$Lread,\%E1);
			update(\%E1,\%A,\%L,$Lread,\%E2);

			foreach my $g ( keys %E){
				$r{$g} = $E1{$g} - $E{$g};
				$r2{$g} = $E2{$g} - $E1{$g};
				$v{$g} = ($E2{$g} - $E1{$g})-$r{$g};
			}
			my $rNorm = sqrt(ss( \%r,\%r));
			my $r2Norm = sqrt(ss( \%r2,\%r2));
			my $vNorm = sqrt(ss( \%v,\%v));
			my $rvNorm = sqrt(abs(ss(\%r,\%v)));
			if ( $vNorm == 0 ){
				print {*STDERR}	"$rNorm,$r2Norm,$vNorm, vnorm == 0\n"; last;	
			}
			my $alphaS=$rNorm/$rvNorm;
			my $alphaS=max($minStep, min($maxStep,$alphaS));

			if ( $rNorm < $OPT_TOL ){
				print {*STDERR} "rNorm=$rNorm \n"; last;
			}
			if ( $r2Norm < $OPT_TOL ){
				print {*STDERR} "r2Norm=$r2Norm \n"; last;
				$Eres = \%E2;
			}
		
			## extraplate 
			foreach my $g (keys %E){
				$Eprime{$g} = max(0.0, $E{$g} + 2 * $alphaS * $r{$g} + ($alphaS*$alphaS)*$v{$g});
			}	
			update(\%Eprime,\%A,\%L,$Lread,\%E);

			## stabilization
			if (abs( $alphaS - 1.0 ) > 0.01 ){
				print {*STDERR} "stablizaation step:\n";
				if ($alphaS == $maxStep){ 
					$maxStep=max($maxStep0, $maxStep/$mStep);
					$alphaS = 1.0;
				}
			}
			if ($alphaS == $maxStep){ $maxStep = $mStep * $maxStep;}
			if ($minStep < 0 && $alphaS == $minStep ){ $minStep = $mStep * $minStep;}
			print "alpha=$alphaS\n";

		}
	'
}
bed.mhits.test(){
echo \
"c	0	100	t1	0	+
c	100	200	t2	0	+
c	300	400	t3	0	+" > genes

echo \
"c	10	21	r1	1	+
c	10	21	r4	1	+
c	300	310	r1	1	+
c	110	120	r2	1	+
c	310	320	r2	1	+
c	13	23	r3	1	+
c	150	170	r3	1	+" > reads

bed.mhits genes reads
rm -rf genes reads
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

bed.flank(){
usage="
FUNCT : extract flanking regions
USAGE : $FUNCNAME <bed> <left> <right> <strand_opt>
	<strand_opt>: [1: strand specific]
" 
if [ $# -ne 4 ];then echo "$usage"; return; fi
	awk -v OFS="\t" -v L=$2 -v R=$3 -v S=$4 '{ 
		if((S == "-s" || S == 1) && $6 == "-"){
			$2=$2-R; $3=$3+L;
		}else{ $2=$2-L; $3=$3+R; } 
		if( $2 < 0){ $2=0;} if( $3 < 0){ $3=1;}
		print $0;
	}' $1;
}

bed.ss(){
 	awk -v OFS="\t" '{ s="-"; if($6=="-"){ s="+";} $6=s; }1' $1
}
bed.3p(){
 	awk -v OFS="\t" '{ if($6=="-"){$3=$2+1;} $2=$3-1; print $0; }' $1
}

bed.5p(){
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

bed.mapq_to_acc(){
#http://samtools.github.io/hts-specs/SAMv1.pdf
	mycat $1 | awk -v OFS="\t" '{ if($5 > 0){ $5 = 1- exp( - $5/10 * log(10));} print $0;}'
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
bed.score(){
usage="
FUNCT: process score fields
USAGE: $FUNCNAME <bed6> <method>
	<method>: 0: do nothing, 1: make it 1, 2: calculate accuracy assuming mapq
"
	if [ $# -ne 2 ]; then echo "$usage"; return; fi
	awk -v OFS="\t" -v ME=$2 '{
		if(ME==1){
			$5=1;
		}else if(ME==2){
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
	modify_score $@
}
bed.sum(){
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
bed.sum.test(){
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

bed.exon(){
usage="
USAGE: $FUNCNAME <bed12>
"
if [ $# -ne 1 ];then echo "$usage"; return; fi
	mycat $1 | awk -v OFS="\t" '{ split($11,sizes,",");split($12,starts,",");
		for(i=1; i<=$10; i++){
	    		s=$2+starts[i]; e=s+sizes[i];
			print $1,s,e,$4,$5,$6;
		}
	}'
	
}


bed.intron(){
usage="
USAGE: $FUNCNAME <bed12>
"
if [ $# -ne 1 ];then echo "$usage"; return; fi
	##  [es   ]s----e[    ee]
	mycat $1 | awk -v OFS="\t" '$10 > 1{
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
	}' 
}

bed.flat(){ 
usage="$FUNCNAME [options] <bed6>
     input:
     [     )
        [      )
     output:
     [  )  [   )
        [  )
"; if [ $# -lt 1 ]; then echo "$usage";exit; fi
	bed.sort $1 | mergeBed -i stdin -d -1 -s -c 2,3,4,5 -o distinct,distinct,distinct,distinct \
	| perl -ne ' chomp; my @a=split /\t/,$_; 
		#print join("\t",@a),"\n";
		my %h=();
		my @s=split/,/,$a[4];
		my @e=split/,/,$a[5];
		if( scalar @s == 1){
			print $a[0],"\t",$a[1],"\t",$a[2],"\t",$a[6],"\t",$a[7],"\t",$a[3],"\n";
			next;
		}
		for( my $i=0;$i<=$#s; $i++){
			$h{$s[$i]}=0;
			$h{$e[$i]-1}=1;
		}
		my @k=sort {$a<=>$b} keys %h;
		for( my $i=0; $i<$#k; $i++){
			my $s= $k[$i] + $h{ $k[$i] };
			my $e= $k[$i+1] + $h{ $k[$i+1] };
			print $a[0],"\t$s\t$e\t",$a[6],"\t",$a[7],"\t",$a[3],"\n";
		}

	' 
}

bed.flat.test(){

echo \
"chr	1	10	r1	1	+
chr	3	13	r2	3	+
chr	15	17	r3	100	+" > tmp.inp
echo \
"chr	1	3	r1,r2	1,3	+
chr	3	10	r1,r2	1,3	+
chr	10	13	r1,r2	1,3	+
chr	15	17	r3	100	+" > tmp.exp

echo "
INPUT:
	123      10 13
	[         )
          [          )  [ )
"
cat tmp.inp

echo "
OUTPUT:
	[ )       [  )   
	  [       )     [ )
";
cat tmp.exp
echo "
RESULT:"
bed.flat tmp.inp > tmp.obs
check tmp.obs tmp.exp
rm tmp.obs tmp.exp tmp.inp

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
	    chomp;
	    my @aa = split /\t/,$_;
	    my ($bin,$name,$chr,$strand,$start,$end,$thickStart,$thickEnd,$blockCount,$blockStarts,$blockEnds,$id,$name2) = split /\t/, $_;
	    my $itemRgb = "255,0,0";
	    my $score = 0;
		print $blackCount,"\n"; exit;

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


