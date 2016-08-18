#!/bin/bash
. $HMHOME/src/root.sh;
. $HMHOME/src/stat.sh
. $HMHOME/src/bed.sh


asp.m3ss_test(){
cat $1 | perl -e 'use strict; my %H=(); 
	my $ctr="'$2'"; my $trt="'$3'";
	while(<STDIN>){chomp; my @a=split/\t/,$_; 
		print join("@",@a);
		if(scalar keys %H==0){ 
			for( my $i=5; $i<= $#a; $i++){
				if($a[$i] =~ /^$ctr/){
					$H{ctr}{$i}=$a[$i];
				}
				if($a[$i] =~ /^$trt/){
					$H{trt}{$i}=$a[$i];
				}
			}
			foreach my $k (keys %{$H{ctr}}){
				print "\t",$H{ctr}{$k},".c1\t",$H{ctr}{$k},".c2"
			}
			foreach my $k (keys %{$H{trt}}){
				print "\t",$H{trt}{$k},".c1\t",$H{trt}{$k},".c2"
			}
			print "\n";
			#print join (",",keys %{$H{ctr}}),"\n";
			#print join (",",keys %{$H{trt}}),"\n";
			next;
		}	
		my @x=split/,/,$a[2]; # 5 ss
		foreach my $i (keys %{$H{ctr}}){
			my ($ta,$tb)=split/;/,$a[$i];
			my @s=split/,/,$ta;
			my @u=split/,/,$tb;
			my $usum=0; my $ssum=0;
			for(my $j=0; $j<=$#x; $j++){
				$ssum += $s[$j];	
				if($x[$j] == $a[1]){
					$usum = $u[$j];
				}
			}
			print "\t",$usum,"\t",$ssum;
		}
		foreach my $i (keys %{$H{trt}}){
			my ($ta,$tb)=split/;/,$a[$i];
			my @s=split/,/,$ta;
			my @u=split/,/,$tb;
			my $usum=0; my $ssum=0;
			for(my $j=0; $j<=$#x; $j++){
				$ssum += $s[$j];	
				if($x[$j] == $a[1]){
					$usum = $u[$j];
				}
			}
			print "\t",$usum,"\t",$ssum;
		}
		print "\n";

	}
' | stat.edger_interact - $2 $3
}
asp.m3ss_test.test(){
echo \
"chrom	3ss	5ss	num.5ss	strand	c4_3	e1_1	wt_1	wt_3	e1_2	wt_2	e1_3	c4_1	c4_2
chr14	102474514	102473444,102474514	2	+	22,0;6,2	10,0;2,1	7,0;0,1	10,0;0,1	15,0;8,1	25,0;4,0	25,0;10,4	17,0;1,1	31,0;8,1
chr14	105930752	105930483,105930752	2	+	29,0;2,0	20,0;6,3	24,0;2,2	25,0;4,3	35,0;6,6	71,0;12,7	137,0;30,13	12,0;6,3	49,0;2,3"\
| asp.m3ss_test - wt c4

}


asp.group3ss(){
	cat $1 | perl -e 'use strict; my %res=();
	while(<STDIN>){chomp;my@a=split/\t/,$_;
		my $k=$a[0]."\t".$a[6];	
		if($a[5] eq "+"){ 
			$res{$a[0]."\t".$a[5]."\t".($a[2]-1)}{$a[1]}{$a[3]}=$a[4];
		}else{ 
			$res{$a[0]."\t".$a[5]."\t".($a[1])}{$a[2]-1}{$a[3]}=$a[4];
		}
		
	}
	print join("\t",("chrom","start","end","alt5s","sp;unsp","strand")),"\n";
	foreach my $k (keys %res){
		my ($chrom,$strand,$start)=split/\t/,$k;
		my @br=$strand eq "+" ? sort {$b<=>$a} keys %{$res{$k}} : sort {$a<=>$b} keys %{$res{$k}};
		my @ss=(); my @uu=();
		foreach my $x (@br){
			my $s=defined $res{$k}{$x}{"s"}? $res{$k}{$x}{"s"} : 0;
			my $u=defined $res{$k}{$x}{"u"}? $res{$k}{$x}{"u"} : 0;
			push @ss,$s;
			push @uu,$u;
		}	
		if(scalar @ss > 1){
			print join("\t",($chrom, $start, $start+1,
			join(",",@br),join(",",@ss[0..($#ss-1)]).";".join(",",@uu),
			$strand)),"\n";
		}
	}
	'
}

asp.group3ss.test(){
echo \
"chr1	100	200	s	1	+
chr1	50	200	s	11	+
chr1	100	300	s	2	+
chr1	100	101	u	3	+
chr1	199	200	u	4	+
chr1	299	300	u	5	+" \
| asp.group3ss - 
}

asp.jcg(){
cat $1 | perl -e 'use strict; my %res=();
	while(<STDIN>){chomp; my @a=split/\t/,$_;
		my @sizes=split/,/,$a[10];
		my @starts=split/,/,$a[11];
		for(my $i=0; $i<$a[9]-1;$i++){
			my $s=$a[1]+$starts[$i]+$sizes[$i];
			my $e=$a[1]+$starts[$i+1];
			my $l=$sizes[$i];
			my $r=$sizes[$i+1];
			$a[5]=~tr/-+/+-/;	
			if($a[5] eq "-"){ my $tmp=$l;$l=$r;$r=$tmp;}
			push @{$res{$a[0]."\t".$s."\t".$e."\t".$a[5]}{l}},$l;
			push @{$res{$a[0]."\t".$s."\t".$e."\t".$a[5]}{r}},$r;
		}
	}
	foreach my $k (keys %res){
		my @a=split/\t/,$k;
		for(my $i=0; $i< scalar @{$res{$k}{l}}; $i++){
			print $a[2]-$a[1],"\t",$res{$k}{l}->[$i],"\t",$res{$k}{r}->[$i],"\n";
		}
		
#		print join("\t",@a[0..2]),
#			"\t",join(",",@{$res{$k}{l}}),
#			"\t",join(",",@{$res{$k}{r}}),"\t",$a[3],"\n";
	}
'
}
asp.jcg.test(){
echo \
"chr1	1000	2000	a	255	+	1000	2000	0,0,0	2	10,20	0,980
chr1	1001	2001	a	255	+	1001	2001	0,0,0	2	9,21	0,979
chr1	1002	2002	a	255	+	1002	2002	0,0,0	2	8,22	0,978
chr1	1009	1019	a	255	+	1009	1019	0,0,0	1	10	0
chr1	1010	1020	a	255	+	1010	1020	0,0,0	1	10	0
chr1	1011	1021	a	255	+	1011	1021	0,0,0	1	10	0
chr1	1979	1989	a	255	+	1979	1989	0,0,0	1	10	0
chr1	1980	1990	a	255	+	1980	1990	0,0,0	1	10	0
chr1	1981	1991	a	255	+	1981	1991	0,0,0	1	10	0" | asp.jcg -
}
asp.lru(){
	echo '
	my $S="'$1'";
	my @tmp=split/ +/,"'${@:2}'";
	my %F=(); 
	my $i=0;
	foreach my $e (@tmp){
		my @es=split/:/,$e;
		if (scalar @es > 1){ $F{$es[0]}{$es[1]}=1;
		}else{ $F{"tag".$i}{$e}=1; 
		}
		$i++;
	}

	#foreach my $k(keys %F){ foreach my $k1 (keys %{$F{$k}}){ print $k," ",$k1,"\n"; }}

	my @T=keys %F;
	my %L=(); my %R=(); my %U=();
	## preprocessing 
	$i=0;
	foreach my $tag (@T){
		foreach my $f (keys %{$F{$tag}}){
			open my $fh,"<",$f or die "cannot open $f";
			while(<$fh>){ chomp;my @a=split/\t/,$_;
				if( $S eq "a"){ $a[5]="+";
				}elsif($S eq "S"){ $a[5]=~tr/+-/-+/; } 
				my $s=$a[1]; my $e=$a[2]-1; if($a[5] eq "-"){ $s=$a[2]-1;$e=$a[1];}
				my $k=$a[0]."\t".$a[5];
					
				if($a[3] eq "s"){
					$L{$k}{$s}{$e}{$tag}+=$a[4];
					$R{$k}{$e}{$s}{$tag}+=$a[4];
				}elsif($a[3] eq "u"){
					$U{$k}{$s}{$tag}+=$a[4];
				}
			}
			close($fh);
		}
	}
'
}
asp.lru.test(){
asp.lru S a:f1 b:f2
}
asp.m3ss(){
usage(){ echo "
USAGE: $FUNCNAME <tag>:<jc.bed> [<tag>:<jc.bed .. ]
	[options]: 
	 -a : ignore read strands 
	 -s : count read in the same strand 
	 -S : swich read strands 
" 
}
if [ $# -lt 1 ];then usage; return; fi
        local OPTIND; local S="a";
        while getopts ":sSa" arg; do
                case $arg in
                        S) S="S";;
                        s) S="s";;
                        \?) echo "Invalid -${OPTARG}"; return;;
                esac
        done
        shift $(( $OPTIND - 1 ))
	local cmd=`asp.lru $S $@`
	cmd="$cmd"'use strict; 
	print "chrom\t3ss\t5ss\tnum.5ss\tstrand\t",join("\t",@T),"\n";
	foreach my $k (keys %R){
		my ($c,$st)=split/\t/,$k;
	foreach my $s (keys %{$R{$k}}){
		my @es=();
		if($st eq "+"){ @es= sort {$a<=>$b} ($s, keys %{$R{$k}{$s}});
		}else{  @es= sort {$b<=>$a} ($s, keys %{$R{$k}{$s}});}
		print $c,"\t",$s,"\t",join(",",@es),"\t",scalar @es,"\t$st";
		foreach my $t (@T){
			my @vs=(); my @us=();
			foreach my $e (@es){	
				my $v=$R{$k}{$s}{$e}{$t}; $v =0 unless defined $v;
				push @vs,$v;
				my $u=$U{$k}{$e}{$t}; $u =0 unless defined $u;
				push @us,$u;
			}
			print "\t",join(",",@vs);
			print ";",join(",",@us);
		}
		print "\n";
	}}
	'
#	bed.3p $1 #| sort -k1,1 -k2,2n | groupBy -g 1,2,6 -c 3,4,5 -o collapse,collapse,collapse
#	bed.sort $1 | mergeBed -c 2,3,4,5 -s -o collapse,collapse,collapse,collapse \
#	| perl -ne 'chomp;my @a=split/\t/,$_;
#		my $n=join(";",( @a[4..6] ));
#		my $m=$a[7];
#
#		## chr start end starts;ends;jcs;jucs,jdcs n.jcs strand
#		print join("\t",(@a[0..2],$n,$m,$a[3])),"\n";
#	'
	perl -e "$cmd";
}

asp.m3ss.test(){
echo \
"chr1	100	200	s	1	+
chr1	50	200	s	11	+
chr1	100	300	s	2	+
chr1	100	101	u	3	+
chr1	199	200	u	4	+
chr1	299	300	u	5	+" > tmp.a

cat tmp.a
echo "result1"
asp.m3ss tmp.a
echo "result2"
asp.m3ss -S tmp.a
rm tmp.a
}

asp.ce(){
usage(){ echo "
FUNCT: generate cassettee exon events using junction files
USAGE: $FUNCNAME [options] <tag>:<jc> [<tag>:<jc> .. ]
	[options]:
	 -a : ignore strand (default) 
	 -S : switch strand 
	 -s : strand sensitive 
"
}
	
	if [ $# -lt 1 ];then usage; return; fi

        local OPTIND; local S="a";
        while getopts ":sSa" arg; do
                case $arg in
                        S) S="S";;
                        s) S="s";;
                        \?) echo "Invalid -${OPTARG}"; return;;
                esac
        done
        shift $(( $OPTIND - 1 ))
	

	local cmd='use strict; my $S="'$S'"; 
	my @tmp=split/ +/,"'$@'";
	my @files=(); my @tags=();
	my $i=0;
	foreach my $e (@tmp){
		my @es=split/:/,$e;
		if (scalar @es > 1){
			push @files,$es[1];
			push @tags,$es[0];
		}else{
			push @files,$e;
			push @tags,"tag$i";
		}
		$i++;
	}
	
	my %L=(); my %R=(); my %U=();
	## preprocessing 
	$i=0;
	foreach my $f (@files){
		my $tag=$tags[$i++];	
		open my $fh,"<",$f or die "cannot open $f";
		while(<$fh>){ chomp;my @a=split/\t/,$_;
			my $k=$a[0];
			if( $S eq "a"){ $k.=",+";
			}elsif($S eq "s"){ $k.=",".$a[5];   
			}elsif($S eq "S"){ $k.=",".$a[5]; $k=~tr/+-/-+/; } 

			if($a[3] eq "s"){
				$L{$k}{$a[1]}{$a[2]-1}{$tag}+=$a[4];
				$R{$k}{$a[2]-1}{$a[1]}{$tag}+=$a[4];
			}elsif($a[3] eq "u"){
				$U{$k}{$a[1]}{$tag}+=$a[4];
			}
		}
		close($fh);
	}

	print "chrom\tstart\tend\tleft.exon\tright.exon\tstrand\t",join("\t",@tags),"\n";
	## find Cassettee exon candidates
	##    a]---[b      c]---[d 
	foreach my $k (keys %L){
	my ($ch,$st)=split/,/,$k;
	foreach my $a (keys %{$L{$k}}){
	foreach my $d (keys %{$L{$k}{$a}}){ 

		## find the smallest non-overlapping exons
		my @bc=(); my @bc1=();
		foreach my $b (keys %{$L{$k}{$a}}){
		foreach my $c (keys %{$R{$k}{$d}}){ 
			if($b<$c){ push @bc,[$b,$c]; } 
		}}

		for(my $i=0; $i<= $#bc; $i++){
			my $inc=0; 
			#print join("-",@{$bc[$i]}),"\n";
			for(my $j=0; $j<= $#bc; $j++){
				next if( $i==$j);
				## i includes j
				if( $bc[$i]->[0] <= $bc[$j]->[0] && $bc[$i]->[1] >= $bc[$j]->[1]){
					#print join("->>",@{$bc[$j]}),"\n";
					$inc=1; last; 
				}
			}
			if($inc==0){ push @bc1,$bc[$i];}
		}

		## print results
		foreach my $e (@bc1){
			my ($b,$c)=@$e;
			print join("\t",($ch,$b,$c+1,$a,$d+1,$st));
			foreach my $tag (@tags){
				my $ds = defined $L{$k}{$a}{$d}{$tag}? $L{$k}{$a}{$d}{$tag} : 0;
				my $bs = defined $L{$k}{$a}{$b}{$tag}? $L{$k}{$a}{$b}{$tag} : 0;
				my $cs = defined $R{$k}{$d}{$c}{$tag}? $R{$k}{$d}{$c}{$tag} : 0;
				my $abc=join(",",($ds,$bs,$cs));
				my $w=defined $U{$k}{$a}{$tag}? $U{$k}{$a}{$tag} : 0;
				my $x=defined $U{$k}{$b}{$tag}? $U{$k}{$b}{$tag} : 0;
				my $y=defined $U{$k}{$c}{$tag}? $U{$k}{$c}{$tag} : 0;
				my $z=defined $U{$k}{$d}{$tag}? $U{$k}{$d}{$tag} : 0;
				my $wxyz=join(",",($w,$x,$y,$z));
				my $na=$ch.":".$a."-".($d+1).$st;
				my $sc=join(";",($abc,$wxyz));
				print "\t",join("\t",($abc.",".$wxyz));
			}
			print "\n";
		}
	}}}
	'
	perl -e "$cmd"
}

asp.ce.test(){
#    100--120 140-------------200 
#    100--------150 151-------200 
#    100---------------160 
#                        180--200
echo \
"chr1	100	200	s	1	+
chr1	100	150	s	10	+
chr1	151	200	s	11	+
chr1	100	120	s	2	+
chr1	140	200	s	3	+
chr1	100	160	s	4	+
chr1	180	200	s	5	+
chr1	100	101	u	6	-
chr1	119	120	u	7	+" > tmp.a 
 asp.ce -S ctr:tmp.a trt:tmp.a
 asp.ce -a ctr:tmp.a trt:tmp.a
}
asp.ce2ei(){
usage(){ echo "
FUNC: make exclusion and inclusion events of cassette exons
 $FUNCNAME <ce> 
"
}
if [ $# -lt 1 ];then  usage; return; fi
	cat $1 | perl -e 'use strict;
		sub ei{
			my ($x) = @_;
			my @y=split/,/,$x; my ($in,$ex)=(0,0);
			for(my $i=0; $i<=$#y; $i++){
				if($i==1 || $i==2){ $in+=$y[$i];
				}elsif($i==0){ $ex+=$y[$i] }
			}
			return( $ex."\t".$in);
		}	
		my $first=1; my @h=(); my @h2=(); 
		my ($nctr,$ntrt)=(0,0);
		while(<STDIN>){ chomp; my @a=split/\t/,$_; 
			if($first){ 
				@h=@a; $first=0;  
				foreach my $hi (@h[6..$#h]){
					push @h2,$hi.".c1";
					push @h2,$hi.".c2";
				}
				print join("@",@h),"\t",join("\t",@h2),"\n";
				next; 
			}
			print join("@",@a);
			foreach my $ai (@a[6..$#a]){
				print "\t",ei($ai);
			}	
			print "\n";
		}
	'
}

asp.ce2ei.test(){
echo \
"chrom	start	end	left.exon	right.exon	strand	ctr1	ctr2	trt
chr1	159	181	100	200	+	1,4,5,6,0,0,0	1,4,5,6,0,0,0	0,4,5,6,0,0,0
chr1	119	141	100	200	+	1,2,3,6,7,0,0	1,4,5,6,0,0,0	1,2,3,6,7,0,0
chr1	149	152	100	200	+	1,10,11,6,0,0,0	1,4,5,6,0,0,0	1,10,11,6,0,0,0" \
| asp.ce2ei -

}

asp.cecmp(){
usage(){ echo "
FUNC: make exclusion and inclusion events of cassette exons
 $FUNCNAME <ce> <ctr> <trt> 
"
}
	local i=1;
	local inp=;	
	for f in `echo $1 | tr "," " "`;do
		inp+="ctr$i:$f "
		i=$(( $i + 1 ));
	done
	i=1;
	for f in `echo $2 | tr "," " "`;do
		inp+="trt$i:$f "
		i=$(( $i + 1 ));
	done
	asp.ce $inp | perl -e 'use strict;
		sub ei{
			my ($x) = @_;
			my @y=split/,/,$x; my ($in,$ex)=(0,0);
			for(my $i=0; $i<=$#y; $i++){
				if($i==1 || $i==2){ $in+=$y[$i];
				}elsif($i==0 || $i==4 || $i==5){ $ex+=$y[$i] }
			}
			return( $ex."\t".$in);
		}	
		my $first=1; my @h=(); my @h2=(); 
		my ($nctr,$ntrt)=(0,0);
		while(<STDIN>){ chomp; my @a=split/\t/,$_; 
			if($first){ 
				@h=@a; $first=0;  
				foreach my $hi (@h){
					if($hi=~/ctr/){ 
						$nctr++;
						push @h2,"ctr".$nctr.".c1";
						push @h2,"ctr".$nctr.".c2";
					}elsif($hi=~/trt/){
						$ntrt++;
						push @h2,"trt".$ntrt.".c1";
						push @h2,"trt".$ntrt.".c2";
					}
				}
				print join("\t",@h),"\t",join("\t",@h2),"\n";
				next; 
			}
			print $_;
			for(my $i=0; $i<=$#a; $i++){
				if($h[$i]=~/ctr/){
					print "\t",ei($a[$i]);
				
				}elsif($h[$i]=~/trt/){
					print "\t",ei($a[$i]);
				}
			}	
			print "\n";
		}
	'
	
}


asp.ce2ei.test2(){
#    100--120 140-------------200 
#    100--------150 151-------200 
#    100---------------160 
#                        180--200
echo \
"chr1	100	200	s	1	+
chr1	100	150	s	10	+
chr1	151	200	s	11	+
chr1	100	120	s	2	+
chr1	140	200	s	3	+
chr1	100	160	s	4	+
chr1	180	200	s	5	+
chr1	100	101	u	6	+
chr1	119	120	u	7	+" > tmp.a 
 asp.ce2ei tmp.a,tmp.a tmp.a

}

asp.jc(){
usage(){ echo "
Usage : $FUNCNAME <bed12> <out>
Output: <out>.jc 
"
}
if [ $# -lt 1 ];then usage; return; fi
	local tmpd=`mymktempd`; mkdir -p $tmpd/a
	bed.split $1 $tmpd/a
	for f in $tmpd/a/*;do
		bed.intron $f \
		| awk -v OFS="\t" '{ $2=$2-1; $3=$3+1; $4="s"; $5=1; }1' \
		| bed.sum - > $tmpd/${f##*/}.j
		
		awk -v OFS="\t" '{ 
			print $1,$2,$2+1,"l",1,$6;
			print $1,$3-1,$3,"r",1,$6; }' $tmpd/${f##*/}.j > $tmpd/${f##*/}.j53

		awk '$10 == 1' $f | cut -f1-6 \
		| intersectBed -a $tmpd/${f##*/}.j53 -b stdin -wa -wb -s \
		| awk -v OFS="\t" '($4=="l" && $9 > $3 || $4=="r" && $8 < $2){ 
			print $1,$2,$3,"u",$5,$6;
		}' | bed.sum -  > $tmpd/${f##*/}.u 

		cat $tmpd/${f##*/}.j
		cat $tmpd/${f##*/}.u
		
	done
	
	rm -rf $tmpd;
}
asp.jc.test(){
echo \
"chr1	1000	2000	a	255	+	1000	2000	0,0,0	2	10,20	0,980
chr1	1001	2001	a	255	+	1001	2001	0,0,0	2	9,21	0,979
chr1	1002	2002	a	255	+	1002	2002	0,0,0	2	8,22	0,978
chr1	1009	1019	a	255	+	1009	1019	0,0,0	1	10	0
chr1	1010	1020	a	255	+	1010	1020	0,0,0	1	10	0
chr1	1011	1021	a	255	+	1011	1021	0,0,0	1	10	0
chr1	1979	1989	a	255	+	1979	1989	0,0,0	1	10	0
chr1	1980	1990	a	255	+	1980	1990	0,0,0	1	10	0
chr1	1981	1991	a	255	+	1981	1991	0,0,0	1	10	0" | asp.jc -
}

asp.jcprep_ef(){
usage=" $FUNCNAME <method> <bed.jc> 
"
if [ $# -ne 2 ];then echo $usage;return; fi
	cat $2 | perl -e 'use strict; my %res=(); my %u=();  my $M="'$1'";
		while(<STDIN>){chomp; my @a=split/\t/,$_;
			if( $a[3] eq "s"){
				if($a[5] eq "+"){
					$res{5}{$a[0].",".$a[5]}{$a[1]}{$a[2]-1}=$a[4];
					$res{3}{$a[0].",".$a[5]}{$a[2]-1}{$a[1]}=$a[4];
				}else{
					$res{5}{$a[0].",".$a[5]}{$a[2]-1}{$a[1]}=$a[4];
					$res{3}{$a[0].",".$a[5]}{$a[1]}{$a[2]-1}=$a[4];
				}	
			}elsif( $a[3] eq "u"){
				$u{$a[0].",".$a[5]}{$a[1]}=$a[4];
			}
		}
		foreach my $t (keys %res){
		foreach my $c (keys %{$res{$t}}){
		my ($ch,$st)=split/,/,$c;
		foreach my $x (keys %{$res{$t}{$c}}){
			my $ux=defined $u{$c}{$x}? $u{$c}{$x}:0;
			my $sum_s = 0;
			my $sum_u = $ux;
			foreach my $y (keys %{$res{$t}{$c}{$x}}){
				my $uy=defined $u{$c}{$y}? $u{$c}{$y}:0;
				my $s=$res{$t}{$c}{$x}{$y};
				$sum_s += $s;
			}
			if($M eq "3p" && $t == 3){
				print join("\t", ( $c.",".$x,$sum_u,$sum_s)),"\n";
			}
		}}}
	'
}
asp.jcprep(){
usage="$FUNCNAME <method> <tag>:<bed.jc> [<tag>:<bed.jc>]
	<method>: 
		-3p: default 
		-S : switch strand
	 
"
if [ $# -lt 2 ];then echo $usage; return; fi

fnsw(){
	awk -v OFS="\t" '{ if($6=="-"){ $6="+";}else{ $6="-";} print $0;}' $1;
}
	local filt="cat"; # filter function
	local method="3p";
	case $1 in 
		*"-S"* ) filt="fnsw";;
		*"-3p"* ) method="3p";;
	esac

	local tmpd=`mymktempd`;
	local files=( ${@:2} );
	for tf in ${files[@]};do
		if [ -f ${tf#*:} ];then
			local tag=${tf%:*}; local f=${tf#*:} 
			echo "id $tag.c1 $tag.c2" | tr " " "\t" > $tmpd/$tag
			eval $filt $f | asp.jcprep_ef $method - >> $tmpd/$tag
		fi
	done
	## at least one of c2 should be positive
	stat.merge $tmpd/* | run_R '
		tt=read.table("stdin",header=T);
		tt1=tt[apply(tt[,grep(".c2",colnames(tt))],1,function(x){ mean(x > 0)}) >= 0.5,];
		write.table(tt1,"stdout",col.names=T,row.names=F,sep="\t",quote=F);
	'
	rm -rf $tmpd;
}
asp.jcprep.test(){
echo \
"chr1	100	200	s	1	+
chr1	100	300	s	2	+
chr1	100	101	u	3	+
chr1	199	200	u	4	+
chr1	299	300	u	5	+" > tmp.1
asp.jcprep -3p ctr1:tmp.1 ctr2:tmp.1

}

asp.cmpjc(){
usage="
FUNCT: calculate p-value of differential splicing rate
USAGE: $FUNCNAME <ctr1>[,<ctr2> .. ] <trt1>[,<trt2> .. ] 
"
if [ $# -lt 4 ];then echo "$usage"; return; fi
	local i=0; 
	local tmpd=`mymktempd`;
	asp.jcprep $1 $3 $tmpd/ctr 
	asp.jcprep $2 $3 $tmpd/trt 
	fa=`echo $tmpd/ctr* | tr " " ","`
	fb=`echo $tmpd/trt* | tr " " ","`
	if [ $4 == "edger" ];then
		stat.prep $fa $fb | stat.edger - ${5:-""} | perl -ne '$_=~s/\.c2/\.sp/g;$_=~s/\.c1/\.un/g;print $_;'
	else
		echo "ERROR: $4 is unknown method"
	fi
	rm -rf $tmpd;
}

f(){
echo \
"chr1	1000	2000	a	1	+
chr1	1000	2000	b	1	-" > tmp.a

echo \
"chr1	900	1100	r1	100	+
chr1	900	1100	r2	101	-
chr1	1900	2100	r3	255	+
chr1	1900	2100	r4	255	+
chr1	1900	2100	r5	255	-
chr1	1900	2100	r6	255	-
" > tmp.b
}



asp.view(){
usage="     
USAGE: $FUNCNAME <genes.bed12> <reads.bed12>
echo \
"chr1	10004	10148	HISEQ:332:c89cbanxx:1:2301:4302:54544	255	+	10004	10148	255,0,0	2	79,42	0,102
chr1	10076	10256	HISEQ:332:c89cbanxx:1:1214:7419:58717	255	+	10076	10256	255,0,0	2	66,18	0,162
chr1	10094	10465	HISEQ:332:c89cbanxx:1:1306:15197:83203	255	-	10094	10465	255,0,0	2	22,68	0,303
chr1	10105	10467	HISEQ:332:c89cbanxx:1:1309:7231:76774	255	+	10105	10467	255,0,0	2	11,70	0,292
chr1	10105	10465	HISEQ:332:c89cbanxx:1:2109:9823:99561	255	-	10105	10465	255,0,0	2	11,68	0,292
chr1	10108	10460	HISEQ:332:c89cbanxx:1:2111:11378:94238	255	-	10108	10460	255,0,0	2	70,44	0,308
chr1	10110	10269	HISEQ:332:c89cbanxx:1:1202:7536:97863	255	+	10110	10269	255,0,0	2	92,18	0,141
chr1	10110	10260	HISEQ:332:c89cbanxx:1:1203:17552:80952	255	-	10110	10260	255,0,0	2	68,16	0,134
chr1	10111	10463	HISEQ:332:c89cbanxx:1:1215:10190:15343	255	+	10111	10463	255,0,0	2	67,53	0,299
chr1	10111	10440	HISEQ:332:c89cbanxx:1:2115:16272:29726	255	-	10111	10440	255,0,0	2	67,30	0,299" \
| asp.jc -
INPUT: genes (bed12) and reads (bed12)

OUTPUT:
            @
           / \
          *   #
         / \ / \	
    [    ]-[ ]-[   ]

"
	local tmpd=`mymktempd`;
	mycat $1 > $tmpd/a
	mycat $2 > $tmpd/b
	bed.intron  $tmpd/b | bed.flank - 1 1 0 > $tmpd/j
	run_R '
		jc=read.table("'$tmpd/j'",header=F);
		y=as.numeric(jc[,3]-jc[,2]);
		x=jc[,2:3];
		c=apply(x,1,mean);
		xlim=c(min(x[,1]),max(x[,2]));
		ylim=c(0,max(y));
		png("out.png")
		
		plot(NULL,ylim=ylim,xlim=xlim);
		segments(x[,1],0,c,y);
		segments(c,y,x[,2],0);
		
		dev.off();
	' log		

}
asp.view.test(){
echo \
"chr1	1000	2000	g1	0	+	1000	2000	0,0,0	3	100,100,100	0,500,900"> tmp.gene
echo \
"chr1   1050    1550    r1      1       +       1050    1550    255,0,0 2       50,50   0,450
chr1   1050    1950    r1      1       +       1050    1950    255,0,0 2       50,50   0,850"> tmp.read
asp.view tmp.gene tmp.read
rm tmp.read tmp.gene
}

asp.plot_scatter(){
usage="$FUNCNAME [options] <spi|3ss> <spi|3ss>
 [options]:
        -s <int> : minimum counts
        -o <str.png> : output png file (default: out.png)
"
        local thre=100;
        local fpng="out.png";
        local OPTIND;
        while getopts ":s:o:" arg; do
                case $arg in
                        s) thre=${OPTARG};;
                        o) fpng=${OPTARG};;
                        \?) echo "Invalid -${OPTARG}"; return;;
                esac
        done
        shift $(( $OPTIND - 1 ))
        if [ $# -ne 2 ];then echo "$usage"; return; fi

        local tmpd=`mymktempd`;
        for f in $@;do
                local n=${f##*/};
                bed.enc $f | cut -f 4,7- | tr "," "\t" \
                | awk -v thre=$thre -v OFS="\t" '$2 + $3 + $4 > thre{ print $1,$5;}'\
                > $tmpd/$n
        done

        myjoin -d ";" $tmpd/* \
        | run_R '
        library(MASS);
        library(ggplot2);
        tt=read.table("stdin",header=T);
        xlab=colnames(tt)[2];
        ylab=colnames(tt)[3];
        print(ylab);
        tt=tt[apply(tt[,-1],1,function(x){ sum(is.na(x))==0;}),];
        tt[ tt < 0]=0;
        x=tt[,2];y=tt[,3];
        DF= data.frame(x,y);
        dens = kde2d(x,y);
        gr = data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z));
        names(gr) = c("xgr", "ygr", "zgr");
        mod = loess(zgr~xgr*ygr, data=gr);
        DF$pointdens = predict(mod, newdata=data.frame(xgr=x, ygr=y));
        p=ggplot(DF, aes(x=x,y=y, color=pointdens)) + geom_point() +  xlab(xlab) + ylab(ylab);
	p=p + geom_hline(aes(yintercept=median(y)))
	p=p + geom_vline(aes(xintercept=median(x)))
        png("'$fpng'");
        print(p);
        dev.off();
        '
        rm -rf $tmpd;
}


asp.prep(){
usage="
FUNCT: convert bed6+ to id+ form
USAGE: $FUNCNAME <bed6+>[,<bed6+>] <evt> <out>
"
if [ $# -lt 2 ];then echo "$usage";return; fi
	local i=0; 
	for f in `echo $1 | tr "," " "`;do
		i=$(( $i + 1 ));
		if [ $2 == "3ss" ];then
			bed.enc $f | awk -v OFS="\t" '{ split($7,a,","); print $4,a[1],a[2];}' > $3.$i;
		elif [ $2 == "spi" ];then
			bed.enc $f | awk -v OFS="\t" '{ split($7,a,","); print $4,a[2]+a[3],a[1];} ' > $3.$i;
		fi
	done
}
asp.prep.test(){
echo \
"chr1	1	2	n1	0	+	1,3,1	0.6" > tmp.3ss.ctr1
echo \
"chr1	1	2	n1	0	+	1,3,1	0.6" > tmp.3ss.ctr2
echo \
"chr1	1	2	n1	0	+	3,1,1	0.3" > tmp.3ss.trt1
echo \
"chr1	1	2	n1	0	+	3,2,1	0.2" > tmp.3ss.trt2
asp.prep tmp.3ss.trt1,tmp.3ss.trt2 3ss tmp.3ss.trt
asp.prep tmp.3ss.ctr1,tmp.3ss.ctr2 3ss tmp.3ss.ctr
head tmp.3ss*
rm tmp.*	
}

asp.cmp(){
usage="
FUNCT: calculate p-value of differential splicing rate
USAGE: $FUNCNAME <ctr1>[,<ctr2> .. ] <trt1>[,<trt2> .. ] <input_type> <method> [<bcd>]
	<input_type>: [3ss|spi|..]
	<method>: edger
"
if [ $# -lt 4 ];then echo "$usage"; return; fi
	local i=0; 
	local tmpd=`mymktempd`;
	asp.prep $1 $3 $tmpd/ctr 
	asp.prep $2 $3 $tmpd/trt 
	fa=`echo $tmpd/ctr* | tr " " ","`
	fb=`echo $tmpd/trt* | tr " " ","`
	if [ $4 == "edger" ];then
		stat.prep $fa $fb | stat.edger - ${5:-""} | perl -ne '$_=~s/\.c2/\.sp/g;$_=~s/\.c1/\.un/g;print $_;'
	else
		echo "ERROR: $4 is unknown method"
	fi
	rm -rf $tmpd;
}

asp.cmp.test(){
echo \
"chr1	1	2	n1	0	+	1,3,1	0.6" > tmp.3ss.ctr1
echo \
"chr1	1	2	n1	0	+	1,3,1	0.6" > tmp.3ss.ctr2
echo \
"chr1	1	2	n1	0	+	3,1,1	0.3" > tmp.3ss.trt1
echo \
"chr1	1	2	n1	0	+	3,2,1	0.2" > tmp.3ss.trt2
asp.cmp tmp.3ss.ctr1,tmp.3ss.ctr2 tmp.3ss.trt1,tmp.3ss.trt2 3ss edger

echo \
"chr1	1	2	n1	0	+	1,1,3	0.6" > tmp.spi.ctr1
echo \
"chr1	1	2	n1	0	+	1,1,3	0.6" > tmp.spi.ctr2
echo \
"chr1	1	2	n1	0	+	3,3,1	0.3" > tmp.spi.trt1
echo \
"chr1	1	2	n1	0	+	3,3,2	0.2" > tmp.spi.trt2
asp.cmp tmp.spi.ctr1,tmp.spi.ctr2 tmp.spi.trt1,tmp.spi.trt2 spi edger
rm tmp.*
}

asp.spi(){
	usage(){
	echo " 
	 FUNCT: Calculate SPI
	 USAGE: $FUNCNAME [options] <intron.bed6> <read.bed12> 

		[options]:
		 -s : only count reads on the same strand of intron.bed6 
		 -S : only count reads on the same opposite strand of intron.bed6 
		 -c <flag> : 
			0 : use as it is (default)
			1 : count as 1
			2 : convert phred score to probability

	 OUTPUT: bed6 a,b,c a/(0.5*(b+c)+a)
			       (a)	
			_/     sp       \_
		[    5'  ]--------------[   3'      ]
		       _(b)_          _(c)_
	"
	return;
	}
	if [ $# -lt 2 ];then echo "$usage"; return; fi
	local OPTIND opt strand count C S
	C=0; S="";
	while getopts ":sSc:" opt; do
		case "${opt}" in
			s) S="-s";;
			S) S="-S";;
			c) C=$OPTARG;;
			*) usage
		esac
	done
	shift $((OPTIND-1))
	
	local tmpd=`mymktempd`; 
	touch $tmpd/a $tmpd/b $tmpd/c $tmpd/r.1 $tmpd/r.2
	#local tmpd=tmpd; rm -rf $tmpd; mkdir -p $tmpd
	mycat $1 | bed.enc - | sort -u > $tmpd/i
	bed.5p $tmpd/i > $tmpd/i.5
	bed.3p $tmpd/i > $tmpd/i.3
	mycat $2 | bed.score - $C \
		| awk -v OFS="\t" -v O=$tmpd/r '{ n=1; if($10 > 1){ n=2;} print $0 >> O"."n; }'

	## count spliced
	bed.intron $tmpd/r.2 \
		| intersectBed -a $tmpd/i -b stdin -wa -wb -f 1 -F 1 $S \
		| cut -f4,11 | stat.sum - | sort -k 1,1  > $tmpd/a 

	## count unspliced
	intersectBed -a $tmpd/i.5 -b $tmpd/r.1 -wa -wb $S \
		| cut -f4,11 | stat.sum - | sort -k 1,1 > $tmpd/b

	intersectBed -a $tmpd/i.3 -b $tmpd/r.1 -wa -wb $S \
		| cut -f4,11 | stat.sum - | sort -k 1,1 > $tmpd/c
	
	join -a 1 -a 2 -e 0 -o 0,1.2,2.2 $tmpd/a $tmpd/b \
	| join -a 1 -a 2 -e0 -o 0,1.2,1.3,2.2 - $tmpd/c \
	| awk '{ print $1,$2","$3","$4,$2/(0.5*($3+$4)+$2);}' \
	| tr "@ " "\t"	

	rm -rf $tmpd;
}
asp.spi.test(){
echo \
"c	100	200	intron1	0	+
c	100	199	intron2	0	+" > tmp.intron
echo \
"c	90	210	r1	10	+	90	210	0,0,0	2	10,10	0,110
c	189	199	r2	1	+	190	200	0,0,0	1	10	0
c	190	200	r3	1	+	190	200	0,0,0	1	10	0
c	191	201	r4	1	+	190	200	0,0,0	1	10	0
c	199	209	r5	1	+	190	200	0,0,0	1	10	0
c	200	210	r6	1	+	190	200	0,0,0	1	10	0
c	201	211	r7	1	+	190	200	0,0,0	1	10	0" > tmp.r
head tmp.*
echo "==> result"
	asp.spi tmp.intron tmp.r 
	rm tmp.intron tmp.r
}

asp.3ss_smart(){ 
usage=" 
 FUNCT: calculate 3SS statistics for each intron
 OUTPU: intron.bed6 a,b,c  1-a/(b+c)
 USAGE: $FUNCNAME <intron.bed6> <read.bed12> <window> <strand> [<count>]
	<window>: windowsize for a and b boundaries 
	<strand>: 0, 1(same strand), 2(opposite strand) 
	<count> : 0: use as it is (default), 1: count as 1, 2: phred score

# >>>>> : read
# >>>>>----->>>>> : splicing read
# [     ] : exon
# ------- : intron
# |-----| : predefined window (e.g., 25bp)
# _/   \_ : junction counts


               _/     (c)       \_
 [  5' exon    ]-----------------[   3' exon     ]
                         |--(a)--|--(b)--|

                     >>>>>      >>>>>           (support a)
                                 >>>>>   >>>>>  (support b)
          >>>>>>----------------->>>>>>         (support c)
             >>>>>>-------------->>>>>>         (support none)
          >>>>>>--------------------->>>>>>     (support none)

"
	if [ $# -lt 4 ];then echo "$usage"; return; fi
	local W=$3;
	local S=${4/0/""}; S=${S/1/"-s"}; S=${S/2/"-S"};

	local tmpd=`mymktempd`;
	touch $tmpd/a $tmpd/b $tmpd/c $tmpd/r.1 $tmpd/r.2

	mycat $1 | bed.enc - | sort -u > $tmpd/i 
	bed.3p $tmpd/i | bed.flank - $(( $W -1 )) 0 1  > $tmpd/i.a
	bed.3p $tmpd/i | bed.flank - -1 $W 1 > $tmpd/i.b 
	mycat $2 | bed.score - ${5:-"0"} \
		| awk -v OFS="\t" -v O=$tmpd/r '{ n=1; if($10 > 1){ n=2;} print $0 >> O"."n; }'
	intersectBed -a $tmpd/i.a -b $tmpd/r.1 $S -wa -wb \
		| cut -f4,11 | stat.sum - | sort -k1,1 > $tmpd/a 	
	intersectBed -a $tmpd/i.b -b $tmpd/r.1 $S -wa -wb \
		| awk '$6=="+" && $8>=$2 || $6=="-" && $9<=$3' \
		| cut -f4,11 | stat.sum - | sort -k1,1 > $tmpd/b 	
	bed.intron $tmpd/r.2 \
		| intersectBed -a $tmpd/i -b stdin $S -wa -wb -f 1 -F 1 \
		| cut -f 4,11 | stat.sum - | sort -k1,1 > $tmpd/c

	join -a 1 -a 2 -e 0 -o 0,1.2,2.2 $tmpd/a $tmpd/b \
	| join -a 1 -a 2 -e 0 -o 0,1.2,1.3,2.2 - $tmpd/c \
	| awk '{v=$3+$4; if(v==0){ v="-inf";}else{ v=1 - $2/v;}  print $1,$2","$3","$4,v;}' \
	| tr "@ " "\t"
	
	rm -rf $tmpd
}
asp.3ss(){ 

usage(){ 
echo "
 FUNCT: calculate 3SS statistics for each 3 prime regions 
 OUTPU: 3\' regions of introns + a,b,c + 1-a/b (note, c is added in b)
 USAGE: $FUNCNAME [options] <intron.bed6> <read.bed12> 

  [options]:
	-w <int> : a window size for upstream (unspliced) and downstream (spliced) regions 
	-s : only count reads on the same strand of intron.bed6 (default none)
	-S : only count reads on the same opposite strand of intron.bed6 (default none)
	-c <int> : 
		0 : use as it is (default)
		1 : count as 1
		2 : convert phred score to probability

               _/  _/  (c)       \_
 [  5\' exon    ]-----------------[   3\' exon     ]
                          |--(a)--|--(b)--|
" 
return; 
}

	if [ $# -lt 2 ];then usage; return; fi
	local OPTIND opt C S W
	C=0; S=""; W=25;
	while getopts ":sSc:w:" opt; do
		case "${opt}" in
			s) S="-s";;
			S) S="-S";;
			c) C=$OPTARG;;
			w) W=$OPTARG;;
			*) usage
		esac
	done
	shift $((OPTIND-1))

	local tmpd=`mymktempd`;
	touch $tmpd/a $tmpd/b $tmpd/c $tmpd/r.1 $tmpd/r.2

	mycat $1 | bed.3p - | sort -u | bed.enc - > $tmpd/i 

	bed.3p $tmpd/i | bed.flank - $(( $W -1 )) 0 1  > $tmpd/i.a
	bed.3p $tmpd/i | bed.flank - -1 $W 1 > $tmpd/i.b 

	mycat $2 | bed.score - $C \
		| awk -v OFS="\t" -v O=$tmpd/r '{ n=1; if($10 > 1){ n=2;} print $0 >> O"."n; }'

	cut -f1-6 $tmpd/r.1 > $tmpd/r.s
	bed.exon $tmpd/r.2 >> $tmpd/r.s

	intersectBed -a $tmpd/i.a -b $tmpd/r.s $S -wa -wb \
		| cut -f4,11 | stat.sum - | sort -k1,1 > $tmpd/a 	
	intersectBed -a $tmpd/i.b -b $tmpd/r.s $S -wa -wb \
		| cut -f4,11 | stat.sum - | sort -k1,1 > $tmpd/b 	
	bed.intron $tmpd/r.2 \
		| intersectBed -a $tmpd/i -b stdin $S -wa -wb \
		| awk '$6=="+" && $3==$9 || $6=="-" && $2==$8' \
		| cut -f 4,11 | stat.sum - | sort -k1,1 > $tmpd/c

	join -a 1 -a 2 -e 0 -o 0,1.2,2.2 $tmpd/a $tmpd/b \
	| join -a 1 -a 2 -e 0 -o 0,1.2,1.3,2.2 - $tmpd/c \
	| awk '{v=$3; if(v==0){ v="-inf";}else{ v=1 - $2/v;}  print $1,$2","$3","$4,v;}' \
	| tr "@ " "\t"
	
	rm -rf $tmpd
}
asp.3ss.test(){ 
echo \
"chr1	100	200	intron1	0	+
chr1	100	199	intron2	0	+" > tmp.intron
echo \
"chr1	90	210	r1	1	+	90	210	0,0,0	2	10,10	0,110
chr1	80	210	r2	1	+	80	210	0,0,0	2	10,10	0,120
chr1	189	199	r3	1	+	189	199	0,0,0	1	10	0
chr1	190	200	r4	1	+	190	200	0,0,0	1	10	0
chr1	191	201	r5	1	+	191	201	0,0,0	1	10	0
chr1	199	209	r6	1	+	199	209	0,0,0	1	10	0
chr1	200	210	r7	1	+	200	210	0,0,0	1	10	0
chr1	201	211	r8	1	+	201	211	0,0,0	1	10	0" > tmp.r
head tmp.*
echo "==> result"
	asp.3ss -s tmp.intron tmp.r 
	rm tmp.intron tmp.r
}
asp.3ss_lessfp(){ 
usage="
FUNNC: calculate 3ss with less false positives caused by:
	1. introns w/o splicing reads 
	2. wrong 5/3 \' junctions
	3. overlapping reads at the 3\' end 

USAGE: $FUNCNAME <intron.bed6> <read.bed12> <window> <strand>
	<window>: windowsize for a and b boundaries 
	<strand>: 0, 1(same strand), 2(opposite strand) 

                     |---a---|---b---|
	[    ]---------------[             ]
	   \                 /       ( wrong 5\' splicing )
	     \                   /   ( wrong 3\' splicing )
	OUTPUT: intron.bed6 + a + b 
"
	if [ $# -lt 4 ];then echo "$usage"; return; fi
	local W=$3;
	local S=${4/1/"-s"}; S=${S/2/"-S"};

	tmpd=`mymktempd`;
	#tmpd=tmpd; mkdir -p $tmpd; 
	mycat $1 | awk -v OFS="\t" '{ $4=$1"@"$2"@"$3"@"$4"@"$5"@"$6;} 1' > $tmpd/a
	mycat $2 > $tmpd/b
	awk '$10==1' $tmpd/b | cut -f1-6 > $tmpd/b.1
	
	## junctio counts
	bed.intron $tmpd/b | intersectBed -a $tmpd/a -b stdin -f 1 -F 1 $S -wa -wb \
	| awk -v OFS="\t" '{ print $4,$11;}' | stat.sum - > $tmpd/c.j

	cat $tmpd/c.j | tr "@" "\t" | awk -v OFS="\t" '{ $4=$1"@"$2"@"$3"@"$4"@"$5"@"$6;} 1' \
	| bed.3p - | bed.flank - $(( $W -1))  $W 1 \
	| intersectBed -a stdin -b $tmpd/b.1 -wa -wb $S \
	| perl -e 'use strict; my $W='$W'; my %res=();
		while(<STDIN>){ chomp; my @a=split/\t/,$_;
			## bed7 + bed6
			my ($s, $e) = ($a[1],$a[2]);
			my ($rs, $re) = ($a[8],$a[9]);

			my $x="x";
			if($s + $W < $rs ){ $x="b";		
			}elsif( $s + $W >= $re ){ $x="a"; }
			if($a[5] eq "-"){
				if($x eq "a"){ $x = "b";}
				elsif($x eq "b" ){ $x = "a";}
			}	
			$res{$a[3]}{"j"} = $a[6];
			$res{$a[3]}{$x} += $a[11];
		}
		foreach my $k (keys %res){
			my $a=defined $res{$k}{"a"} ? $res{$k}{"a"} : 0;
			my $b=defined $res{$k}{"b"} ? $res{$k}{"b"} : 0;
			my $c=defined $res{$k}{"j"} ? $res{$k}{"j"} : 0;
			print join("\t", split /@/,$k),"\t",$a,"\t",$b+$c,"\n";
		}
	'
	rm -rf $tmpd
}
aso.cosi(){
	## BED: 2-3: intron boundary, 5th: exon coordinates
	##      /    a     \         /   b     \ 
	##               __c__     __d__       
	##      /               e              \
	##     ]------------[       ]-----------[
    BED=$1; BAM=$2;
    chroms=( `cut -f1 $BED | sort -u` )
    TMP=`make_temp`
    for CHROM in ${chroms[@]}
    do
        echo " $CHROM .." >&2
        awk -v CHROM=$CHROM '$1==CHROM' $BED > $TMP
        samtools view -bq 255 $BAM $CHROM | bamToBed -bed12 \
		| intersectBed -a $TMP -b stdin -wa -wb \
		| awk -v OFS="\t" '{ split($17,l,",");split($18,s,",");split($5,ex,",");
			is=$2-1;ie=$3; es=ex[1];ee=ex[2]-1; 
			a=0;b=0;c=0;d=0;e=0;
			for(i=1;i<=$16;i++){ ## count unsplicing events
					rs=$8+s[i]; re=rs+l[i]-1;
					if(rs < es && re > es){ c=1;
					}else if(rs < ee && re > ee){ d=1;}
			}
			for(i=2;i<=$16;i++){ ## count splicing events
					rs=$8+s[i-1]+l[i-1]-1; re=$8+s[i];
					if(rs==is && re== es){ a=1;
					}else if(rs==ee && re == ie){ b=1;
					}else if(rs ==is && re == ie){ e=1;}
			}
			if(a+b+c+d+e > 0){
				print $1,$2,$3,$4,$5,$6,a,b,c,d,e;
			}
		}'  | groupBy -g 1,2,3,4,5,6 -c 7,8,9,10,11 -o sum,sum,sum,sum,sum
	done
}

#
#gen_bed(){
#	echo "hi" | awk -v OFS="\t" -v s=$1 -v e=$2 -v l=$3 '{
#		for(i=s; i< e; i+=1){
#			print "c",i,i+l,"p",1,"+";	
#		}	
#		for(i=s; i< e; i+=1){
#			print "c",i,i+l,"n",1,"-";	
#		}	
#	}'
#
#}
#test__3ss(){
#echo \
#"c	100	200	intron1	0	+
#c	100	200	intron1	0	-" > tmpa
#echo \
#"c	100	200	intron1	0	+	25	25
#c	100	200	intron1	0	-	25	25" > exp
#
#gen_bed 0 250 50 | 3ss tmpa - 25 -Sm > obs
#check obs exp
#rm -rf tmpa exp obs
#}
##test__3ss
#
#test__count_us(){
#fig="
#Count unsplicing for exons:
#            100         200
#	    [          ]    : exon
#       ======      ======   : unsplicing
#            ===== ======    : not unsplicing
#"
#echo \
#"chr1	100	200	e	0	+" > a
#echo \
#"chr1	50	101	r1	1	+
#chr1	150	201	r2	10	+
#chr1	100	150	r3	100	+
#chr1	150	200	r4	1000	+" > b
#echo \
#"chr1	100	200	e	0	+	11" > exp
#count_us a b > obs
#check exp obs
#rm -f a b exp obs
#}
##test__count_us
#
#bed12_to_jc(){
#	bed12_to_intron $1 \
#	| awk -v OFS="\t" '{ $2=$2-1;$3=$3+1;} 1' \
#	| bed_sum - 
#}
#bede_to_bed6(){
### copy 5th,7th to the end to the name field, set 5th to be 0
#	cat $1 | perl -ne 'chomp; my @a=split/\t/,$_;
#	print join("\t",@a[0..3]),"@",$a[4],"@",join("@",@a[6..$#a]),"\t0\t",$a[5],"\n"; '
#}
#bed6_to_bede(){
#	cat $1 | perl -ne 'chomp; my @a=split/\t/,$_; my @b=split/@/,$a[3]; $a[3] = shift @b; $a[4] = shift @b; 
#	print join("\t",@a),"\t",join("\t",@b),"\n";	'
#}
#
#filter_intron_skip(){
#	awk '$2-1==$8 && $3+1==$9'
#}
#
#
#count_is(){
#	opt=${3:-""};
#	intersectBed -a ${1/-/stdin} -b ${2/-/stdin} -wa -wb $opt \
#	| awk -v OFS="\t" '$2==$8+1 && $3==$9-1{ print $1,$2,$3,$4"|"$5,$11,$6;}' \
#	| bed_sum - \
#	| awk -v OFS="\t" '{ split($4,a,"|");print $1,$2,$3,a[1],a[2],$6,$5;}' 
#}
#test__count_is(){
#echo "c	100	200	.	1000	+" > a
#echo "c	99	201	.	1	+
#c	99	200	.	10	+
#c	99	201	.	2	-" >b
#	count_is a b -s > obs
#	count_is a b -S >> obs
#	count_is a b >> obs
#echo "c	100	200	.	1000	+	1
#c	100	200	.	1000	+	2
#c	100	200	.	1000	+	3" > exp
#cat obs
#	check exp obs
#	rm -f obs exp a b
#
#}
##test__count_is;
#
#
#
#count_jei(){
#usage="
#usage: $FUNCNAME <target.bed6> <read.bed12>
#function: count exclusion and inclusion of junction reads
#<read.bed12> : use modify_bed .. to manage score the column
#
#
#output: bed6@ count_exc count_inc
#       /                  \     exclusion
#                 /         \    exclusion
#            \          /        inclusion
#  -----------[        ]----------
#"
#if [ $# -lt 2 ]; then echo "$usage"; return; fi
#        bed12_to_junction $2 \
#        | intersectBed -a $1 -b stdin -wa -wb \
#        | awk -v OFS="@" '{
#                e=0;i=0; ## exclusion and inclusion
#                if( $8 == $3 -1  || $9 == $2 + 1){ i=$11;}
#                else{ e=$11;}
#                print $1,$2,$3,$4,$5,$6"\t"e"\t"i;
#        }' | sort -k1,1 | groupBy -g 1 -c 2,3 -o sum,sum
#}
#
#
#make_temp(){
#	mktemp 2>/dev/null || mktemp -t $0
#}
#
#
#fdr_thre(){
#	F=$1; T=$3;
#	cmd='
#		tt = read.table("stdin",header=F);
#		x= tt[,C];
#		ix = p.adjust(x,method="fdr") <= T;
#		cat(paste("res=",max(x[ix]),"\n",sep=""));
#	'
#	cmd=${cmd/C/$C};
#	cmd=${cmd/T/$T};
#
#	tmp=`make_temp`
#	echo "$cmd" > $tmp
#	cat $F | R --no-save -f $tmp | perl -ne 'chomp;if($_=~/res=([\d|\.]+)/){ print $1,"\n";}'
#}
#
#count_a53ss_jc(){
#	##          /         a         \        
#	##               /    b         \
#	##     [   |   ]----------------[       ]
#    BED=$1; JBED=$2;
#	intersectBed -a $BED -b $JBED -wa -wb -f 1 -r \
#	| awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$10; }'
#}
#
#count_a53ss(){
#	##          /         a         \        
#	##               /    b         \
#	##     [   |   ]----------------[       ]
#    BED=$1; BAM=$2;
#    chroms=( `cut -f1 $BED | sort -u` )
#    TMP=`make_temp`
#    for CHROM in ${chroms[@]}
#    do
#        echo " $CHROM .." >&2
#        awk -v CHROM=$CHROM '$1==CHROM' $BED > $TMP
#        samtools view -bq 255 $BAM $CHROM | bamToBed -bed12 \
#		| intersectBed -a $TMP -b stdin -wa -wb \
#		| awk -v OFS="\t" '{ alen=split($7,a,",");blen=split($8,b,",");split($19,sizes,",");split($20,starts,",");
#			for(i=1 ; i<= alen; i++){ ai=a[i]-1;
#			for(j=1 ; j<= blen; j++){ bj=b[j]+1;
#			for(k=2 ; k<= $18; k++){
#				s=$10 + starts[k-1] + sizes[k-1]-1;		
#				e=$10 + starts[k] + 1;		
#				if(ai == s && bj == e ){
#					print $1,$2,$3,$4,ai "," bj,$6,1;
#				}
#			}}}
#		}' | sort -k1 -k2,3n | groupBy -g 1,2,3,4,5,6 -c 7 -o sum | awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7;}'
#	done
#}
##head -n 100 a53ss.bed  > tmp
##count_a53ss tmp ../Tophat/Wt1/accepted_hits.bam -wa -wb
#
#
#count_bed(){
#    BED=$1; BAM=$2;
#    chroms=( `cut -f1 $BED | sort -u` )
#    TMP=`make_temp`
#    for CHROM in ${chroms[@]}
#    do
#        echo " $CHROM .." >&2
#        awk -v CHROM=$CHROM '$1==CHROM' $BED > $TMP
#        samtools view -bq 255 $BAM $CHROM | bamToBed -bed12 \
#		| intersectBed -a $TMP -b stdin -wa -wb \
#		| awk -v OFS="\t" '{ split($17,l,",");split($18,s,","); 
#			hit=0;
#			if($16==1){ hit=1;
#			}else{
#				for(i=1;i<=$16;i++){
#					if(i==1 && $8+s[i]+l[i] == $3){  ##  [  _]/
#						hit=1;
#					}else if(i==$16 && $8+s[i] == $2){ ##  \[_  ]
#						hit=1;
#					}else if($2 == $8+s[i] && $3==$8+s[i]+l[i]){ ## \[__]/
#						hit=1;
#					}
#				}
#			}
#			if(hit){ print $1,$2,$3,$4,$5,$6,1;}
#		}'  | groupBy -g 1,2,3,4,5,6 -c 7 -o sum
#	done
#}
#
#count_nj(){
#    BED=$1; BAM=$2;
#    chroms=( `cut -f1 $BED | sort -u` )
#    TMP=`make_temp`
#    for CHROM in ${chroms[@]}
#    do
#        echo " $CHROM .." >&2
#        awk -v CHROM=$CHROM '$1==CHROM' $BED > $TMP
#        samtools view -bq 255 $BAM $CHROM | bamToBed -split \
#        | intersectBed -a $TMP -b stdin -wa -wb \
#        | awk '{ if($2 > $8 || $3 < $9){ print $1,$2,$3,$4,$5,$6;}}' \
#        | uniq -c | awk -v FS=" " -v OFS="\t" '{ print $2,$3,$4,$5,$6,$7,$1;}'
#    done
#}
#
#intron_excinc(){
#    NJ=$1; JC=$2;
#    intersectBed -a $NJ -b $JC -wa -wb  \
#    | awk -v OFS="\t" '{ if($2-1==$9 && $3+1==$10){ print $1,$2,$3,$4,$5,$6,$11,$7;}}'
#}
#
#
#make_retainable_introns(){
### [    ]—-intron--[   ]
### [                   ]
#	intersectBed -a $1 -b $2 -wa -wb \
#	| awk -v OFS="\t" '{ split($5,a,",");
#		if($8 == a[1] && $9 == a[2]){ ## match boundary
#			print $1,$2,$3,$4,0,$6;
#		}
#	}' | sort -u  
#	## unique retained intron candidates
#}
#
#
#
#
#
#
#
### changed from bed12_to_exon
#bed12_to_exonEvents(){
#	## [s1  ]e1----[s   ]e----[s2   ]e2
#	awk -v OFS="\t" '{
#		## take introns
#		split($11,sizes,",");
#		split($12,starts,",");
#		for(i=1;i<= $10;i++){
#			s = $2 + starts[i];
#			e = $2 + starts[i]+sizes[i];	
#			s1 = -1; e1 = -1;
#			s2 = -1; e2 = -1;
#			if( i > 1){ s1 = $2 + starts[i-1]; e1 = s1 + sizes[i-1]; }
#			if( i < $10){ s2 = $2 + starts[i+1]; e2 = s2 + sizes[i+1]; }
#			split($4,gene,"::");
#			## gene base
#			print $1,s,e,gene[1],s1 "," e1 "," s2 "," e2,$6;
#		}	
#	}' | sort -u
#}
#
#
#
#
#each_chrom(){
#	## lambda function accepts chrom and size parameters
#	## lambda(){ 
#	##   # .. handle $1 $2 
#	## } 
#	chrom_size_file=$1; lambda_func=$2;
#	tmp=( `cat $chrom_size_file` )
#	for (( i=0; i< ${#tmp[@]}; i+=2 ))
#	do
#		chrom=${tmp[$i]}; size=${tmp[$i+1]};
#		$lambda_func $chrom $size
#	done
#}
#
#
#bam_to_junctions(){
#	## split bam by chrom to reduce memory usage
#	bam=$1; chromsize=$2; quality=$3;
#	lambda(){
#		samtools view -bq $quality $bam $1 | bamToBed -bed12 | awk '$10 > 1' | bed12_to_junction
#	}
#	each_chrom $chromsize lambda
#}
##bam_to_junctions $1 $2 255
#exon_excinc_junction(){
#    ## input: exon and junction_count
#    ## output: exon skipping and junction counts 
#	## [   ]2----[5.1  5.2]---3[   ]
#    EXON=$1;JUNCTION=$2;
#	cat $EXON | awk -v OFS="\t" '{ 
#		## take neighbor intron boundary
#   		split($5,a,","); ## [a1  ]a2--[2   ]3--[a3  ]a4
#		s=$2;e=$3;
#		if(a[2] > 0){ s= a[2];} # intron start
#		if(a[3] > 0){ e= a[3];} # intron end
#		print $1,s, e, $4, $2 "," $3,$6; 
#	}' | sort -u | intersectBed -a stdin -b $JUNCTION -wa -wb  | awk -v OFS="\t" '{
#    	## [  ]e1----[e2   e3]---e4[   ]
#    	split($5,a,","); e1=$2-1; e2=a[1]; e3=a[2]-1; e4=$3; # exon
#		js = $8; je = $9-1; jc=$10; # junction 
#        x=0;y=0;z=0;
#        if(js==e1 && je==e2 ){ ## left splicing  
#            x=jc;
#        }else if(js==e3 && je==e4){ ## right splicing 
#            y=jc;
#        }else if(js==e1 && je==e4){ ## [  ]/  [  ] \[  ]skipping
#            z=jc;
#        }
#    	if( x+y+z > 0){ 
#			## output: bed6: [   2]--[5.1   5.2]---[3   ], counts:  7,8,9
#			#print $1,$2,$3,$4,$5,$6,x,z,y; ## left splicing, skipping, right splicing
#			print $1,$2,$3,$4,$5,$6,z,x+y;  ## exclusion, inclusion
#		}
#    }' | sort -k1,1 -k2,3n -k4,6 | groupBy -g 1,2,3,4,5,6 -c 7,8 -o sum,sum 
## 	| awk -v OFS="\t" '{
##    	split($5,a,","); # [   2]---[a1  a2]---[3  ]
##		if($8>0 && $7> 0){
##			print $1,$2,$3,$2,a[1],$6,$8,$7;
##		}
##		if($8>0 && $9>0){
##			print $1,$2,$3,a[2],$3,$6,$8,$9;
##		}
##	}' | sort -u
#}
#
#count_nonjunction_events(){
#	#   --  --  ---  : left unspliced, within, right unspliced
#    #   -----------  : within = 0, left unspliced == right unspliced 
#	#    [       ]
#	bed6=$1; bam=$2; chromsize=$3; quality=$4
#	lambda(){
#		samtools view -bq $quality $bam $1| bamToBed -bed12 \
#		| intersectBed -a $bed6 -b stdin -wa -wb \
#		| awk -v OFS="\t" '{ 
#			L = 0; C = 0; R = 0;
#			if($8 < $2 && $9-1 > $2){ L += 1;}  # left unspliced
#			if($8 < $3-1 && $9 > $3){ R += 1;}  # right unspliced 
#			if($8 > $2 && $9 < $3){ C += 1;}    # within 
#			if( L + C + R > 0){
#				print $1,$2,$3,$4,$5,$6,L,C,R;
#			}
#		}' | sort -k1,1 -k2,3n -k4,6 | groupBy -g 1,2,3,4,5,6 -c 7,8,9 -o sum,sum,sum  
#	}
#	each_chrom $chromsize lambda
#}
##count_nonjunction_events $1 $2 $3 $4
#
#count_intron_junction_events(){
#	intron=$1; junction=$2;
#	intersectBed -a $intron -b $junction -wa -wb  | awk -v OFS="\t" '{
#		##     /       \
#		## [   ]-------[    ]
#		if($2-1 == $8 && $3+1 == $9){ 
#			print $1,$2,$3,$4,$5,$6,$10;
#		}
#	}'
#}
#		
#count_intron_events(){
#	intron=$1;bam=$2;junction=$3;chromsize=$4;quality=$5;
#	tmp1=`make_temp`
#	tmp2=`make_temp`
#	count_intron_junction_events  $intron $junction  > $tmp1
#	count_nonjunction_events $intron $bam $chromsize $quality > $tmp2
#	intersectBed -a $tmp2 -b $tmp1 -wa -wb -f 1 -r -s | awk -v OFS="\t" '{
#		if($4 == $13){ ## shared by different genes
#			## exclusion and inclusion
#			print $1,$2,$3,$4,$5,$6,$16,$7+$9;	
#		}
#	}'
#}
##count_intron_events Events/rintrons.bed ../Tophat/Wt1/accepted_hits.bam Events/Wt1/jc.bed Data/chrom.size 255
#
#quote(){ 
#	perl -ne 'chomp; my @a = map{ "\"$_\"" } split /,/,$_; print join ",",@a; '; 
#}
#
#testI(){
#	## INPUT: comma separated control and treatment EI file ( bed6 + exclusion + inclusion)
#	## OUTPUT: bed6 + logFC + pvalue 
#	a=`echo $1 | quote`
#	b=`echo $2 | quote`
#	rcmd='
#	fa=c(FILEA)
#	fb=c(FILEB)
#	#fa=c("Events/Wt1/exons_jc.bed","Events/Wt2/exons_jc.bed")	
#	#fb=c("Events/C41/exons_jc.bed","Events/C42/exons_jc.bed")	
#	out="OUT"
#	group=factor(c(rep(1,length(fa)),rep(2,length(fb))));
#	D=NULL;
#	i=1;
#	for( f in c(fa,fb)){
#		tt=read.table(f,header=F);
#		colnames(tt)=c("chr","start","end","name","score","strand",paste(i,c("exc","inc"),sep="."))
#		if(is.null(D)){ D=tt;
#		}else{ D=merge(D,tt,by=1:6,all=T); }
#		i=i+1;
#	}
#	D[is.na(D)]=0;
#	ix = apply(D[,grep("inc",colnames(D))],1, min) > 2;
#	d=D[ix,1:6];
#	m=D[ix,grep("inc",colnames(D))];
#	o=order(d$name);
#	## make exclusion event from the gene sum
#	d=d[o,]; m=m[o,]; s=apply(m,2,function(x){ ave(x,d$name,FUN=sum)});
#	s = s-m; colnames(s)=paste(1:ncol(s),"exc");
#	m = cbind(s, m);
#
#
#	## test
#    library(edgeR)
#    #j=ncol(m)/2; #y=DGEList(counts=m[,1:j]+m[,(j+1):ncol(m)],group=group)
#    y=DGEList(counts=m,group=rep(group,2))
#    y=calcNormFactors(y);
#
#    event.this=factor(rep(1:2,each=length(group)));
#    group.this=factor(rep(group,2));
#    H1 <- model.matrix(~ event.this + group.this + event.this:group.this, data=y$samples )
#    H0 <- model.matrix(~ event.this + group.this )
#    coef <- (ncol(H0)+1):ncol(H1)
#    #y=estimateCommonDisp(y)
#    #y=estimateTagwiseDisp(y, trend="movingave")
#    y = estimateGLMCommonDisp(y,H1);
#    y = estimateGLMTrendedDisp(y,H1);
#    y = estimateGLMTagwiseDisp(y,H1);
#
#    fit=glmFit(m, H1, y$tagwise.dispersion,offset=0,prior.count=0)
#    llh=glmLRT(fit,coef=coef)
#
#    ex.h0=apply( m[,group.this == 1 & event.this == 1], 1, sum);
#    in.h0=apply( m[,group.this == 1 & event.this == 2], 1, sum);
#
#    res=data.frame(d[,1:6], logIR=log( in.h0/ ex.h0), logFC=llh$table$logFC, pval=llh$table$PValue)
#    ## chrom start end logFC pval
#    write.table(res, out, col.names=T,row.names=F,sep="\t",quote=F);
#	'
#	tmp=`mktemp`
#	rcmd=${rcmd/FILEA/$a}
#	rcmd=${rcmd/FILEB/$b}
#	rcmd=${rcmd/OUT/$tmp}
#	echo "$rcmd" | R --no-save >&2
#	cat $tmp
#}
#
#test_exin_fisher(){
#	## input : bed6 + exclusion + inclusion counts
#	intersectBed -a $1 -b $2 -wa -wb -f 1 -r -s \
#	| awk -v OFS="@" '{ print $1,$2,$3,$4,$5,$6"\t"$7"\t"$8"\t"$15"\t"$16;}' \
#	| fisher_test -
#}
#test_exin_edger(){
#	## INPUT: comma separated control and treatment EI file ( bed6 + exclusion + inclusion)
#	## OUTPUT: bed6 + logFC + pvalue 
#	a=`echo $1 | quote`
#	b=`echo $2 | quote`
#
#	rcmd='
#	fa=c(FILEA)
#	fb=c(FILEB)
#	#fa=c("Events/Wt1/exons_jc.bed","Events/Wt2/exons_jc.bed")	
#	#fb=c("Events/C41/exons_jc.bed","Events/C42/exons_jc.bed")	
#	out="OUT"
#
#	group=factor(c(rep(1,length(fa)),rep(2,length(fb))));
#	D=NULL;
#	i=1;
#	for( f in c(fa,fb)){
#		tt=read.table(f,header=F);
#		colnames(tt)=c("chr","start","end","name","score","strand",paste(i,c("exc","inc"),sep="."))
#		if(is.null(D)){ D=tt;
#		}else{ D=merge(D,tt,by=1:6,all=T); }
#		i=i+1;
#	}
#	D[is.na(D)]=0;
#	ix=apply(D[,7:ncol(D)], 1, min) > 0 & apply(D[,7:ncol(D)],1,max) > 10
#	d=D[ix,]
#	#d[d==0]=0.5
#	m=cbind(d[,grep("exc",colnames(d))],d[,grep("inc",colnames(d))])
#
#	library(edgeR)
#	#j=ncol(m)/2; #y=DGEList(counts=m[,1:j]+m[,(j+1):ncol(m)],group=group)
#	y=DGEList(counts=m,group=rep(group,2))
#	y=calcNormFactors(y);
#
#	event.this=factor(rep(1:2,each=length(group)));
#	group.this=factor(rep(group,2));
#	H1 <- model.matrix(~ event.this + group.this + event.this:group.this, data=y$samples )
#	H0 <- model.matrix(~ event.this + group.this )
#	coef <- (ncol(H0)+1):ncol(H1)
#	#y=estimateCommonDisp(y)
#	#y=estimateTagwiseDisp(y, trend="movingave")
#	y = estimateGLMCommonDisp(y,H1);
#	y = estimateGLMTrendedDisp(y,H1);
#	y = estimateGLMTagwiseDisp(y,H1);
#
#	fit=glmFit(m, H1, y$tagwise.dispersion,offset=0,prior.count=0)
#	llh=glmLRT(fit,coef=coef)
#
#	ex.h0=apply( m[,group.this == 1 & event.this == 1], 1, sum);
#	in.h0=apply( m[,group.this == 1 & event.this == 2], 1, sum);
#
#	res=data.frame(d[,1:6], logIR=log( in.h0/ ex.h0), logFC=llh$table$logFC, pval=llh$table$PValue)
#	## chrom start end logFC pval
#	write.table(res, out, col.names=T,row.names=F,sep="\t",quote=F);
#	'
#
#	tmp=`mktemp`
#	rcmd=${rcmd/FILEA/$a}
#	rcmd=${rcmd/FILEB/$b}
#	rcmd=${rcmd/OUT/$tmp}
#	echo "$rcmd" | R --no-save >&2
#	cat $tmp
#}
#
##data='chr7	99647389	99662661	ENSG00000166529::ENST00000543588	0	+	99647389	99662661	0,0,0	6	75,133,435,193,192,956	0,1711,7204,7931,14021,14316
##chr7	99647396	99662661	ENSG00000166529::ENST00000456748	0	+	99647396	99662661	0,0,0	5	68,495,193,192,956	0,7137,7924,14014,14309
##chr7	99654523	99662660	ENSG00000166529::ENST00000379635	0	+	99654523	99662660	0,0,0	4	349,81,193,1250	0,424,797,6887
##chr7	99647416	99655464	ENSG00000166529::ENST00000438937	0	+	99647416	99655464	0,0,0	4	48,133,495,144	0,1684,7117,7904
##chr7	99654533	99661878	ENSG00000166529::ENST00000477297	0	+	99654533	99661878	0,0,0	3	495,193,173	0,787,7172
##chr7	99647396	99662659	ENSG00000166529::ENST00000292450	0	+	99647396	99662659	0,0,0	4	68,495,193,1249	0,7137,7924,14014';
#
#bed12_to_exon(){
#	awk -v OFS="\t" '{
#		## take introns
#		split($11,sizes,","); split($12,starts,",");
#		for(i=1;i<= $10;i++){ ## 1-base index not perl
#			ord=i; if($5 eq "-"){ ord = $10-ord+1; }
#			s = $2 + starts[i]; e = $2 + starts[i]+sizes[i];	
#			print $1,s,e,$4,0,$6;
#		}	
#	}' $1
#}
#
#bed12_to_a53ss(){
#	## event id is composed of start,end,gene
#	## 
#	##     /  splicing event    \
#	## [a   b]--------------------[ c   ]d : numbers represent columns
#	## event_id:  gene:a-d:c or gene:a-d:b
#	awk -v OFS="\t" '{ split($11,sizes,",");split($12,starts,",");split($4,gt,"::");
#		for(i=2;i<=$10;i++){
#			a=$2+starts[i-1]; b=a+sizes[i-1]-1;
#			c=$2+starts[i]; d=c+sizes[i];
#			print $1 "," gt[1] ":" a "-" d ":" b "," $6,b,c+1;
#			print $1 "," gt[1] ":" a "-" d ":" c "," $6,b,c+1;
#		}
#	}'  | sort -uk 1 | groupBy -g 1 -c 2,3 -o collapse,collapse \
#	| perl -ne ' chomp; my @a = split /\t/,$_;
#		my @starts = split /,/,$a[1];
#		next if scalar @starts < 2;
#		my ($chrom,$eid,$strand) = split /,/,$a[0];
#		my ($event_gene, $event_range, $event_pivot) = split /:/,$eid;
#		my @ends = split /,/,$a[2];
#		if($starts[0] == $event_pivot){ ## [   ]-----[  ][  ]
#			my $type = $strand eq "+" ? "A3" : "A5";
#			my @ix = sort { $ends[$a] cmp $ends[$b] } 0 .. $#ends;
#			for my $i (0 .. $#ix){
#				print join "\t",($chrom,$starts[$ix[$i]],$ends[$ix[$i]],$eid,$type.":".$i,$strand),"\n";	
#			}
#		}else{ ## [   ][   ]-------[   ]
#			my $type = $strand eq "+" ? "A5" : "A3";
#			my @ix = sort { $starts[$b] cmp $starts[$a] } 0 .. $#starts;
#			for my $i (0 .. $#ix){
#				print join "\t",($chrom,$starts[$ix[$i]],$ends[$ix[$i]],$eid,$type.":".$i,$strand),"\n";	
#			}
#		} 
#	'
#}
#
#
#bed12_to_rintrons(){
#	local tmp=`make_temp`	
#	local tmp1=`make_temp`
#	local tmp2=`make_temp`
#	cat > $tmp
#	cat $tmp | bed12_to_exon > $tmp1
#	cat $tmp | bed12_to_dumbbells > $tmp2
#	intersectBed -a $tmp2 -b $tmp1 -r -f 1  -s -wa 
#}
#
#bed12_to_a53ss1(){
## 5' splicing      [$2                 a1]----[a2 $3]
##                  [    ]                     [     ]
## 3' splicing      [$2 a1]---------[a2           $3 ]
#	bed12_to_2exons | groupBy -g 1,2,3,4,6 -c 5,5 -o collapse,count | awk '$7 > 1' #\
##	 | intersectBed -a stdin -b $tmp1 -wa -f 1 -r -s -v \
##	 | perl -ne ' chomp; my ($c,$s,$e,$n,$st,$tmp) = split /\t/,$_;
##		my %h5=();
##		my %h3=();
##		my %ri_filter=();
##		foreach my $e (split /,/,$tmp){
##			my ($es, $ee) = split /:/,$e;
##			$h5{$ee}{$es} = 5; # 5 prim splicing [   |   ]------[   ]
##			$h3{$es}{$ee} = 3; # 3 prim splicing [   ]--------[  |  ]
##		}	
##		foreach $k (keys %h5){
##			my @a=keys %{$h5{$k}};
##			if(scalar @a > 1){
##				my $st1 = ($st eq "+" ? 5 : 3);
##				print "$c\t$s\t$e\t$n\t$st1\t$st\t";
##				print join( ",",@a),"\t",$k,"\n";
##			}
##		}
##		foreach $k (keys %h3){
##			my @a=keys %{$h3{$k}};
##			if(scalar @a > 1){
##				my $st1 = ($st eq "+" ?  3 : 5);
##				print "$c\t$s\t$e\t$n\t$st1\t$st\t";
##				print $k,"\t",join( ",",@a),"\n";
##			}
##		}
##	' 
#} 
### test
##echo -e "$data" | bed12_to_a53ss
#
#
