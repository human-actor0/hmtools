. $HMHOME/src/bed.sh;


splicing.toy(){
usage="$FUNCNAME <intput.txt>"
if [ $# -lt 1 ];then echo "$usage"; return; fi
cat $1 | perl -ne 'use strict; chomp; 
	my @S=split//," ".$_." "; ## add pads for eacy calc.
	my @sizes=(); my @starts=();
	my $type="";
	my $start=0;
	for(my $i=1; $i< scalar @S; $i++){ 
		my $a=$S[$i-1];
		my $b=$S[$i];
		my $p=$i-1; ## original position before the padding
		if( ($a eq " " || $a eq "-") && $b =~ /[eErR]/){
			$type=$b;
			if( $a eq " "){ $start=$p;}
			push @starts,$p-$start;
		}elsif($a =~ /[eErR]/ && ($b eq " " || $b eq "-")){
			push @sizes, $p - $starts[$#starts] - $start;
		}
	}
	if( $type ne ""){
		my $strand = $type eq uc $type ? "+" : "-";
		my $end = $start + $starts[$#starts] + $sizes[$#sizes];
		print join("\t",(
			"chr1",
			$start, $end,
			$type, 0, $strand,
			$start, $end,
			"0,0,0",scalar @starts,
			join(",",@sizes),join(",",@starts)
		)),"\n";
	}
'
}

splicing.toy.test(){
echo \
"01234567890123456789012345678901234567890123456789
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
 EEEE-------------EEEEEEEE------------EEEEEEEEEEEE
  RRR-------------RR
   rr-------------rrrrrr
                RRRR 
                                    RRRRR
" | splicing.toy -
}

splicing.count_3SS(){
usage="
FUNCT: count splicing events occuring at 3 prime splicing sites
USAGE: $FUNCNAME <intron3p.bed> <read.bed12> <window> <min_splicing_dist> [options]
 [options]:
	-s : count reads on the same strand 
	-S : count reads on the opposite strand 
"

if [ $# -lt 4 ];then echo "$usage"; return; fi
        intersectBed -a ${1/-/stdin} -b ${2/-/stdin} -wa -wb ${@:5} \
        | perl -e 'use strict; my %res=(); my $W='$3'; my $M='$4';
	sub same3{
		my ($s1,$e1,$s2,$e2,$strand) = @_;
		my $res=0;
		if($strand eq "+" && $e1==$e2 || $strand eq "-" && $s1==$s2){ $res=1;}
		print join(",",@_),"=>",$res,"\n";
		return $res;
	}
	sub within{
		my ($s1,$e1,$s2,$e2) = @_;
		if( $s1 >= $s2 && $e1 <= $e2){ return 1;}
		return 0;
	}
	sub cross{
		my ($s1,$e1,$s2,$e2) = @_;
		if( $s1 < $s2 && $e1 > $e2){ return 1;}
		return 0;
	}
        while(<STDIN>){ chomp; my @a=split/\t/,$_;
                my $tstart=$a[1]; my $tend=$a[2]; my $tstrand=$a[5]; 
                my @sizes=split/,/,$a[16]; my @starts=split/,/,$a[17];
                ## handle crossing reads
                my $us=0; ## unsplicing events
		my @tags=();
            #    for( my $i=0; $i< $a[15]; $i++){
            #            my $s1=$a[7]+$starts[$i];
            #            my $e1=$s1+$sizes[$i];
	    #    	my $tag="";
	    #    	print join(" ",($s1,$e1,$tstart,$tend)),"\n";
	    #    	if( same3($s1,$e1,$tstart,$tend)){ 	
	    #    		if($e1-$s1 > $W){ $tag="S";
	    #    		}else{$tag="U";}
	    #    	}
	    #    	if($tag ne ""){
	    #    		$res{ join("\t",@a[0..5]) }{$tag} ++;
	    #    	}
            #    }
                for( my $i=1; $i< $a[15]; $i++){ ### extranc intronic to count splicing reads 
                        my $s1=$a[7]+$starts[$i-1];
                        my $e1=$s1+$sizes[$i-1];
                        my $s2=$a[7]+$starts[$i];
                        my $e2=$s2+$sizes[$i];
			if( $e2-$s1 > $M &&
				same3($e1,$s2,$tstart,$tend,$tstrand)){ 
				push @tags,"S"; 
			}elsif( cross($s1,$e1,$tstart,$tend) ){
				push @tags,"U"; 
			}
                }
		foreach my $t (@tags){
	    		$res{ join("\t",@a[0..5]) }{$t} ++;
		}
        }
        foreach my $k (keys %res){
		print $k,"\n";
		foreach my $tag (keys %{$res{$k}}){
			print $tag," ",$res{$k}{$tag},"\n";
		}
               # my $U=defined $res{$k}{U} ? $res{$k}{U} : 0;
               # my $S=defined $res{$k}{S} ? $res{$k}{S} : 0;
               # my $u=defined $res{$k}{u} ? $res{$k}{u} : 0;
               # my $s=defined $res{$k}{s} ? $res{$k}{s} : 0;
               # print join("\t",($S,$U,$s,$u)),"\n";
        }
        '
}

splicing.count_3SS.test(){
echo \
"chr1	150	210	U	0	+	150	210	0,0,0	2	10,10	0,50
chr1	180	210	S	0	+	180	210	0,0,0	2	10,10	0,20
chr1	100	210	r1	0	+	100	210	0,0,0	2	10,10	0,100
chr1	0	210	r2	0	+	0	210	0,0,0	2	10,10	0,200
chr1	191	201	r3	0	+	191	201	0,0,0	1	10	0
chr1	189	199	r3	0	+	189	199	0,0,0	1	10	0
chr1	192	202	r4	0	+	192	202	0,0,0	1	10	0
chr1	199	209	r5	0	+	199	209	0,0,0	1	10	0
chr1	200	210	r6	0	+	200	210	0,0,0	1	10	0
chr1	300	310	r7	0	+	300	310	0,0,0	1	10	0" > tmp.read

echo \
"chr1	199	200	intron3p	0	+
chr1	110	111	intron3p	0	-" > tmp.intron3p

splicing.count_3SS tmp.intron3p tmp.read 50 20

}

splicing.junction(){
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
splicing.junction.test(){
echo \
"chr1	1000	2000	a	255	+	1000	2000	0,0,0	2	10,20	0,980
chr1	1001	2001	a	255	+	1001	2001	0,0,0	2	9,21	0,979
chr1	1002	2002	a	255	+	1002	2002	0,0,0	2	8,22	0,978
chr1	1009	1019	a	255	+	1009	1019	0,0,0	1	10	0
chr1	1010	1020	a	255	+	1010	1020	0,0,0	1	10	0
chr1	1011	1021	a	255	+	1011	1021	0,0,0	1	10	0
chr1	1979	1989	a	255	+	1979	1989	0,0,0	1	10	0
chr1	1980	1990	a	255	+	1980	1990	0,0,0	1	10	0
chr1	1981	1991	a	255	+	1981	1991	0,0,0	1	10	0" | splicing.junction  -
}



junction.count(){
usage(){ echo "
Usage : $FUNCNAME <junction.bed> <read.bed12>  [options]
 [options]:
	-s : strand specific count
	-S : count overlapping reads with an opposite strand 
"
}
if [ $# -lt 2 ];then usage; return; fi
	bed.exon $2 \
	| intersectBed -a $1 -b stdin  -wa -wb ${3:""} \
	| perl -e 'use strict; my %res=();
	while(<STDIN>){chomp; my @a=split/\t/,$_;
		if( $a[7] < $a[1] ){
			$res{ join("\t",@a[0..5]) }{"L"} ++;
		}elsif( $a[8] > $a[2] ){
			$res{ join("\t",@a[0..5]) }{"R"} ++;
		}else{
			$res{ join("\t",@a[0..5]) }{"W"} ++;
		}
	}
	foreach my $k (keys %res){
		my $l= defined $res{$k}{"L"} ? $res{$k}{"L"} : 0;	
		my $w= defined $res{$k}{"W"} ? $res{$k}{"W"} : 0; 
		my $r= defined $res{$k}{"R"} ? $res{$k}{"R"} : 0; 
		print $k,"\t",join(",",( $l,$w,$r)),"\n";
	}
	' 
}

junction.count.test(){
echo \
"chr1	100	200	j1	1	+
chr1	50	200	j2	1	-" > tmp.a
echo \
"chr1	99	101	r1	1	+	100	200	0,0,0	1	2	0	
chr1	99	101	r5	1	+	100	200	0,0,0	1	2	0	
chr1	199	201	r2	1	-	199	200	0	1	2	0
chr1	100	101	r3	1	+	100	200	0,0,0	1	2	0	
chr1	199	200	r4	1	-	199	200	0	1	2	0
" > tmp.b
junction.count tmp.a tmp.b
#rm tmp.a tmp.b
	
}

splicing.junction_to_SUtable(){
usage="$FUNCNAME <ctr.junction>[,<ctr.junction..] <trt.junction>[,<trt.junction>] [options]
 [options]
	-s : strand specific (default)
	-S : switch strand 
 	-a : ignore strand (convert - to +)
"
perl -e 'use strict;
	my $option="'${3:-"-s"}'";
	my @fctr=map{"$_"} split /,/,"'$1'";
	my @ftrt=map{"$_"} split /,/,"'$2'";

	sub readf{
		my ($f,$r,$tag,$opt)=@_;
		my %D=(); my %U=();
		open(my $fh,"<",$f) or die "$! : $f";
		while(<$fh>){chomp; my @a=split/\t/,$_;
			my $v=$a[4]; $a[4]=0;
			$a[5]="+" if $opt eq "-a";
			$a[5]=~tr/+-/-+/ if $opt eq "-S";
			if($a[3] eq "s"){
				my $k=join("\t",@a);
				$r->{$k}->{$tag.".c2"} += $v;
				$D{$a[1]}=$k; $D{$a[2]-1}=$k; ## link to edges
			}elsif($a[3] eq "u"){
				$U{$a[1]} += $v;
			}
		}	
		close($fh);
		foreach my $k (keys %U){
			if(defined $D{$k}){
				$r->{$D{$k}}->{$tag.".c1"} += $U{$k};
			}
		}
	}
	my %res=();
	my $i=0;
	my @cols=();
	foreach my $f (@fctr){
		readf($f,\%res,"ctr_".$i,$option);
		push @cols,"ctr_".$i.".c1";
		push @cols,"ctr_".$i.".c2";
		$i++;
	}
	$i=0;
	foreach my $f (@ftrt){
		readf($f,\%res,"trt_".$i,$option);
		push @cols,"trt_".$i.".c1";
		push @cols,"trt_".$i.".c2";
		$i++;
	}
	print join("\t",("chrom","start","end","name","score","strand")),"\t";
	print join("\t",@cols),"\n";
	foreach my $k (keys %res){
		my $output=$k;
		my $nz=0;
		foreach my $c (@cols){
			my $v=defined $res{$k}{$c} ? $res{$k}{$c} : 0;
			if($v > 0){ $nz++;}
			$output .= "\t".$v;
		}
		if ( $nz > (scalar @cols)/3 ){ ### filtering zero-inflated rows
			print $output,"\n";
		}
	}
'

}

splicing.junction_to_SUtable.test(){
	echo \
"chr1	100	200	s	1	+
chr1	100	101	u	2	+
chr1	50	200	s	1	-
chr1	200	201	u	10	-
chr1	199	200	u	3	+" > tmp.a

splicing.junction_to_SUtable tmp.a,tmp.a tmp.a,tmp.a -a

}
splicing.junction_to_table(){
usage="$FUNCNAME <ctr.junction>[,<ctr.junction..] <trt.junction>[,<trt.junction>] [options]
 [options]
	-s : strand specific (default)
	-S : switch strand 
 	-a : ignore strand (convert - to +)
"
perl -e 'use strict;
	my $option="'${3:-"-s"}'";
	my @fctr=map{"$_"} split /,/,"'$1'";
	my @ftrt=map{"$_"} split /,/,"'$2'";

	sub readf{
		my ($f,$r,$tag,$opt)=@_;
		open(my $fh,"<",$f) or die "$! : $f";
		while(<$fh>){chomp; my @a=split/\t/,$_;
			my $v=$a[4];
			$a[4]=0;
			$a[5]="+" if $opt eq "-a";
			$a[5]=~tr/+-/-+/ if $opt eq "-S";
			$r->{join("\t",@a)}{$tag}=$v;
		}	
		close($fh);
	}
	my %res=();
	my $i=0;
	my @cols=();
	foreach my $f (@fctr){
		readf($f,\%res,"ctr_".$i,$option);
		push @cols,"ctr_".$i;
		$i++;
	}
	$i=0;
	foreach my $f (@ftrt){
		readf($f,\%res,"trt_".$i,$option);
		push @cols,"trt_".$i;
		$i++;
	}
	print join("\t",("chrom","start","end","name","score","strand")),"\t";
	print join("\t",@cols),"\n";
	foreach my $k (keys %res){
		print $k;
		foreach my $c (@cols){
			my $v=0; $v=$res{$k}{$c} if defined $res{$k}{$c};
			print "\t",$v;
		}
		print "\n";
	}
'

}
