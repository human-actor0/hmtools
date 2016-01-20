. $HMHOME/src/root.sh

<<<<<<< HEAD



mhits.em(){
## EM algorithm is referred to
# https://github.com/mhammell-laboratory/tetoolkit/blob/master/TEToolkit/EMAlgorithm.py

usage="
FUNCT: redistribute multi-hits using EM algorithm
USAGE: $FUNCNAME <target> <reads> <strand_opt> <norm_opt> [options] 
	<target>    : bed6 w/ unique names at 4th column
	<reads>     : bed6 w/ <uniq_id>.<multihits> format (produced by bed.n2i)
	<strand_opt>: [0: all, 1:count reads in the same strand, 2: count .. opposite strand]
	<norm_opt>  : [0: no addtional normalization, 1: weight by effective length]
"; if [ $# -lt 4 ];then echo "$usage"; return; fi
local MAXITER=100;
local DELTA="1.e-07";
local STRAND_OPT=$3; STRAND_OPT=${STRAND_OPT/0/}; STRAND_OPT=${STRAND_OPT/1/\-s}; STRAND_OPT=${STRAND_OPT/2/\-S};
local NORM_OPT=$4;

	intersectBed -a ${1/\-/stdin} -b $2 -wa -wb $STRAND_OPT \
	| perl -e 'use strict; 
		my $NORM_OPT='$NORM_OPT';
		my $MAXITER='$MAXITER';
		my $OPT_TOL='$DELTA';

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
				if ($NORM_OPT==1){
					my $len = $L->{$g} - $Lread + 1;
					if($len >0){ $E->{$g} /= $len;
					}else{ $E->{$g} = 0;}
				}
				$s += $E->{$g};
			}
			my $norm=0; if( $s > 0){ $norm=1/$s;}
			## normalize gene expressions
			foreach my $g (keys %$E){ $E->{$g} *= $norm;}
		}
		sub updateE{
			my ($Ein, $A,      $Eout) = @_;
 			foreach my $k (keys %$Eout) { delete $Eout->{$k}; }

			## handle multi hits independent of unique hits
			foreach my $r (keys %$A){
				my @G=keys %{$A->{$r}};
				my $s=0; my $norm=0;
				## multi hits onnly
				if ( scalar @G < 2){ next;}
				foreach my $g (@G){ $s += $Ein->{$g}; }
				if($s > 0){ $norm = 1/$s; }
				foreach my $g (@G){ $Eout->{$g} += $Ein->{$g} * $norm; }
			}

			## handle unique hits
			foreach my $r (keys %$A){
				my @G=keys %{$A->{$r}};
				if( scalar @G == 1){
					foreach my $g(@G){ $Eout->{$g} ++; }
				}
			}
		}
		sub update{
			my ($Ein, $A,$L,$Lread,     $Eout) = @_;

			## estimate expression by reads
			updateE( $Ein, $A, $Eout);
			normE( $Eout, $L, $Lread);
		}

		my $minStep0 = 1.0;
		my $minStep = 1.0;
		my $maxStep0 = 1.0;
		my $maxStep = 1.0;
		my $mStep = 4.0;


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
		my %r=(); my %r2=(); my %v=();
		foreach my $iter ( 1..$MAXITER){
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
				foreach my $k (keys %E2){ $E{$k}=$E2{$k}; }
			}
		
			## extrapolate 
			foreach my $g (keys %E){
				$Eprime{$g} = max(0.0, $E{$g} + 2 * $alphaS * $r{$g} + ($alphaS*$alphaS)*$v{$g});
			}	
			update(\%Eprime,\%A,\%L,$Lread,\%E);

			## stabilization
			if (abs( $alphaS - 1.0 ) > 0.01 ){
				print {*STDERR} "stablizaation step: iteration=$iter\n";
				if ($alphaS == $maxStep){ 
					$maxStep=max($maxStep0, $maxStep/$mStep);
					$alphaS = 1.0;
				}
			}
			if ($alphaS == $maxStep){ $maxStep = $mStep * $maxStep;}
			if ($minStep < 0 && $alphaS == $minStep ){ $minStep = $mStep * $minStep;}

		}
		updateE( \%E, \%A, \%E1);
		foreach my $k (keys %E1){ print $k,"\t",$E1{$k},"\n";}
	'
}
mhits.em.test(){
echo \
"c	0	100	t1	0	+
c	100	250	t2	0	+
c	300	500	t3	0	+" > tmp.genes

echo \
"c	10	21	r1	1	+
c	10	21	r2	1	+
c	110	120	r2	1	+
c	310	320	r2	1	+" \
> tmp.reads

mhits.em tmp.genes tmp.reads 1 0
rm tmp.genes tmp.reads
}

mhits.cgraph(){
#obtain connected graph
	cat $1 | perl -e 'use strict;
		my %tmp=();
		my %U=(); # unvisited
		my %V=(); # visited
		my %N=();

		while(<STDIN>){ chomp; 
			my ($t,$r) = split/\t/,$_;
			$tmp{$r}{$t}=1;
			$U{$t}=1;
		}
		## make a network
		foreach my $r (keys %tmp){
			my @neig=keys %{$tmp{$r}};
			for( my $i=0; $i< $#neig  ; $i++){
			for( my $j=$i+1; $j <= $#neig; $j++){
				$N{$neig[$i]}{$neig[$j]}=1;
			}}
		}
		sub df(){
			my ($U,$G,$C) = @_;
			foreach my $v (keys %{$G->{$u}}){
				df(	
				delete $U->{$v};	
			}	
		}
		
				
	';
}

mhits.cgraph.test(){
=======
mhits.clique(){
	cat $1 | perl -e 'use strict;
		my %res=();
		while(<STDIN>){ chomp; 
			my ($t,$r) = split/\t/,$_;
			$res{$r}{$t}=1;
		}
		foreach my @$r (keys %res){
			if($#r > 0){
			}
		}
	';
}

mhits.clique.test(){
>>>>>>> 25773133df552fe5b6d4ecae6e912ef7f7fcbcf5
echo "t1	r1
t1	r2
t2	r1
t2	r2
<<<<<<< HEAD
t3	r1" | mhits.cgraph -
=======
t3	r1" | mhits.clique -
>>>>>>> 25773133df552fe5b6d4ecae6e912ef7f7fcbcf5
}
