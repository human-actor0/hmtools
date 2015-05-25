#!/bin/bash
THIS=${BASH_SOURCE##*/};
W=0;
usage="
$THIS <543> [<543> ..]
 	<543> : id, group, feature, count
	[options]:
	-w <int>: smoothing <int> neighbors (default $W);
"
while getopts "w:" arg; do
	case $arg in
		w) W=${OPTARG};;
	esac
done
shift $(( OPTIND -1 ));
if [ $# -lt 1 ];then
	echo "$usage"; exit 1;
fi
cmd='use strict; my %data=(); my %cols=(); my $window=W;
	for(my $i=0; $i <= $#ARGV; $i++){
		open(F,$ARGV[$i]) or die "$ARGV[$i]";	
		while(<F>){ chomp;
			my ($id,$group,$fea,$count) = split /\t/,$_;
			my $cn="$group.$fea";
			$cols{$cn}=1;
			$data{$id}{$i}{$cn} = $count;
		}
		close(F);
	}
	sub mysort(a,b){
		my ($a1,$a2) = split /\./,$a;
		my ($b1,$b2) = split /\./,$b;
		if($a1 eq $b1){ 
			return $a2 <=> $b2;
		}
		return $b1 <=> $a1;
	}
	
	my @cols = sort mysort keys %cols;
	print "id\tsample\t",join( "\t",@cols),"\n";
	foreach my $id (keys %data){    
	for(my $i=0; $i <= $#ARGV; $i++){
		print $id,"\t",$i;
		for(my $j=0; $j<=$#cols; $j++){
			my $c=$cols[$j];
			my $v=0;
			$v = $data{$id}{$i}{$c} if defined $data{$id}{$i} && defined $data{$id}{$i}{$c};
			if($window > 0){
				my ($l,$r)=($j-$window,$j+$window);
				$l=0 if($l < 0);
				$r=$#cols if ($r > $#cols);
				for(my $k=$l; $k<=$r; $k++){
					$v += $data{$id}{$i}{$cols[$k]} if($k != $j && defined $data{$id}{$i} && defined $data{$id}{$i}{$cols[$k]});
				}
			}
			print "\t$v";
		}
		print "\n";
	}}
'
cmd=${cmd//W/$W};
perl -e "$cmd"  $@

