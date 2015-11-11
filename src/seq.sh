#!/bin/bash
. $HMHOME/src/root.sh

kmers(){
## obtained from Cook malcolm
## http://comments.gmane.org/gmane.comp.lang.perl.bio.general/18242
k=$1;
s=$( printf "%${k}s" ); # a string with $k blanks
s=${s// /{A,T,G,C\}};   # substitute '{A,T,G,C}' for each of the k blanks
echo 'kmers using bash to expand:' $s > /dev/stderr
bash -c "echo  $s";     # let brace expanion of inferior bash compute the cross product
}

mutate(){
usage="$FUNCNAME <seq> <mutations>";
	echo $1 | perl -ne 'chomp; my @seq=split //,$_;
		my @b=("A","C","G","T");
		print $_;
		for(my $i=0; $i<=$#seq; $i++){
			foreach my $bi (@b){ if($bi ne  $seq[$i]){
				my @s=@seq;
				$s[$i]=$bi;
				print " ",join( "",@s);
			}}		
		}
		print "\n";
	'	
}



test__mutate(){
echo \
"AAAA CAAA GAAA TAAA ACAA AGAA ATAA AACA AAGA AATA AAAC AAAG AAAT" > exp
mutate "AAAA" > obs
check obs exp
rm obs exp
}
#test__mutate


