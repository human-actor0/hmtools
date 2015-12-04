#!/bin/bash
. $HMHOME/src/root.sh
. $HMHOME/src/bed.sh

seq.read(){
usage="
USAGE: $FUNCNAME <fa> <bed> [options]
 [options]: 
	-s : reverse complements for the negative strand
	-i : ignore case 

"
local option=${@:3};
if [ $# -lt 2 ]; then echo "$usage"; return; fi
	local tmpd=`mymktempd`;
	mycat $2 > $tmpd/a;	
	mycat $1 | perl -e 'use strict; 
		sub revc{
			my ($seq) = @_;
			$seq =~ tr/ACGTacgt/TGCAtgca/;
			return join("",reverse(split //,$seq));
		}

		my $file="'$tmpd/a'"; my $opt="'"$option"'";
		my %ref=(); my $chrom="";
		while(<STDIN>){ chomp;
			if($_=~/>(\S+)/){ $chrom=$1; next; }		
			$ref{$chrom} .= $_;
		}
		#foreach my $c (keys %ref){ print $c," ",substr($ref{$c},0,10),"\n"; }
		open(my $fh, "<", $file) or die "$!";
		while(<$fh>){chomp;
			my @a=split/\t/,$_;
			my $seq="NULL";
			if(defined $ref{$a[0]}){
				$seq=substr($ref{$a[0]},$a[1],$a[2]-$a[1]);
				if($opt =~ /s/ && $#a >=5 && $a[5] eq "-"){
					$seq=revc($seq);
				}
			}
			if($opt=~/lc/){ $seq = lc $seq; }elsif ($opt=~/uc/){ $seq = uc $seq; }
			print join("\t",@a),"\t",$seq,"\n";
		}
		close($fh);
		
	'
	rm -rf $tmpd;
}
seq.read.test(){
echo	"chr22	17265299	17265350	ENSG00000172967	0	-
chr22	17280914	17280954	ENSG00000172967	0	-
chr2	217280752	217280979	ENSG00000138375	0	+" \
| seq.read $HMHOME/data/hg19_chr22.fa.gz - -s -lc
}

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


