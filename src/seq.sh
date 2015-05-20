#!/bin/bash 
. $HMHOME/src/root.sh
usage="
USAGE: $0 [options] <bed> <fasta|dir> 
	[options]:
	 -l <int>: left flank
	 -r <int>: right flank 
	 -s      : strand specific
"
LEFT=0;
RIGHT=0;
STRAND=0;
while getopts "hsl:r:" arg; do
	case $arg in
		l) LEFT=${OPTARG};;
		r) RIGHT=${OPTARG};;
		s) STRAND=1;;
		?) echo "$usage"; exit 1;;
	esac
done
shift $(( OPTIND - 1 ))
if [ $# -ne 2 ]; then
	echo "$usage"; exit 1;
fi
BED=$1; FA=$2;
cmd='
	use strict;
	sub revComp{
		my ($seq) = @_;
		$seq =~ tr/ACGTacgt/TGCAtgca/;
		return join("",reverse(split //,$seq));
	}
	my %seq=();
	my $chrom="";
   	open (F, "FA") or die "$!";
    while(<F>){ chomp;
		if($_=~/>([\w|\d]+)/){ $chrom=$1; next;}
        $seq{$chrom} .= $_;
    } close(F);
	while(<STDIN>){ chomp; my @a=split /\t/,$_;
		my $sseq = substr($seq{$a[0]},$a[1],$a[2]-$a[1]);	
		my $sseq_left="";my $sseq_right="";
		if(STRAND==0 || $a[5] eq "+"){
			$sseq_left = substr($seq{$a[0]},$a[1]-LEFT,LEFT) if( LEFT > 0);	
			$sseq_right = substr($seq{$a[0]},$a[2],RIGHT) if(RIGHT > 0);	
		}else{
			$sseq = revComp($sseq);
			$sseq_left = revComp(substr($seq{$a[0]},$a[2],LEFT)) if( LEFT > 0);	
			$sseq_right = revComp(substr($seq{$a[0]},$a[1]-RIGHT,RIGHT)) if(RIGHT > 0);	
		}
		print $_,"\t",uc $sseq_left,",",uc $sseq,",",uc $sseq_right,"\n";
	}
	'			
	cmd=${cmd//LEFT/$LEFT};
	cmd=${cmd//RIGHT/$RIGHT};
	cmd=${cmd//STRAND/$STRAND};

	tmpd=tmpd; mkdir -p $tmpd
	mycat $FA | awk -v O=$tmpd '{
		if(substr($0,1,1)==">"){ CHROM=substr($0,2,length($0));} 
		fout=O"/"CHROM; print $0 >> fout;
	}'

	
#	FA=${FA%\/}
#	if [ -f $FA ]; then
#		cmd=${cmd//BED/$BED};
#		cmd=${cmd//FA/$FA};
#		cat $BED | perl -e "$cmd"
#	elif [ -d $FA ]; then
#		tmpd=`make_tempdir`;
#		for f in `split_by_chrom $BED $tmpd`;do
#			chrom=${f##*/};
#			echo "$f $chrom"
#			fa=$FA/$chrom.fa
#			if [ -f $fa ]; then
#				cmd1=$cmd;
#				cmd1=${cmd1//BED/$tmp_bed};
#				cmd1=${cmd1//FA/$fa};
#				awk -v CHROM=$chrom '$1==CHROM' $f \
#				| perl -e "$cmd1" 
#			fi
#		done
#	fi
