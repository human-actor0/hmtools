#!/bin/bash 
#. $HMHOME/src/bed.sh
## Author Hyunmin Kim
## contact human.gim@gmail.com
mycat(){
        if [[ ${1##*.} = "gz" ]];then
                gunzip -dc $1;
        elif [[ ${1##*.} = "bz2" ]];then
                bunzip2 -dc $1;
        elif [[ ${1##*.} = "zip" ]];then
                unzip -p $1;
        else
                cat $1;
        fi
}
split_by_chrom(){
	awk -v OFS="\t" -v O=$2 '{
		fout=O"/"$1;
		print $0 >> fout;
	}' $1
	echo `ls $2/*`;	## return a list of splited files
}

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
THIS_RUN="$BASH_SOURCE $@";
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
	    }
	close(F);
	open (F,"BED") or die "$!";
	while(<F>){ chomp; my @a=split /\t/,$_;
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
	close(F);
	'			
	cmd=${cmd//LEFT/$LEFT};
	cmd=${cmd//RIGHT/$RIGHT};
	cmd=${cmd//STRAND/$STRAND};

	mkdir -p tmpd
	tmpd="tmpd";
	echo "$THIS_RUN " >&2
	for f in `split_by_chrom $BED $tmpd`;do
		chrom=${f##*/};
		echo " .. process $chrom .." >&2
		fa=`ls  $FA* | grep $chrom.fa`;
		if [[ $fa != "" &&  -f $fa ]];then
			fa1=$tmpd/$chrom.fa
			mycat $fa > $fa1
			cmd1=$cmd;
			cmd1=${cmd1//BED/$f};
			cmd1=${cmd1//FA/$fa1};
			perl -e "$cmd1"
		else
			echo "$fa not exist" >&2
		fi
	done
	rm -rf $tmpd
