#!/bin/bash

util.grep(){
        local OPTIND C D; C=-1; D="";
        while getopts ":c:d:" arg; do
                case $arg in
                        c) C=$OPTARG;;
                        d) D=$OPTARG;;
                        \?) echo "Invalid -${OPTARG}"; return;;
                esac
        done
        shift $(( $OPTIND - 1 ))

if [ $# -lt 2 ];then echo "
$FUNCNAME [options] <file> <pattern_file>
 -c <int>: column 
 -d <delim> : secondary delimiter
"; fi

	cat $1 | perl -e 'use strict; my $f2="'$2'"; my $C='$C'; my $D="'$D'";
		open my $fh, "<",$f2 or die "cannot open $f2";
		my %P=(); while(<$fh>){chomp; $P{$_}++;}
		close($fh);

		while(<STDIN>){chomp; my @a=split/\t/,$_;
			my $hit=0;
			if($C < 0){
				foreach my $ai (@a){
					if(defined $P{$ai}){ $hit=1; last; }	
				}
			}elsif(defined $P{$a[$C-1]}){
				$hit=1;
			}elsif($D ne ""){
				foreach my $s (split /$D/,$a[$C-1]){
					if(defined $P{$s}){ $hit=1;}
				}
			}
			
			if($hit){
				print $_,"\n";
			}
		}
	'
}

util.grep.test(){
echo \
"a
b
c" > tmp.a

echo \
"1	a	3
a	2	2
a	b	2;b
3	b	1" > tmp.b 
util.grep -c 2 tmp.b tmp.a 
util.grep -c 1 tmp.b tmp.a 
echo "d=;"
util.grep -c 3 -d ";" tmp.b tmp.a 
rm tmp.a tmp.b
}
