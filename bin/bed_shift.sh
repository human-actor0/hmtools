#!/bin/bash

usage="
USAGE: $0 [options] <bed> 
        [options]:
         -l <int>: left flank
         -r <int>: right flank
         -s      : strand specific
"
LEFT=0; RIGHT=0; STRAND=0;
while getopts "hsl:r:" arg; do
        case $arg in
                l) LEFT=${OPTARG};;
                r) RIGHT=${OPTARG};;
                s) STRAND=1;;
                ?) echo "$usage"; exit 1;;
        esac
done
shift $(( OPTIND - 1 ))
if [ $# -ne 1 ]; then echo "$usage"; exit 1; fi

bed_shift(){
	cat $1 | awk -v OFS="\t" -v L=$LEFT -v R=$RIGHT -v S=$STRAND '{
		if(S == 1 && $6=="-"){ $2 = $2 - R; $3 = $3 + L;
		}else{ $2 = $2 - L; $3 = $3 + R; }
		if( $3 >= $2){
			print $0;
		}else{
			print "ERROR: ",$0 | "cat >&2"; 
		}
	}'
}
if [ $1 = "test" ];then
LEFT=10; RIGHT=20; STRAND=1;
echo \
"chr1	100	200	.	0	+
chr1	100	200	.	0	-
chr1	200	100	.	0	-" | bed_shift -
else
	bed_shift -
fi

