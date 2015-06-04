#!/bin/bash
. $HMHOME/src/root.sh
. $HMHOME/src/polya.sh
Q=10; M=count; D=0; O=out
usage="
USAGE: $THIS [options] <genomic_range> <type>:<name>:<bam> [<type>:<name>:<bam>] 
 <type>: gene, polya, chip
 [options]:
  -q <int> : quality threshold (default $Q)
  -m <method> : count or phred score (default $M)
  -d <mind> : clustering parameter (default $D)
  -o <str> : output directory (default out)
"
while getopts "q:c:d:m:o:" arg
do
	case $arg in 
		q) Q=${OPTARG};;
		m) M=${OPTARG};;
		d) D=${OPTARG};;
		o) O=${OPTARG};;
		?) echo "$usage" ; exit 1;;
	esac
done
shift $(( OPTIND - 1 ))
if [ $# -lt 3 ];then
	echo "$usage"; exit 1;
fi
datq=( $@ );
for f in ${@:2};do
	a=( `echo "$f" | tr ":" "\t"` )
	odir=$O/$1/${a[0]}/${a[1]}; mkdir -p $odir 
	echo "${a[*]} =>  $odir .. ";
	if [ ${a[0]} = "gene" ];then
		#samtools view -b ${a[2]} $1 | bamToBed -bed12 > $odir/a.bed 
		echo "$1" | perl -ne '
			if( $_ =~ /(.+):(\d+)-(\d+)/ ){ 
				print $1,"\t",$2,"\t",$3,"\n";
			}
		' | intersectBed -a ${a[2]} -b stdin -wa -u  > $odir/a.bed
	elif [ ${a[0]} = "polya" ]; then
		eval "
		samtools view -bq $Q ${a[2]} $1 | bamToBed | modify_score - $M | point - > $odir/a.bed 
		"
	fi
done
