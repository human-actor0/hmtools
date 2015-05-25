#!/bin/bash 
. util.sh
THIS=${BASH_SOURCE##*/}
echo "#$THIS $@" >&2;
FIVE=1000;
THREE=1000;
B=200; ## Binsize
usage="
USAGE: $THIS [options] <target.bed> <feature.bed>
 [options]:
	-5 <int>: 5' flank size (default $FIVE)
	-3 <int>: 3' flank size (default $THREE)
	-b <int>: binsize   (default $B)
"
while getopts "q:b:5:3:" arg; do 
	case $arg in 
		5) FIVE=${OPTARG};;
		3) THREE=${OPTARG};;
		b) B=${OPTARG};;
		?) echo "$usage"; exit 1;;
	esac
done
shift $(( OPTIND - 1 ))
if [ $# -ne 2 ];then
	echo "$usage"; exit 1;
fi
#tmpd='tmp'; mkdir -p $tmpd;
tmpd=`makeTempDir`
cat $1 > $tmpd/a
cat $2 > $tmpd/b
getChroms $tmpd/a | sort > $tmpd/a.chrom
getChroms $tmpd/b | sort > $tmpd/b.chrom 

for chrom in `join $tmpd/a.chrom $tmpd/b.chrom `;do
	echo " running $THIS .. $chrom" >&2
	awk -v C=$chrom '$1==C' $tmpd/a > $tmpd/a.$chrom
	awk -v C=$chrom '$1==C' $tmpd/b > $tmpd/b.$chrom
	#samtools view -bq $Q $2 $chrom | bamToBed | frag_shift - $F \
	windowBed -a $tmpd/a.$chrom -b $tmpd/b.$chrom -sw -l $FIVE -r $THREE \
	| awk -v OFS="\t" -v B=$B '{ 
		## relative position 
		l=$3-$2; s=$8-$2; e = $9-$2; f=$10; r=4; 
		bs=int(s/B); be=int((e-1)/B);
		for(i=bs; i<=be; i++){
			b=i*B;	
			if(b<0){ r=5; b=b-B;
			}else if(b >= l || l<B){ r=3;}
			n=$1":"$2"-"$3"("$6")"$4;
			print n"@"r"@"f"@"b;
		}
	}'| sort | uniq -c | awk -v OFS="\t" '{ split($2,a,"@"); print a[1],a[2],a[3],a[4],$1;}'
	exit 1
done

