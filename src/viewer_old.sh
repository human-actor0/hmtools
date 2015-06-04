#!/bin/bash
. $HMHOME/src/root.sh
. $HMHOME/src/polya.sh
BSIZE=100;
TYPE="polya";
echo "#$THIS $@" >&2
usage="
USAGE: $THIS [options] <genomic_range> <gene.bam> <bam> 
 [options]:
	-m <experiment type> : polya (default)
	-x <float> : zoom in/out
	-b <int>: bin size (default $BSIZE)
"
ZOOM=1;
while getopts "b:x:m:" arg
do
	case $arg in 
		b) BSIZE=${OPTARG};;
		x) ZOOM=${OPTARG};;
		m) TYPE=${OPTARG};;
		?) echo "$usage" ; exit 1;;
	esac
done
shift $(( OPTIND - 1 ))
if [ $# -ne 3 ];then
	echo "$usage"; exit 1;
fi
GRANGE=$1;
genes=$2;
chips=$3;
tmpd=`make_tempdir`
echo $GRANGE;
GRANGE=`perl -e '$ARGV[0]=~/(\w+):(\d+)-(\d+)/;
	my $x=$ARGV[1]; my $l=int(($3-$2)/2); my $c=int(($2+$3)/2);
	my $s= $c - $l*$x; my $e= $c + $l*$x;
	print $1,":",$s,"-",$e;' $GRANGE $ZOOM`
samtools view -b $genes $GRANGE | bamToBed -bed12 > $tmpd/genes
samtools view -b $chips $GRANGE | bamToBed | modify_score - count | point - | sum_score - > $tmpd/chips

#cat $tmpd/chips | sort -k2,2n | awk -v OFS="\t" '{ print $1,$2,$5;}' \
# | groupBy -g 1,2 -c 3 -o sum \
# | awk -v OFS="\t" '{ print $1,$2,$3;}' > tmp1
# cat tmp1 | groupBy -g 1 -c 2,3 -o collapse,collapse > tmp
# cat tmp | cluster_1d_kmeans0.sh -k 10 -c 2,3 - \
# | perl -ne 'chomp; my @a=split /\t/,$_;
#	my @s=split /,/,$a[$#a - 3];
#	my @e=split /,/,$a[$#a - 2];
#	my @c=split /,/,$a[$#a - 1];
#	my @h=split /,/,$a[$#a];
#	for(my $i=0;$i<=$#s;$i++){
#		print $a[0],"\t",$s[$i],"\t",$s[$i] + 100,"\tcluster\t",$h[$i],"\t+\n";
#	}
#' > $tmpd/clusters
cat $tmpd/clusters
cmd='
	a=strsplit("GRANGE",":")[[1]];
	chrom=a[1];
	grange=as.integer(strsplit(a[2],"-")[[1]]);
	genes=read.table("GENES",header=F);
	#clusters=read.table("CLUSTERS",header=F);
	#clusters;
	tt=read.table("CHIPS",header=F);
	my=min(tt$V5); My=max(tt$V5);
	ylim=c(my - 0.3* (My-my), My);
	xlim=grange;
	png("out.png");
	plot(main="GRANGE",tt[tt$V6=="+",2],tt[tt$V6=="+",5],type="l",ylim=ylim,xlim=xlim);
	#w=clusters[,3]-clusters[,2];
	#abline(v=clusters[,2]);
	lines(tt[tt$V6=="-",2],tt[tt$V6=="-",5],lty=2)

	## make gene track
	gene_track.top = my;
	gene_track.bot = ylim[1];
	gene_track.n = nrow(genes);
	gene_track.h = ((gene_track.top-gene_track.bot)/gene_track.n);
	for( i in 1:gene_track.n){
		s=genes[i,2];
		e=genes[i,3];
		st=genes[i,6];
		ybot=gene_track.top - i*gene_track.h;
		ytop=ybot+gene_track.h/2;
		if ( st == "+"){
			rect(s, ybot,e,ytop);
		}else{
			rect(s, ybot,e,ytop,lty=2);
		}
	}
	dev.off();
'
cmd=${cmd/GENES/$tmpd\/genes};
cmd=${cmd/CHIPS/$tmpd\/chips};
cmd=${cmd/BSIZE/$BSIZE};
cmd=${cmd/CLUSTERS/$tmpd\/clusters};
cmd=${cmd//GRANGE/$GRANGE};
echo "$cmd" |  R --no-save -q 

rm -rf $tmpd
