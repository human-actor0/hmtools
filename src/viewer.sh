#!/bin/bash
. $HMHOME/src/root.sh
. $HMHOME/src/polya.sh

usage="
USAGE: $THIS [options] <tracks> <outdir>
 [options]:
"
while getopts "h" arg
do
	case $arg in 
		?) echo "$usage" ; exit 1;;
	esac
done
shift $(( OPTIND - 1 ))
if [ $# -ne 2 ];then
	echo "$usage"; exit 1;
fi

using_sushi(){
echo '
	library("Sushi");
'
}

## make bedGraph style 
draw_point(){
cmd='
	a=strsplit("GRANGE",":")[[1]];
	chrom=a[1]; grange=as.integer(strsplit(a[2],"-")[[1]]);

	tt=read.table("INPUT",header=F); d=tt;
	d[tt[,6]=="-",5] = - d[tt[,6]=="-",5]; ## negative strand
	x=c(d[,2],d[,3],d[,2]-0.1,d[,3]+0.1);
	y=c(d[,5],d[,5],rep(0,nrow(d)),rep(0,nrow(d)));
	ix=sort(x,index=T)$ix;
	x=x[ix];y=y[ix];

	png("OUTPUT",width = 480, height = 200);
	plot(x,y,type="l", main="GRANGE",xlim=grange)
	dev.off();
'
	cmd=${cmd//INPUT/$1};
	cmd=${cmd//GRANGE/$2};
	cmd=${cmd//OUTPUT/$3};
	echo "$cmd";
	echo "$cmd" | R --no-save -q
}

## make gene track style
draw_genetrack(){
cmd='
	genes=read.table("INPUT",header=F);
	a=strsplit("GRANGE",":")[[1]];
	chrom=a[1]; grange=as.integer(strsplit(a[2],"-")[[1]]);

	gene_track.n = nrow(genes);
	gene_track.top = gene_track.n;
	gene_track.bot = 0;
	gene_track.h = ((gene_track.top-gene_track.bot)/gene_track.n);
	
	ylim=c(gene_track.bot, gene_track.top);
	xlim=grange;
	png("OUTPUT",width =480, height=200);
	plot(NULL,main="GRANGE",ylim=ylim,xlim=xlim);
	for( i in 1:gene_track.n){
		s=genes[i,2];
		e=genes[i,3];
		st=genes[i,6];
		ybot=gene_track.top - i*gene_track.h;
		ytop=ybot+gene_track.h/2;
		lengths=as.integer(strsplit(as.character(genes[i,11]),",")[[1]])
		starts=as.integer(strsplit(as.character(genes[i,12]),",")[[1]])
		for (j in 1:genes[i,10] ){
			start=s+starts[j];	
			end=start+lengths[j];
			rect(start,ybot,end,ytop);
			if (j < genes[i,10]){
				if ( st == "+"){
					lines(c(end, starts[j+1] +s),c(ybot,ybot));
				}else{
					lines(c(end, starts[j+1] +s),c(ytop,ytop));
				}
			}
		}
	}
	dev.off();
'
	cmd=${cmd//INPUT/$1};
	cmd=${cmd//GRANGE/$2};
	cmd=${cmd//OUTPUT/$3};
	echo "$cmd";
	echo "$cmd" | R --no-save -q
}
O=${2/%\/};
for f in ${1/\//}/*\:*-*;do
	interval=${f##*/};
	for t in $f/*; do
		track_type=${t##*/};
		for n in $t/*; do
			if [ ! -f $n/a.bed ];then
				"ERROR: $n/a.bed not exists" >&2
				continue;
			fi
			track_name=${n##*/};
			if [ "$track_type" = "gene" ];then
				draw_genetrack $n/a.bed $interval $O/$interval.$track_type.$track_name.png
			elif [ "$track_type" = "polya" ]; then
				draw_point $n/a.bed $interval $O/$interval.$track_type.$track_name.png
			fi
		done
	done
done
