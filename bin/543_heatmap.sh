#!/bin/bash 
. $HMHOME/src/root.sh
THIS=${BASH_SOURCE##*/}
echo "#$THIS $@" >&2;

TOP=1000;
K=4;
O="out";
C="";
N="";
usage="
USAGE: $THIS [options] <matrix>
	<matrix>: id, sample,  feature1, ...
 [options]:
	-t <int> : use top-N entries (default $TOP)
	-k <int> : number of clusters (default $K);
	-c <colors> : colors in file orders (default $C)
	-n <names> : comma separated sample names
	-o <dir> : output dir (default $O)
"
while getopts "ht:k:o:c:n:" arg; do 
	case $arg in 
		t) TOP=${OPTARG};;
		k) K=${OPTARG};;
		o) O=${OPTARG};;
		c) C=${OPTARG};;
		n) N=`echo ${OPTARG} | perl -ne 'chomp;my @a=split/,/,$_; print join (",", map{"\"$_\""} @a)'`;;
		?) echo "$usage"; exit 1;;
	esac
done
shift $(( OPTIND - 1 ))
echo $#
if [ $# -lt 1 ];then
	echo "$usage"; exit 1;
fi
#FILES=`echo "$@" | perl -ne 'chomp;@a=split/\s+/,$_; print join(",", map {"\"$_\""} @a);'`
mkdir -p $O
rcmd='
set.seed(416);
library(gplots);
odir="OUT"; 
tt=read.table("stdin",header=T);
#tt=read.table("tt",header=T);
experiments=unique(tt$sample);
colors=1:length(experiments); COLORS
names=1:length(experiments); NAMES
m=tt[,-c(1:3,ncol(tt))];

## normalize each sample by total counts
n=ave(apply(m,1,sum), tt$sample,FUN=sum); 
write.table(file="stdout",minn)
m=m/n*10^6;## RPM

## use sample 0 as feature expression 
eden=apply(m[tt$sample==experiments[1],],1,sum);
ix=sort(eden,index=T,decreasing=T)$ix[1:TOP]
mk=t(apply(m[tt$sample==experiments[1],][ix,], 1, function(x){ 
	res=rep(3,length(x)); q=quantile(x,prob=c(.5,.7)); res[ x >= q[[1]] ]=1; res[x>= q[[2]] ]=2; res[x <=q[[1]] ] = 0; res } ));
clusters=kmeans(mk,K)$cluster
tmp=tt[tt$sample==experiments[1],][ix,] 
write.table(data.frame(id=tmp$id, group=clusters),file=paste(odir,"/clusters.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

for( k in unique(clusters)){
	m1=c();
	for( s in experiments){
		tmp = m[tt$sample==s,][ix,][clusters==k,];
		m1=rbind(m1,apply(tmp,2,mean));
		pdf(file=paste(odir,"/heatmap_s",s,"g",k,".pdf",sep=""))
		xx=as.matrix(tmp)
		col=colorpanel(256, "black", "grey", "white")
		xx=t(apply(xx,1,function(x){ exp(x-median(x))}));
		heatmap(xx, main=paste("sample=",s," n=",nrow(xx)," cluster=",k), Rowv=NA,Colv=NA,col=col);
		dev.off();
	}

    ylim=c(min(m1),max(m1));
    x=colnames(m1);
    png(paste(odir,"/plot_g",k,".png",sep=""));
    plot(as.numeric(m1[1,]),xaxt="n",col=colors[1], ylim=ylim,type="l",lwd=3,main=paste("class=",k," n=",sum(clusters==k)));
	if(nrow(m1)> 1){
		for(j in 2:nrow(m1)){
			lines(as.numeric(m1[j,]),col=colors[j],lwd=3);
		}
	}
    abline(v=which(x=="X3.0" | x=="X4.0"),lty=2,col=3)
    axis(1,1:ncol(m1),colnames(m1))
    dev.off();

}
png(paste(odir,"/legend.png",sep=""));
	plot(1,col="white");
	legend("center",legend=names,lty=rep(1,length(colors)),col=colors);
dev.off();
for( s in experiments){
    #pdf(file=paste(odir,"/heatmap_s",s,".pdf",sep=""))
    png(file=paste(odir,"/heatmap_s",s,".png",sep=""))
	oc=order(clusters);
    xx=as.matrix(m[tt$sample==s,][ix,][order(clusters),])
#	for( c in unique(clusters)){
#		xx[min(which(oc==c)),]=max(xx);
#	}
	#xx=t(apply(xx,1,function(x){ exp(x/sum(x)); }));
	col=colorpanel(256, "black", "grey", "white")
	#xx=t(apply(xx,1,function(x){ exp(x-median(x))}));
    heatmap(xx, main=paste("sample=",s," n=",nrow(xx)," clusters=",K), Rowv=NA,Colv=NA,scale="row",labRow=NA,col=col);
    #image(xx, main=paste("sample=",s," n=",nrow(xx)," clusters=",K),col=col);
    #heatmap.2(xx, main=paste("sample=",s," n=",nrow(xx)," clusters=",K), Rowv=NA,Colv=NA,trace="none",labRow=NA,scale="row",col=redgreen(256))
    dev.off();
}
'
rcmd=${rcmd//TOP/$TOP}
rcmd=${rcmd//K/$K}
rcmd=${rcmd//OUT/$O}
if [ ! $C = "" ];then 
	rcmd=${rcmd//COLORS/"colors=c($C)"}; 
else
	rcmd=${rcmd//COLORS/}; 
fi
if [ ! $N = "" ];then 
	rcmd=${rcmd//NAMES/"names=c($N)"}; 
else
	rcmd=${rcmd//NAMES/}; 
fi
cat $1 | run_R "$rcmd"
#hc.rows=hclust(dist(m[,-1]));
#plot(hc.rows);
#
#mtscaled=as.matrix(scale(tt[1:10,-1]));
#
#foreachChrom "$2"
#montage -geometry 100% -tile 2x2 ./*.png outfile.pdf
