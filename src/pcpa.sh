. $HMHOME/src/polya.sh;
. $HMHOME/src/bed.sh;
. $HMHOME/src/stat.sh;

pcpa.prep(){
usage(){ echo "
 $FUNCNAME [options] <gene.bed> <3utr.bed> <name>:<file> [<name>:<file>..]
	<file>: this can be filtered polyA sites or clusters
		 (any of intermediate files are acceptable,  pa.point .. | pa.filter .. | pa.cluster .. )	
	[options]: 
	 -d <int> : maximum distance to be merged (default = 10)
"
}
local OPTIND d=10 g= u= tmpd=`mymktempd` files= names=
while getopts ":d:g:u:n:" opt; do
case $opt in
	d) d=$OPTARG;;
	\?) echo "Invalid option: -$OPTARG" >&2; usage ;;
esac
done
shift $(( OPTIND -1 ))
if [ $# -lt 3 ];then usage; return; fi

for f in ${@:3}; do
	files="$files ${f#*\:}"
	names="$names ${f%\:*}"
done
local tmp=( $names );
echo "chrom start end geneinfo score strand ${tmp[@]/#/count.} ${tmp[@]/#/total.}" | tr " " "\t"
pa.cluster -d $d ${files[*]} | bed.group - $1 | perl -ne '
	chomp; my @a=split /\t/,$_;
	my $g=$a[$#a]; $g=~s/@/;/g;
	next if($g eq "NA");
	$a[3]=$g;
	my $id=join("@",@a[0..5]);
	print "$g\t$id\t",join("\t",@a[6..($#a-1)]),"\n";
' | stat.gix2gixx -  | cut -f2- | tr "@" "\t" \
  | intersectBed -a stdin -b $2 -s -v 

rm -rf $tmpd
}

pcpa.prep.test(){
echo \
"c	100	200	gene1	0	+
c	300	400	gene2	0	+" > tmp.gene

echo \
"c	150	200	3utr1	0	+" > tmp.3utr

echo \
"c	10	11	pa0	1	+
c	110	111	pa1	1	+
c	120	121	pa2	1	+
c	150	151	pa3	1	+
c	310	311	pa4	1	+
c	399	400	pa5	1	+" > tmp.ctrpa

echo \
"c	110	111	pa1	10	+
c	120	121	pa2	10	+
c	310	311	pa3	10	+
c	300	301	pa5	10	+
c	399	400	pa4	10	+" > tmp.trtpa
pcpa.prep -d 1 tmp.gene tmp.3utr ctr:tmp.ctrpa trt:tmp.trtpa
rm tmp.*
}

pcpa.test(){
usage(){
	echo "
	$FUNCNAME <pcpa.prep> <method> <ctr> <trt>
	<method>: proportion
	<ctr> : <columnname of ctr> (default ctr)
	<trt> : <columnname of trt> (default trt)
	"
}
	if [ $# -ne 4 ];then usage;return; fi
	if [ $2 = "proportion" ]; then
	cat $1 | run_R '
		ctr="'$3'";trt="'$4'";
		tt=read.table("stdin",header=T);tt.cn=colnames(tt);
		j=c(
			which(paste("count.",ctr,sep="") == tt.cn),
			which(paste("count.",trt,sep="") == tt.cn),
			which(paste("total.",ctr,sep="") == tt.cn),
			which(paste("total.",trt,sep="") == tt.cn)
		)
		## global total counts
		tt.u=unique(tt[,c(4,j[3],j[4])])
		ctr.n=sum(tt.u[,2])
		trt.n=sum(tt.u[,3])

		pval=apply(tt[,j],1,function(x){
			if(sum( x[1:2] == x[3:4] ) == 2 || sum(x[3:4] - x[1:2]) == 2 || sum(x[3:4] == 0) > 0
			){ p=1;
			}else{ p=prop.test(x[1:2],n=x[3:4])$p.value;}
			return(p);
		})

		relFC=apply(tt[,j],1,function(x){
			x1=x+0.5;##pseudo counting
			return(log2( x1[2]*x1[3]/(x1[1]*x1[4])));
		})
		totFC=log2( (tt[,j[4]]+0.5)*ctr.n/(tt[,j[3]]+0.5)/trt.n);

		tt=tt[,c(1:6,j)];
		tt$log2RFC=relFC;
		tt$log2TFC=totFC;
		tt$pval=pval;
		write.table(tt,file="stdout",col.names=T,row.names=F,quote=F,sep="\t");
	' 
	else
		echo "method, $2 is not supported"
		usage; return;
	fi
}
pcpa.test.test(){
echo  \
"chrom	start	end	geneinfo	score	strand	count.ctr	count.trt	total.ctr	total.trt
c	310	311	c;300;400;gene2;0;+	11.00	+	1	10	2	30
c	300	301	c;300;400;gene2;0;+	10.00	+	0	10	2	30
c	399	400	c;300;400;gene2;0;+	11.00	+	1	10	2	30
c	110	111	c;100;200;gene1;0;+	11.00	+	1	10	3	20
c	120	121	c;100;200;gene1;0;+	11.00	+	1	10	3	20" \
| pcpa.test - proportion ctr trt 
}

pcpa.distance(){
usage(){ echo "
	find distance to the closes TSS
	usage: $FUNCNAME pcpa.bed gene.bed
"
}
	if [ $# -ne 2 ]; then usage; return; fi
	local tmpd=`mymktempd`; 
	mycat $1 > $tmpd/a;
	local nf=`head -n 1 $tmpd/a | awk '{print NF;}'`;
	bed.enc $2 | intersectBed -a $tmpd/a -b stdin -wa -wb -s \
	| perl -e 'use strict; my %res=(); my $nf='$nf';
		while(<STDIN>){ chomp; my @a=split/\t/,$_;
			my $d= $a[5] eq "+" ? $a[1]-$a[$nf+1] : $a[$nf+2]-$a[2];
			$res{ join("\t",@a[0..5]) }{ $d } = $a[$nf+3];
		}
		foreach my $k (keys %res){
			print $k,"\t";
			foreach my $d (sort {$a<=>$b} keys %{$res{$k}}){
				print $d,"\t",$res{$k}{$d},"\n";
				last ;
			}
		}
	'
	#| awk -v OFS="\t" '{ s=$2-$8; e=$9-$3; print $0,s,e;}'
	rm -r $tmpd
}

pcpa.distance.test(){
echo \
"chr1	100	200	t1	0	+
chr1	150	250	t2	0	+
chr1	90	300	t3	0	+" > tmp.a

echo \
"chr1	80	90	p1	0	+
chr1	110	120	p2	0	+
chr1	270	280	p3	0	+" > tmp.b

pcpa.distance tmp.b tmp.a


}
pcpa.summary(){
	local gene=$1; local utr=$2; local cls=$3;	
	cat $cls > tmp.0; cls=tmp.0;
	cat $cls | cut -f7- | awk '{print "#total\t",$0;}' | stat.sum -	
	intersectBed -a $cls -b $gene -wa -u -s > tmp.gene
	intersectBed -a tmp.gene -b $utr -v -s > tmp.pcpa
	cut -f7- tmp.gene | awk '{ print "#gene\t",$0;}' | stat.sum -
	cut -f7- tmp.pcpa | awk '{ print "#pcpa\t",$0;}' | stat.sum -
}

pcpa.cor(){
usage(){ echo "
	$FUNCNAME <file> [<method>]
	<method>: pearson, kendall, spearman (default)
"
}
	if [ $# -lt 1 ];then usage; return; fi
	local m=${2:-pearson}
	cat $1 | run_R '
	tt=read.table("stdin",header=F);
	J=as.integer((ncol(tt)-1)/2);
	ij=combn(1:J,2);
	x=tt[,2:(J+1)];
	y=tt[,(J+2):(2*J+1)];
	d=x/y;
	head(d)
	r.cor=apply(ij, 2, function(x){
		xx=d[,x[1]]; yy=d[,x[2]];
		return(cor(xx, yy, method="'$m'",use="pair"));
	})
	o=cbind(t(ij),r.cor)
	write.table(file="stdout", o ,sep="\t",quote=F,col.names=F,row.names=F);
	' 
}
pcpa.cor.test(){
echo \
"A	1	2	3	4	4	7	
B	3	4	5	6	5	5
C	5	6	7	8	9	10" | pcpa.cor -
		
}
pcpa.countbygene(){
usage(){ echo "
 $FUNCNAME [options] -g <gene.bed> -u <3utr.bed> <bed> [<bed> ..]
	[options]:
	 -n <name1>[,<name2> ..]  
	 -d <maxd> : maximum distance of peaks to be merged 
"
}
local OPTIND d=10 g u tmpd=`mymktempd` names=
while getopts ":d:g:u:n:" opt; do
case $opt in
	d) d=$OPTARG;;
	g) g=$OPTARG;;
	u) u=$OPTARG;;
	n) names=$OPTARG;;
	\?) echo "Invalid option: -$OPTARG" >&2; usage ;;
esac
done
shift $(( OPTIND -1 ))
if [ $# -lt 1 ];then usage; return; fi
	f(){
	 perl -e 'use strict; my %res=();
		while(<STDIN>){ chomp; my @a=split/\t/,$_;
			my $g=join("@",@a[0..5]);
			my $i=0;
			foreach my $e (@a[12..$#a]){
				$res{$g}[$i] += $e;
				$i++;
			}
		}
		foreach my $g (keys %res){
			print $g,"\t",join("\t",@{$res{$g}}),"\n";
		}'
	}
	echo "# files=[ $@ ]";
	pa.cluster -d $d $@ > $tmpd/a
	intersectBed -a $g -b $tmpd/a -wa -wb -s  | f - | sort -k1,1 > $tmpd/g
	intersectBed -a $tmpd/a -b $u -v -s  \
	| intersectBed -a $g -b stdin -wa -wb -s  | f - | sort -k1,1 > $tmpd/pcpa
	join -e NA $tmpd/pcpa $tmpd/g   | tr " " "\t"

	rm -rf $tmpd
#	bed.group tmp.pcpa $gene \
#	| perl -ne 'chomp; my@a=split/\t/,$_;print $a[$#a],"\t",join("\t",@a[6..($#a-1)]),"\n";' | stat.sum -  > tmp.a
#
#	bed.group tmp.gene $gene \
#	| perl -ne 'chomp; my@a=split/\t/,$_;print $a[$#a],"\t",join("\t",@a[6..($#a-1)]),"\n";' | stat.sum -  > tmp.b
#	
#	run_R 'tt=read.table("tmp.a",header=F);
#		tt=merge(tt, read.table("tmp.b",header=F),by=1,all.x=T);
#		tt[ is.na(tt) ] = 0;
#		write.table(tt,file="stdout",sep="\t",quote=F,col.names=F,row.names=F);
#	'

	#bed.group $cls $gene | head #perl -ne 'chomp;my@a=split/\t/,$_; print join("\t",@a[7..($#a-1)]),"\n";' | head #awk '$(NF) != "NA"{ for(i=0;i<NR;i++){ sprintf("\t' | cut -f 7,-1 | head #awk '{print "gene\t"$0;}' | stat.sum -	
	#intersectBed -a $cls -b $utr -s -v \
	#bed.group - $gene | cut -f7,-1 | awk '{print "pcpa\t"$0;}' | stat.sum -	
}




f(){
local OPTIND g u
while getopts ":g:u:i:" opt; do
case $opt in
	g) g="$OPTARG";;
	u) u="$OPTARG";;
	\?) echo "Invalid option: -$OPTARG" >&2; usage ;;
esac
done
shift $(( OPTIND -1 ))
if [ $# -lt 1];then usage; fi 

	local tmpd=`mymktempd`;
	local fa=( ${3//,/ } ); local fb=( ${4//,/ } )

	perl -e 'my $na='${#fa[@]}'; my $nb='${#fb[@]}';
		print "chrom\tstart\tend\tname\tscore\tstrand";
		for(my $i=1;$i<= $na; $i++){ print "\tctr".$i.".c1";}
		for(my $i=1;$i<= $nb; $i++){ print "\ttrt".$i.".c1";}
		print "\tgene\tpcpa\n";
	'
	pa.cluster $d ${fa[*]} ${fb[*]} | bed.group - $1 | bed.group - $2 \
	| awk -v OFS="\t" '{ if($(NF-1) != "NA" && $(NF)=="NA"){ $(NF)=1;}else{ $(NF)=0;}} 1' | sort -u
	rm -rf $tmpd	
}

pcpa.filter(){
usage(){ echo "
 FUNCTION: find PCPAs by gene not in 3\'UTRs
 OUTPUT: chrom start end gene_name pcpa_count strand total_gene_count
 USAGE: $FUNCNAME [options] <polya.bed> <gene.bed> <3utr.bed>
    [options]:
	-s : count pcpa in the same strand with gene.bed
	-S : count pcpa in the same strand with gene.bed
" 
}
local OPTIND s="-s";
while getopts ":sS" opt; do
case $opt in
	s) s="-s";;
	S) s="-S";;
	\?) echo "Invalid option: -$OPTARG" >&2; usage;return ;;
esac
shift $(( OPTIND -1 ))
done
if [ $# -ne 3 ];then usage; return; fi
	intersectBed -a $1 -b $2 $s -wa -wb \
	| perl -e 'use strict; my %res=(); my %tot=();
	while(<STDIN>){ chomp; my @a=split/\t/,$_;
		my $v=join("\t",@a[0..5]);	
		my $k=join("@",@a[6..$#a]);	
		$tot{$k} += $a[4];
		$res{$k}{$v}=1;
	}
	foreach my $k (keys %res){
		foreach my $v (keys %{$res{$k}}){
			print $v,"\t",$k,"\t",$tot{$k},"\n";	
		}	
	}
	' | intersectBed -a stdin -b $3 $s -v
}
pcpa.filter.test(){
echo \
"c	100	200	gene1	0	+
c	300	400	gene2	0	+" > tmp.gene

echo \
"c	150	200	3utr1	0	+" > tmp.3utr

echo \
"c	10	11	pa0	1	+
c	110	111	pa1	1	+
c	120	121	pa2	1	+
c	150	151	pa3	1	+
c	310	311	pa4	1	+
c	399	400	pa5	1	+" > tmp.ctrpa

echo \
"c	110	111	pa1	10	+
c	120	121	pa2	10	+
c	310	311	pa3	10	+
c	399	400	pa4	10	+" > tmp.trtpa

pcpa.filter tmp.ctrpa tmp.gene tmp.3utr 
rm tmp.*
}

pcpa.relfreq(){
usage(){ 
echo " 
$FUNCNAME [options] <pcpa.prep> 
	[options]:
	-a <num> : minimum value of the sum of all expriments (default 10)
	-b <num> : minimum absolute value of the relative log2 fold change  (default 1)
"
return;
}
local OPTIND a b
a=10;b=1;
while getopts ":" opt; do
case $opt in
	\?) echo "Invalid option: -$OPTARG" >&2; usage ;;
esac
done
shift $(( OPTIND -1 ))
if [ $# -ne 1 ];then usage; fi
cat $1 | run_R 'thre.a='$a'; thre.b='$b';
	f=function(d){
		if(is.null(dim(d))){ return(d);
		}else{ return(apply(d,1,sum)); }
	}
		
	tt=read.table("stdin",header=T);
	cn=colnames(tt);
	x=f(tt[,grep("ctr",cn)]);
	y=f(tt[,grep("trt",cn)]);
	tt$ctr.g=ave(x, tt$gene,FUN=sum);
	tt$trt.g=ave(y, tt$gene,FUN=sum);

	a=x+y;
	b=log2((y+1)/(tt$trt.g+1)*(tt$ctr.g+1)/(x+1));
	tt$log2fc=b;
	write.table(tt, file="stdout",sep="\t",col.names=T,row.names=F,quote=F);
' 
}
#pcpa.filter1(){
#usage(){ 
#echo " 
#$FUNCNAME [options] <pcpa.prep> 
#	[options]:
#	-a <num> : minimum value of the sum of all expriments (default 10)
#	-b <num> : minimum absolute value of the relative log2 fold change  (default 1)
#"
#return;
#}
#local OPTIND a b
#a=10;b=1;
#while getopts ":a:b:" opt; do
#case $opt in
#	a) a=$OPTARG;;
#	b) b=$OPTARG;;
#	\?) echo "Invalid option: -$OPTARG" >&2; usage ;;
#esac
#done
#shift $(( OPTIND -1 ))
#if [ $# -ne 1 ];then usage; fi
#cat $1 | run_R 'thre.a='$a'; thre.b='$b';
#	f=function(d){
#		if(is.null(dim(d))){ return(d);
#		}else{ return(apply(d,1,sum)); }
#	}
#		
#	tt=read.table("stdin",header=T);
#	cn=colnames(tt);
#	x=f(tt[,grep("ctr",cn)]);
#	y=f(tt[,grep("trt",cn)]);
#	tt$ctr.g=ave(x, tt$gene,FUN=sum);
#	tt$trt.g=ave(y, tt$gene,FUN=sum);
#
#	a=x+y;
#	b=log2((y+1)/(tt$trt.g+1)*(tt$ctr.g+1)/(x+1));
#	tt$log2fc=b;
#	ix.p = tt$pcpa == 1 & a >= thre.a & b >= thre.b;
#	ix.n = tt$pcpa == 1 & a >= thre.a & b <= -thre.b;
#
#	tt$sig=0;
#	tt[ix.p,]$sig = 1;
#	tt[ix.n,]$sig = -1;
#	write.table(tt, file="stdout",sep="\t",col.names=T,row.names=F,quote=F);
#' 
#}
#
#pcpa.prep.test
