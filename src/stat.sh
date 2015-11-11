#!/bin/bash
. $HMHOME/src/root.sh
stat.sum(){
	cat $1 | perl -e 'use strict; my %res=();
	while(<STDIN>){ chomp; my @a=split/\t/,$_;
		$res{ $a[0] } += $a[1]; 
	}
	foreach my $k (keys %res){
		print $k,"\t",$res{$k},"\n";
	}'
}

stat.prep(){
usage="
FUNCT: make a merged table with headers 
USAGE: $FUNCNAME <trt1>[,<trt2 ..]  <ctr1>[,<ctr2..]
"; if [ $# -ne 2 ];then echo "$usage"; return; fi

        trt=`echo $1 | tr "," " " | quote -`;
        ctr=`echo $2 | tr "," " " | quote -`;
        run_R '
                trt=c('$trt'); ctr=c('$ctr');
                d=NULL;
                for( i in 1:length(trt)){
                        tt=read.table(trt[i],header=F);
                        colnames(tt)=c("id",paste("trt",i,".c1",sep=""), paste("trt",i,".c2",sep=""));
                        if( is.null(d)){ d=tt;
                        }else{ d=merge(d,tt,by="id",all=T); }
                }
                for( i in 1:length(ctr)){
                        tt=read.table(ctr[i],header=F);
                        colnames(tt)=c("id",paste("ctr",i,".c1",sep=""), paste("ctr",i,".c2",sep=""));
                        d=merge(d,tt,by="id",all=T);
                }
		d[ is.na(d) ] = 0;
		write.table(file="stdout",d, col.names=T,quote=F,sep="\t", row.names=F);
        ' 
}

stat.prep.test(){
echo \
"a	1	10
b	2	20
c	3	30" > trt1
echo \
"a	11	10
b	22	20
c	33	30" > trt2
echo \
"a	11	10
b	22	20
c	33	30" > ctr1

echo \
"id	trt1.c1	trt1.c2	trt2.c1	trt2.c2	ctr1.c1	ctr1.c2
a	1	10	11	10	11	10
b	2	20	22	20	22	20
c	3	30	33	30	33	30" > exp

	stat.prep trt1,trt2 ctr1  > obs
	check exp obs
	rm -rf trt1 trt2 ctr1 exp obs
}
add(){
	perl -e 'use strict; my %res=();
	foreach my $f (@ARGV){
		my $fh;
		open($fh, $f) or die "$!";
		while(<$fh>){chomp; my ($k,$v)=split/\t/,$_;
			$res{$k} += $v;
		}
		close($fh);	
	}
	foreach my $k (keys %res){ print $k,"\t",$res{$k},"\n";}
	' $@;
}
test__add(){
echo "a	1" > tmpa
echo "a	10" > tmpb
add tmpa tmpa tmpb tmpb tmpb
rm -f tmpa
}
max(){
	awk 'NR == 1 { max=$1;} { if($1 > max) max=$1;} END { print max;}'	
}
min(){
	awk 'NR == 1 { min=$1;} { if($1 < min) min=$1;} END { print min;}'	
}

cor(){
        cat $1 | R --no-save -q -e 'tt=read.table("stdin",header=F);cor(tt[,1],tt[,2],method="spearman" );' \
        | perl -ne 'chomp;if($_=~/\[1\] ([\-|\d|\.]+)/){print $1;}'
}

cor_to_hclust(){
	cmd='
		tt=read.table("stdin",header=F);
		x=as.character(unique(c(as.character(tt[,1]),as.character(tt[,2]))));
		tt;
		#tt=rbind(tt,data.frame(V1=x,V2=x,V3=1));
		tt[,3]=1-tt[,3];
		d=as.dist(xtabs(tt[, 3] ~ tt[,2] + tt[,1]));
		d;
		png("cor.png"); plot(hclust(d));dev.off();
	'
	cut -f 1 $1  > a
	cut -f 2 $1 >> a
	rm -f b;
	for f in `sort -u a`;do
		echo -e "$f\t$f\t1" >> b;
	done
	cat $1 b | sort -k1,2| run_R "$cmd"
	rm -rf a b
}

pcombine_fisher(){
usage="
usage: $FUNCNAME <file> 
"
if [ $# -ne 1 ]; then echo "$usage"; return; fi
cmd='
	fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)
	tt=read.table("stdin",header=F);
	cpval=unlist(apply(tt[,-1],1,fishersMethod));
	tt$cpval=cpval;
        write.table(file="stdout",tt,row.names=F,col.names=F,quote=F,sep="\t");
'
	cmd=${cmd/PCOL/$2};
	cat $1 | run_R "$cmd"
}
igx_to_igxnx(){
usage="
usage: $FUNCTNAME <file>
 <file>: id, group, counts delimited by [tab]s
"
cat $1 | perl -e 'use strict;
        my %d=();
        my %dt=();
        while(<STDIN>){ chomp;
                my ($i,$g,@x) = split /\t/,$_;
		if( !defined $d{$g} ){
			@{$d{$g}}=();
			@{$dt{$g}}=();
		}
                push @{$d{$g}}, [$i, $g, @x];
		for( my $i=0; $i < scalar @x; $i++){
			$dt{$g}[$i] += $x[$i];
		}
        }
        foreach my $g (keys %d){
	foreach my $ae (@{$d{$g}}){
		my @b=@{$dt{$g}};
		my @a=@$ae;
		for(my $i=0; $i <= $#b; $i++){
			$b[$i] -= $a[$i+2];
		}
		print join("\t",@a),"\t",join("\t",@b),"\n";
        }}
'
}
test__igx_to_igxnx(){
echo \
"1	g1	1	11
2	g1	2	22
3	g2	3	33" \
| igx_to_igxnx - > obs
echo \
"3	g2	3	33	0	0
1	g1	1	11	2	22
2	g1	2	22	1	11" > exp
check obs exp
rm obs exp
}
#test__igx_to_igxnx
stat.edger_test(){
	## INPUT: comma separated control and treatment EI file ( bed6 + x + nx )
	## OUTPUT: bed6 + logFC + pvalue 
	cat $1 | run_R '
	tt=read.table("stdin",header=T);
	cn=colnames(tt)[2:ncol(tt)];
	group=rep(1,length(cn)); group[ grep("trt",cn)]=2;
	event=rep(1,length(cn)); event[ grep(".c2",cn)]=2;
	
		
	#ix=apply(D[,7:ncol(D)], 1, min) > 0 & apply(D[,7:ncol(D)],1,max) > 10
	library(edgeR)
	y=DGEList(counts=tt[,2:ncol(tt)],group=factor(group));
	y=calcNormFactors(y);

	event.this=factor(event);
	group.this=factor(group);
	H1 <- model.matrix(~ event.this + group.this + event.this:group.this, data=y$samples )
	H0 <- model.matrix(~ event.this + group.this )
	coef <- (ncol(H0)+1):ncol(H1)
	#y=estimateCommonDisp(y)
	#y=estimateTagwiseDisp(y, trend="movingave")
	y = estimateGLMCommonDisp(y,H1);
	y = estimateGLMTrendedDisp(y,H1);
	y = estimateGLMTagwiseDisp(y,H1);

	fit=glmFit(y$counts, H1, y$tagwise.dispersion,offset=0,prior.count=0)
	llh=glmLRT(fit,coef=coef)

	res=data.frame(tt, logFC=llh$table$logFC, pval=llh$table$PValue)
	## chrom start end logFC pval
	write.table(res,file="stdout", col.names=T,row.names=F,sep="\t",quote=F);
	' $2
}
stat.edger_test.test(){
echo \
"id	trt1.c1	trt1.c2	trt2.c1	trt2.c2	ctr1.c1	ctr1.c2
a	1	10	11	10	11	10
b	2	20	22	20	22	20
c	3	30	33	30	33	30" | stat.edger_test -
}

padjust(){
usage="
usage: $FUNCNAME <file> <pvalue_index>
"
if [ $# -ne 2 ]; then echo "$usage"; return; fi
cmd='
	pcol=PCOL;
	tt=read.table("stdin",header=F);
	if( pcol < 0){
		pcol = ncol(tt) + pcol + 1;
	}
	tt$fdr=p.adjust(tt[,pcol],method="fdr");
        write.table(file="stdout",tt,row.names=F,col.names=F,quote=F,sep="\t");
'
	cmd=${cmd/PCOL/$2};
	cat $1 | run_R "$cmd"
}

lineartrend_test(){
cmd='
from scipy.stats.stats import pearsonr
from scipy.stats import chi2
import sys

def test_lineartrend(x1,y1,x2,y2):
	if len(x1) < 2 or len(x2) < 2 or sum(y1) == 0 or sum(y2) == 0:
		raise ValueError
	sx = [];
	sy = [];
	for i in range(len(x1)):
			for j in range(y1[i]):
				sx.append(x1[i]);
				sy.append(1); 
	for i in range(len(x2)):
			for j in range(y2[i]):
				sx.append(x2[i]);
				sy.append(2);
	if len(set(sx)) == 1 :
		return 0,1;
	r,pval = pearsonr(sx,sy);
	s = r*r*(len(sx)-1);
	pval = 1-chi2.cdf(s,1)
	return r,pval

def s2a(x):
	return [ int(float(x)) for x in x.split(",")];

for line in sys.stdin:
	if line[0] == "#":
		print line.rstrip();
		continue;
	a = line.rstrip().split("\t");
	x1,y1,x2,y2 = map(s2a, a[1:]);
	try:
		#r,p = test_lineartrend(x1,y1,x2,y2)
		r,p = test_lineartrend(x2,y2,x1,y1)
		print "\t".join(map(str,a + [r,p]))
	except ValueError:
		1
'
	tmpd=`mymktempd`;	
	echo "$cmd" > $tmpd/cmd
	cat $1 | python $tmpd/cmd 
	rm -rf $tmpd
	#| awk -v OFS="\t" 'BEGIN{print "id","x1","y1","x2","y2","r","pval";}{ print $0;}' \
}


fisher_test(){
## input: id x y nx ny 
## output: id x y nx ny or pvalue
##      or = x/(x+nx) * (y+ny)/y
cmd='
	con = file("stdin","r")
	out = file("stdout","w");
	line =readLines(con,n=1);
	while( length(line) > 0){ 
		tmp=strsplit(line,"\t")[[1]];
		x=unlist(lapply(tmp[2:length(tmp)], as.numeric));
		r=fisher.test(matrix(x,nrow=2))
		or=r$estimate[[1]];
		p=r$p.value;
		oline=sprintf("%s\t%g\t%g",line,or,p);
		writeLines(oline,out);
		line =readLines(con,n=1);
	}
	close(out);
	close(con);
'
	cat $1 | run_R "$cmd"
}
fisher_test2(){
	local tmpd=`mymktempd`;
#	local tmpd=tmpd;mkdir -p $tmpd;	
	cat $1 > $tmpd/a	
	awk -v OFS="\t" '{ print 1,$2,$3,$4,$5;}' $tmpd/a | sort -u > $tmpd/b 
	echo "N= "`cat $tmpd/b | wc -l`>&2;

	
	fisher_test $tmpd/b \
	| awk -v OFS="," '{ print $2,$3,$4,$5"\t"$6,$7;}' > $tmpd/ta

	cat $tmpd/a | awk -v OFS="," '{print $1"\t"$2,$3,$4,$5;}' | sort -k2,2 \
	| join -o 1.1,0,2.2 -1 2 -2 1 - $tmpd/ta \
	| perl -ne 'chomp; my @a=split/ /,$_; $a[1]=~s/,/\t/g;$a[2]=~s/,/\t/g;
		print join("\t",@a),"\n";'
	rm -rf $tmpd
}
test_fisher_test(){
echo \
'one	1	20	10	2
two	20	1	20	100'\
| fisher_test - > obs

echo \
'one	1	20	10	2	0.0134052	7.22344e-06
two	20	1	20	100	95.7161	2.53047e-12' > exp
	check obs exp
	rm obs exp
}

readline(){
cmd='
	con = file("stdin","r")
	line =readLines(con,n=1);
	while( length(line) > 0){ 
		tmp=strsplit(line,"\t")[[1]];
		g=strsplit(tmp[2],",")[[1]]
		x=as.numeric(strsplit(tmp[3],",")[[1]])
		y=as.numeric(strsplit(tmp[4],",")[[1]])
		line =readLines(con,n=1);
	}
	close(con);
'
	tmpd=`make_tempdir`
	echo "$cmd" > $tmpd/cmd
	cat $1 | R --no-save -f $tmpd/cmd
	cat $tmpd/out
	rm -rf $tmpd
}

test(){
	echo "testing test_lineartrend .. "
	echo \
	"id1	1,2,3	10,20,30	1,2,3	300,200,100" > inp
	test_lineartrend inp > obs	
	echo \
"id	x1	y1	x2	y2	r	pval	fdr
id1	1,2,3	10,20,30	1,2,3	300,200,100	-0.249029122546	1.62848179386e-10	1.62848179386e-10" > exp
	check exp obs	
	rm inp exp obs

echo "testing test_fisherexact .. "
echo \
"id	group	ctr_count	trt_count
id1	1	10	200
id2	1	20	100
id3	2	10	20" | test_fisherexact - > obs
cat obs
echo \
'id	group	ctr_count	trt_count	ctr_total	trt_total	log2fc	pval	fdr
id1	1	10	200	30	300	0.933212908788798	0.000524699917979567	0.000549410715429099
id2	1	20	100	30	300	-1.02842840832652	0.000524699917979567	0.000549410715429099
id3	2	10	20	40	320	-2.03476541816068	0.000549410715429099	0.000549410715429099' > exp
check exp obs
rm exp obs
}
