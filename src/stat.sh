#!/bin/bash
. $HMHOME/src/root.sh

max(){
	awk 'NR == 1 { max=$1;} { if($1 > max) max=$1;} END { print max;}'	
}
min(){
	awk 'NR == 1 { min=$1;} { if($1 < min) min=$1;} END { print min;}'	
}

cor(){
        cat $1 | R --no-save -q -e 'tt=read.table("stdin",header=F);cor(tt[,1],tt[,2],method="spearman" );' \
        | perl -ne 'chomp;if($_=~/\[1\] ([\d|\.]+)/){print $1;}'
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

test_edger(){
	## INPUT: comma separated control and treatment EI file ( bed6 + x + nx )
	## OUTPUT: bed6 + logFC + pvalue 
	local rcmd='
	tt=read.table("stdin",header=T);
	cn=colnames(tt)[3:ncol(tt)];
	group=rep(1,length(cn)); group[ grep("t.",cn)]=2;
	event=rep(1,length(cn)); event[ grep(".x",cn)]=2;
	
		
	#ix=apply(D[,7:ncol(D)], 1, min) > 0 & apply(D[,7:ncol(D)],1,max) > 10
	library(edgeR)
	y=DGEList(counts=tt[,3:ncol(tt)],group=factor(group));
	y=calcNormFactors(y);

	event.this=factor(event);
	group.this=factor(group);
	H1 <- model.matrix(~ event.this + group.this + event.this:group.this, data=y$samples )
	H0 <- model.matrix(~ event.this + group.this )
	print (H1)
	coef <- (ncol(H0)+1):ncol(H1)
	#y=estimateCommonDisp(y)
	#y=estimateTagwiseDisp(y, trend="movingave")
	y = estimateGLMCommonDisp(y,H1);
	y = estimateGLMTrendedDisp(y,H1);
	y = estimateGLMTagwiseDisp(y,H1);

	fit=glmFit(y$counts, H1, y$tagwise.dispersion,offset=0,prior.count=0)
	llh=glmLRT(fit,coef=coef)

	ex.h0=apply( y$counts[,group.this == 1 & event.this == 1], 1, sum);
	in.h0=apply( y$counts[,group.this == 1 & event.this == 2], 1, sum);

	out="OUT";
	res=data.frame(tt, logIR=log( in.h0/ ex.h0), logFC=llh$table$logFC, pval=llh$table$PValue)
	## chrom start end logFC pval
	write.table(res, out, col.names=T,row.names=F,sep="\t",quote=F);
	'

	local tmpd=`make_tempdir`;
	rcmd=${rcmd/OUT/$tmpd/out}
	echo "$rcmd" > $tmpd/a
	cat $1 | R --no-save -f $tmpd/a >&2
	cat $tmpd/out
	rm -rf $tmpd
}
test__test_edger(){
echo \
"id	group	c.x	c.x	t.x	c.n	c.n	t.n
id1	g1	1	1	20	3	30	20
id2	g2	2	1	22	4	31	20
id4	g2	4	1	24	6	31	20" \
| test_edger - 
}
#test__test_edger

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
test_lineartrend(){
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
		r,p = test_lineartrend(x1,y1,x2,y2)
		print "\t".join(map(str, a + [r,p]))
	except ValueError:
		1
'
	tmpd=`make_tempdir`;	
	echo "$cmd" > $tmpd/cmd
	cat $1 | python $tmpd/cmd \
	| padjust - -1
	rm -rf $tmpd

	#| awk -v OFS="\t" 'BEGIN{print "id","x1","y1","x2","y2","r","pval";}{ print $0;}' \
}

test_fisherexact(){
usage="
usage: $FUNCNAME <file> [pseudo_zero]
 <file>: columns of id, group, ctr_count, trt_count
"
cmd='
        tt=read.table("stdin",header=T);
	G=tt$group;
        gs=ave(1:length(G),G,FUN=length); ## group sum

        G=G[gs>1]; tt=tt[gs>1,];
        m=cbind(tt$ctr_count,tt$trt_count);
        M=apply(m,2,function(x){ ave(x,G,FUN=sum)}) ## group sum

	PM= apply(m,2,sum); M[gs==1,1]=PM[1]; M[gs==1,2]=PM[2];

        pval=unlist(apply(cbind(m,M-m),1,function(x){ fisher.test(matrix(x,byrow=F,nrow=2))$p.value }))
        fdr=p.adjust(pval,method="fdr")
        log2fc= log2((0.5+m[,2])/(m[,1]+0.5)*M[,1]/M[,2]);

	tt$ctr_total=M[,1];
	tt$trt_total=M[,2];
        tt$log2fc=log2fc; tt$pval=pval; tt$fdr=fdr;
        write.table(file="stdout",tt,row.names=F,col.names=T,quote=F,sep="\t");
'
	
	cmd=${cmd/EPS/${2:-0}};
	cat $1 \
	| awk -v OFS="\t" 'BEGIN{print "id","group","ctr_count","trt_count";}{ print $0;}' \
	| run_R "$cmd"
	#tmpd=`make_tempdir`
	#cmd=${cmd//stdout/$tmpd/out};
	#echo "$cmd" > $tmpd/cmd
	#cat $1 | R --no-save -f $tmpd/cmd &> $tmpd/log
	#cat $tmpd/out
	#rm -rf $tmpd
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
