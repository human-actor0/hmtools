#!/bin/bash
. $HMHOME/src/root.sh

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
        write.table(file="stdout",tt,row.names=F,col.names=T,quote=F,sep="\t");
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

        #G=G[gs>1]; tt=tt[gs>1,];
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
## input: id x  nx  y ny
## output: id x nx y ny or pvalue
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
"id1	1	10	200
id2	1	20	100
id3	2	10	20" | test_fisherexact - > obs
echo \
'id	group	ctr_count	trt_count	ctr_total	trt_total	log2fc	pval	fdr
id1	1	10	200	30	300	0.933212908788798	0.000524699917979567	0.000549410715429099
id2	1	20	100	30	300	-1.02842840832652	0.000524699917979567	0.000549410715429099
id3	2	10	20	40	320	-2.03476541816068	0.000549410715429099	0.000549410715429099' > exp
check exp obs
rm exp obs
}
