#!/bin/bash
. $HMHOME/src/root.sh

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
	cat $1 | python $tmpd/cmd

}

test_fisherexact(){
GRS=2; CLS="3,4";
cmd='
        cols=c(CLS);
        grp=c(GRS);
        tt=read.table("stdin",header=F);
        if(length(grp) == 1){ G=tt[,grp];
        }else{
                G=apply(tt[,grp],1,function(x){ paste(x,collapse="|");});
        }       
        gs=ave(1:length(G),G,FUN=length); ## group sum
        G=G[gs>1]; tt=tt[gs>1,];
        m=tt[,cols];
        M=apply(m,2,function(x){ ave(x,G,FUN=sum)}) ## group sum
        p=unlist(apply(cbind(m,M-m),1,function(x){ fisher.test(matrix(x,byrow=F,nrow=2))$p.value }))
        fdr=p.adjust(p,method="fdr")
        log2fc= log2((0.5+m[,2])/(m[,1]+0.5)*M[,1]/M[,2]);
        tt$log2FC=log2fc; tt$pval=p; tt$FDR=fdr;
        write.table(file="stdout",tt,row.names=F,col.names=T,quote=F,sep="\t");
'
	tmpd=`make_tempdir`
	cmd=${cmd//GRS/$GRS};
	cmd=${cmd//CLS/$CLS};
	cmd=${cmd//stdout/$tmpd/out};
	echo "$cmd" > $tmpd/cmd
	cat $1 | R --no-save -f $tmpd/cmd &> $tmpd/log
	cat $tmpd/out
	rm -rf $tmpd
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
	echo \
	"id1	1,2,3	10,20,30	1,2,3	300,200,100	-0.249029122546	1.62848179386e-10" > exp
	test_lineartrend inp > obs	
	check exp obs	
	rm inp exp obs
	echo "testing test_fisherexact .. "
	echo \
"id1	1	10	200
id2	1	20	100" | test_fisherexact -
}
