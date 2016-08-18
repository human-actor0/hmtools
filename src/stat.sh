#!/bin/bash
. $HMHOME/src/root.sh

stat.prop_test(){
usage(){
	echo "
	$FUNCNAME <txt> <header>
	<header> : T or F
	"
}; if [ $# -lt 1 ];then usage;return; fi
	cat $1 | run_R '
		H='${2:-"F"}';
		tt=read.table("stdin",header=H);
		if(!H){ colnames(tt)=c("ID","CTR.cnt","TRT.cnt","CTR.tot","TRT.tot"); }
		# global total counts
		pval=apply(tt[,-1],1,function(x){
			p=1;
			if( sum(x > 0) == 4 && sum(x[3:4]==x[1:2]) < 2){
				p=prop.test(x[1:2],n=x[3:4])$p.value;
			}
			return(p);
		})

		relFC=apply(tt[,-1],1,function(x){
			x1=x+0.5;##pseudo counting
			return(log2( x1[2]*x1[3]/(x1[1]*x1[4])));
		})
		#totFC=log2( (tt[,5]+0.5)*sum(tt$CTR.N)/(tt[,4]+0.5)/sum(tt$TRT.N));
		tt[[ "log2_TRT/CTR" ]]=relFC;
		#tt$log2TFC=totFC;
		tt$pval=pval;
		write.table(tt,file="stdout",col.names=T,row.names=F,quote=F,sep="\t");
	'
}
stat.prop_test.test(){
echo  \
"ID	CTR.cnt	TRT.cnt	CTR.tot	TRT.tot
A	1	2	10	30
B	0	2	10	30
C	1	2	1	2
D	0	0	2	10
E	1	2	1	19" \
| stat.prop_test - T  
}
stat.edger_test.test(){
echo \
"a	1	10	11	10	11	10
b	2	20	22	20	22	20
c	3	30	33	30	33	30" | stat.edger_test -
}

stat.lineartrend(){
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
		#r,p = test_lineartrend(x2,y2,x1,y1)
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
stat.lineartrend.test(){
echo	"id1	1,2,3	10,20,30	1,2,3	300,200,100" \
| stat.lineartrend - 
}

stat.rfuncs(){
echo '
	calc_ks=function(x,y){
		x[is.na(x)]=0; y[is.na(y)]=0;
		if( sum(x)==0 || sum(y)==0){ return(0);}	
		x=cumsum(x);x.n=x[length(x)];
		y=cumsum(y);y.n=y[length(y)]; 
		z=max(abs(x/x.n - y / y.n));
		return(z);
	}
	sam=function(x,n){
		return(matrix( unlist(lapply(x,function(x){ rpois(n,x);})), nrow=n,byrow=F));
	}
	sim=function(x,y,n){
		x.s=sam(x,n);
		y.s=sam(y,n);
		apply(cbind(x.s,y.s),1,function(x){
			n=length(x);	
			ks=calc_ks(x[1:(n/2)],x[(n/2+1):n]);
			return(ks);
		})
	}
	ks_sim=function(x,y,N){
		ks0=calc_ks(x,y);
		a=sum(x)/(sum(x)+sum(y));
		ks1=sim((x+y)*a,(x+y)*(1-a),N); ## null model
		pval=(sum( ks1 >= ks0)/N); ## false positives
		return(list(ks=ks0,pval=pval,power=NA));
	}
	
	mymerge=function(files,names,FUN){
		res=NULL;
		for(i in 1:length(files)){
			tt=read.table(files[i],header=T);
			if(!is.null(FUN)){
				tt=FUN(tt);
			}
			cn=colnames(tt);
			colnames(tt)=c(cn[1], paste(names[i],cn[2:length(cn)],sep="."));
			if(is.null(res)){       
				res=tt;
			}else{
				res=merge(res,tt,by=1,all=T);
			}
		}
		return(res);    
	}

'
}

stat.bc(){
usage="
FUNCT: compare bin count patterns
USAGE: $FUNCNAME <ctr.file> <trt.file>
	*.file : three columns (start end count)
"
	cmd=`stat.rfuncs`'
		method="'$2'"
		N=10000;
		tt=read.table("stdin",header=T);
		d=tt[sort(tt$start,decreasing=F,index=T)$ix,];
		res=ks_sim(d$count.1,d$count.2,N)
		write.table(file="stdout",as.data.frame(res),col.names=F,row.names=F,sep="\t",quote=F);
	' 
	run_R "$cmd" 
}
stat.bc.test(){
echo \
"start	end	count.1	count.2
1	2	1	3
3	4	2	2
5	6	3	1" | stat.bc - "ks"
echo \
"start	end	count.1	count.2
1	2	10	30
3	4	20	20
5	6	30	10" | stat.bc - "ks"
echo \
"start	end	count.1	count.2
1	2	1	2
3	4	2	3
5	6	3	4" | stat.bc - "ks"
echo \
"start	end	count.1	count.2
1	2	10	20
3	4	20	30
5	6	30	40" | stat.bc - "ks"
echo \
"start	end	count.1	count.2
1	2	100	200
3	4	200	300
5	6	300	400" | stat.bc - "ks"

}

stat.sum(){
	cat $1 | perl -e 'use strict; my %res=(); my $H='${2:-"0"}';
	while(<STDIN>){ chomp; my @a=split/\t/,$_;
		if($H > 0 ){ print $_,"\n"; $H--;next;}
		if(!defined $res{$a[0]}){
			$res{$a[0]} = ();
		}	
		for(my $i=0; $i< $#a; $i++){
			$res{ $a[0] }[$i] += $a[$i+1]; 
		}
	}
	foreach my $k (keys %res){
		print $k,"\t",join("\t",@{$res{$k}}),"\n";
	}'
}
stat.sum.test(){
echo \
"a	1	2
b	3	4
a	10	20
b	30	40" | stat.sum - 
}

stat.merge(){
	cmd='files=c('`quote $@ | tr " " ","`');
	d=NULL;
	for( f in files){
		tt=read.table(f,header=T);
		if(is.null(d)){
			d=tt;
		}else{
			d=merge(d,tt,by=1,all=T);
		}
	}
	d[is.na(d)]=0;
	write.table(d,file="stdout",col.names=T,row.names=F,quote=F,sep="\t");
	'
	run_R "$cmd";

}

stat.prep(){
usage="
FUNCT: make a merged table with headers 
USAGE: $FUNCNAME <ctr1>[,<ctr2>..] <trt1>[,<trt2 ..] [log]
"; if [ $# -lt 2 ];then echo "$usage"; return; fi
	IFS=,; local trts=`quote $1`; local ctrs=`quote $2`; 
	IFS=$' \t\n';
        cmd='
                trt=c('"$trts"'); ctr=c('"$ctrs"');
                d=NULL;
                for( i in 1:length(trt)){
                        tt=read.table(trt[i],header=F);
			j=ncol(tt)-1;
                        colnames(tt)=c("id",paste(paste("ctr",i,sep=""),".c",1:j,sep="")); 
                        if( is.null(d)){ d=tt;
                        }else{ d=merge(d,tt,by="id",all=T); }
                }
                for( i in 1:length(ctr)){
                        tt=read.table(ctr[i],header=F);
                        colnames(tt)=c("id",paste(paste("trt",i,sep=""),".c",1:j,sep="")); 
                        d=merge(d,tt,by="id",all=T);
                }
		d[ is.na(d) ] = 0;
		write.table(file="stdout",d, col.names=T,quote=F,sep="\t", row.names=F);
        '
	run_R "$cmd"  ${3:-""} 
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
stat.gix2gixx(){
usage(){ echo "
 $FUNCNAME <file> 
  <file> : group_id id count1 [count2 ..]
"
}
if [ $# -ne 1 ];then  usage; return; fi
cat $1 | perl -e 'use strict;
        my %d=(); my %dt=();
        while(<STDIN>){ chomp;
                my ($g,$i,@x) = split /\t/,$_;
		if( !defined $d{$g} ){
			@{$d{$g}}=();
			@{$dt{$g}}=();
		}
                push @{$d{$g}}, [$g, $i, @x];
		for( my $i=0; $i < scalar @x; $i++){
			$dt{$g}[$i] += $x[$i];
		}
        }
        foreach my $g (keys %d){
	my @b=@{$dt{$g}};
	foreach my $ae (@{$d{$g}}){
		my ($gi,$id,@x)=@$ae; 
		print $gi,"\t",$id,"\t",join("\t",@x),"\t",join("\t",@b),"\n";
        }}
'
}
stat.gix2gixx.test(){
echo \
"g1	id1	1	11
g1	id2	2	22
g2	id3	3	33" \
| stat.gix2ixx - 
}
stat.gix2gixnx(){
usage="
FUNCT: Group, Id, Count to Group Id Count GroupCount - Count
USAGE: $FUNCTNAME <file>
"
cat $1 | perl -e 'use strict;
        my %d=(); my %dt=();
        while(<STDIN>){ chomp;
                my ($g,$i,@x) = split /\t/,$_;
		if( !defined $d{$g} ){
			@{$d{$g}}=();
			@{$dt{$g}}=();
		}
                push @{$d{$g}}, [$g, $i, @x];
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
stat.gix2gixnx.test(){
echo \
"g1	1	1	11
g1	2	2	22
g2	3	3	33" \
| stat.gix2gixnx - > obs
echo \
"g2	3	3	33	0	0
g1	1	1	11	2	22
g1	2	2	22	1	11" > exp
check obs exp
rm -rf obs exp
}

stat.de(){
	data=`echo "$@" | quote`;
	cmd='files=c('"$data"');
		d=NULL;
		for( f in files ){
			tt=read.table(f,header=F,stringsAsFactors=F);
			colnames(tt)=c("id",paste("DE.",f,sep=""), paste("X.",f,sep=""));

			if(is.null(d)){ d=tt;
			}else{ d=merge(d,tt,by=c(1),all=T); }
		}
		colnames(d)=gsub("tmpd\\/","",colnames(d));
		d1=d[,grep("DE.",colnames(d))]; d1[ is.na(d1) ] = "I";
		d2=d[,grep("X.",colnames(d))];

		write.table(cbind(d$id,d1,d2),file="stdout",col.names=T,row.names=F,quote=F,sep="\t");
		cb=t(combn(2:ncol(tt),2));

	'
	run_R "$cmd" 
}
stat.de.test(){
echo \
"a	U	0
b	P	2
c	N	0.5" > a
echo \
"a	P	4
b	P	4" > b
stat.de a b
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

stat.cor(){
usage(){ echo \
" $FUNCNAME <file> [<method>]
	<file>: two columns of values
	<method>: pearson, kendall, spearman (default)
"
}
if [ $# -lt 1 ];then usage; return; fi
	local m=${2:-spearman};
	cat $1 | run_R ' 
		tt=read.table("stdin",header=F);
		cor(tt[,1],tt[,2],method="'$m'" );' \
	|perl -ne 'chomp;if($_=~/\[1\] ([\-|\d|\.]+)/){print $1,"\n";}'
}

stat.cor.test(){
echo \
"1	2
2	4
3	5"  | stat.cor -
}

stat.hclustwcor(){
usage(){
echo " $FUNCNAME <file> <out>
	<file> : id1 id2 correlation

"
}
	if [ $# -ne  2 ];then usage; return; fi
	local tmpd=`mymktempd`;
	cmd='
		tt=read.table("stdin",header=F);
		x=as.character(unique(c(as.character(tt[,1]),as.character(tt[,2]))));
		tt;
		#tt=rbind(tt,data.frame(V1=x,V2=x,V3=1));
		tt[,3]=1-tt[,3];
		d=as.dist(xtabs(tt[, 3] ~ tt[,2] + tt[,1]));
		d;
		png("'$2'"); plot(hclust(d));dev.off();
	'
	cat $1 > $tmpd/i
	cut -f 1 $tmpd/i  > $tmpd/a
	cut -f 2 $tmpd/i >> $tmpd/a
	for f in `sort -u $tmpd/a`;do
		echo -e "$f\t$f\t1" >> $tmpd/b;
	done
	head $tmpd/*
	cat $tmpd/i $tmpd/b | sort -k1,2| run_R "$cmd"
	rm -rf $tmpd
}
stat.hclustwcor.test(){
echo \
"A	B	1
A	C	0.5
B	C	0.3" | stat.hclustwcor - tmp.png 
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

stat.edger(){
usage="
FUNCT: calculate p-values 
USAGE: $FUNCNAME <stat.prep> [<bcv>]
"
	## INPUT: comma separated control and treatment EI file ( bed6 + x + nx )
	## OUTPUT: bed6 + logFC + pvalue 
	local BCV=${2:-0.1}; # 0.4: human data, 0.1: genetically identical, 0.01: technical replicates (edgeRUsersGuide)
	cat $1 | run_R '
	tt=read.table("stdin",header=T);
	cn=colnames(tt)[2:ncol(tt)];
	group=rep(1,length(cn)); group[ grep("^trt",cn)]=2;
	event=rep(1,length(cn)); event[ grep(".c2",cn)]=2;
	library(edgeR)
	if( length(grep("trt",cn)) == 1 || length(grep("ctr",cn)) == 1){
		bcv='"$BCV"';
		y = DGEList(counts=tt[,2:3],group=c(2,1))
		et = exactTest(y, dispersion=bcv^2)	
		write.table(file="stdout",cbind(tt,et$table),quote=F,sep="\t",row.names=F);
	}else{
	
	#ix=apply(D[,7:ncol(D)], 1, min) > 0 & apply(D[,7:ncol(D)],1,max) > 10
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
	}' 
}
stat.edger.test(){
echo \
"id	trt1.c1	trt1.c2	trt2.c1	trt2.c2	ctr1.c1	ctr1.c2
a	1	10	11	10	11	10
b	2	20	22	20	22	20
c	3	30	33	30	33	30" | stat.edger -

echo \
"id	trt1.c1	trt1.c2	ctr1.c1	ctr1.c2
a	1	10	11	10
b	2	20	22	20
c	3	30	33	30" | stat.edger -
}



stat.padjust(){
usage="
usage: $FUNCNAME <file> 
"
if [ $# -ne 1 ]; then echo "$usage"; return; fi
cmd='
	tt=read.table("stdin",header=T);
	tt$fdr=p.adjust(tt$pval,method="fdr");
        write.table(file="stdout",tt,row.names=F,col.names=T,quote=F,sep="\t");
'
	cat $1 | run_R "$cmd"
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
