. $HMHOME/src/root.sh;
edger.rep(){
usage(){ echo "
$FUNCNAME <table> <ctr_col> <trt_col> [option]
 [option]
	nooptions 
"
}
	if [ $# -lt 3 ];then usage;return; fi
	cat $1 | run_R '
	library(edgeR)
	CTR="^'$2'"; TRT="^'$3'"; 
	tt=read.table("stdin",header=T,check.names=F);
	cn=colnames(tt);
	m=tt[,c(grep(CTR,cn), grep(TRT,cn))];

	m.cn=colnames(m);
	group=rep(1,length(m.cn)); 
	group[ grep(CTR,m.cn)]=1;
	group[ grep(TRT,m.cn)]=2;
	coef=2;
	
	y=DGEList(counts=m,group=factor(group));
	y=calcNormFactors(y);

	group.this=factor(group);
	H1 <- model.matrix(~ group.this )
	y = estimateGLMCommonDisp(y,H1);
	y = estimateGLMTrendedDisp(y,H1);
	y = estimateGLMTagwiseDisp(y,H1);

	fit=glmFit(y$counts, H1, y$tagwise.dispersion,offset=0,prior.count=0)
	llh=glmLRT(fit,coef=coef)
	nn=paste("logFC",TRT,"/",CTR,sep="")
	tt[[nn]]=llh$table$logFC;
	tt[["pval"]]=llh$table$PValue;
	##res=data.frame(tt, nn=llh$table$logFC, pval=llh$table$PValue)
	## chrom start end logFC pval
	write.table(tt,file="stdout", col.names=T,quote=F,row.names=F,sep="\t");
	' 
}
edger.rep.test(){
echo \
"x@trt@ctr	trt1	trt2	ctr1	ctr2
a	1	2	10	20
b	4	0	50	0" \
| edger.rep - ctr trt


}
edger.interact(){
usage(){ echo "
$FUNCNAME <table> <ctr_col> <trt_col> [option]
"
}
	if [ $# -lt 3 ];then usage;return; fi
	cat $1 | run_R '
	library(edgeR)
	CTR="^'$2'"; TRT="^'$3'"; 
	tt=read.table("stdin",header=T,check.names=F);
	cn=colnames(tt);
	m=tt[,c(grep(CTR,cn), grep(TRT,cn))];

	m.cn=colnames(m);
	group=rep(1,length(m.cn)); 
	group[ grep(CTR,m.cn)]=1;
	group[ grep(TRT,m.cn)]=2;
	event=rep(1,length(m.cn)); event[ grep(".c2",m.cn)]=2;
	
	y=DGEList(counts=m,group=factor(group));
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
	nn=paste("logFC",TRT,"/",CTR,sep="")
	tt[[nn]]=llh$table$logFC;
	tt[["pval"]]=llh$table$PValue;
	##res=data.frame(tt, nn=llh$table$logFC, pval=llh$table$PValue)
	## chrom start end logFC pval
	write.table(tt,file="stdout", col.names=T,quote=F,row.names=F,sep="\t");
	' 
}
edger.interact.test(){
echo \
"x@trt@ctr	trt1.c1	trt1.c2	trc.c1	trt2.c1	trt2.c2	ctr1.c1	ctr1.c2
a	1	10	NA	11	10	11	10
b	2	20	NA	22	20	22	20
c	3	30	NA	33	30	33	30" | edger.interact - ctr trt log

}




edger.norep(){
if [ $# -lt 2 ];then
	echo \
"USAGE $FUNCNAME <txt> <BCV> [<header>]
	<BCV>: 0.4: human data, 0.1: genetically identical, 0.01: technical replicates (edgeRUsersGuide)
"
return;
fi
	cat $1 | run_R 'H='${3:-"F"}';
		library(edgeR);
		tt=read.table("stdin",header=H);
		if(!H){ colnames(tt)=c("ID","CTR","TRT");}
		bcv = '$2';
		y = DGEList(counts=tt[,2:3],group=c(1,2))
		y = calcNormFactors(y);
		et = exactTest(y, dispersion=bcv^2)	
		write.table(file="stdout",cbind(tt,et$table),quote=F,sep="\t",row.names=F);
	'
}
edger.norep.test(){
echo \
"a	2	3
b	4	19" | edger.norep - 0.2
}
