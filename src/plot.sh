. $HMHOME/src/root.sh;

plot.scatter(){
## accept id  exp1 exp2 
usage="
	$FUNCNAME <data> <out.png>
"
if [ $# -ne 2 ];then
	echo "$usage"; return;
fi

        cat $1 | run_R '
        library(MASS);
        library(ggplot2);
        tt=read.table("stdin",header=T);
        xlab=colnames(tt)[2];
        ylab=colnames(tt)[3];
        print(ylab);
        tt=tt[apply(tt[,-1],1,function(x){ sum(is.na(x))==0;}),];
        tt[ tt < 0]=0;
        x=tt[,2];y=tt[,3];
        DF= data.frame(x,y);
        dens = kde2d(x,y);
        gr = data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z));
        names(gr) = c("xgr", "ygr", "zgr");
        mod = loess(zgr~xgr*ygr, data=gr);
        DF$pointdens = predict(mod, newdata=data.frame(xgr=x, ygr=y));
        p=ggplot(DF, aes(x=x,y=y, color=pointdens)) + geom_point() +  xlab(xlab) + ylab(ylab);
	p=p + geom_hline(aes(yintercept=median(y)))
	p=p + geom_vline(aes(xintercept=median(x)))
        png("'$2'");
        print(p);
        dev.off();
        '
}
