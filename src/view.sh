#!/bin/bash

. $HMHOME/src/root.sh
. $HMHOME/src/polya.sh

#datq=( $@ );
#for f in ${@:2};do
#        a=( `echo "$f" | tr ":" "\t"` )
#        odir=$O/$1/${a[0]}/${a[1]}; mkdir -p $odir
#        echo "${a[*]} =>  $odir .. ";
#        if [ ${a[0]} = "gene" ];then
#                #samtools view -b ${a[2]} $1 | bamToBed -bed12 > $odir/a.bed 
#                echo "$1" | perl -ne '
#                        if( $_ =~ /(.+):(\d+)-(\d+)/ ){ 
#                                print $1,"\t",$2,"\t",$3,"\n";
#                        }
#                ' | intersectBed -a ${a[2]} -b stdin -wa -u  > $odir/a.bed
#        elif [ ${a[0]} = "polya" ]; then
#                eval "
#                samtools view -bq $Q ${a[2]} $1 | bamToBed | modify_score - $M | point - > $odir/a.bed 
#                "
#        fi
#done
#
#}

polya_to_track(){
usage="
	$FUNCNAME <bed> <ucsc_track_head>
"; if [ $# -lt 2 ];then echo "$usage"; return; fi;
	tmpd=`make_tempdir`
	sort_bed $1 | cut -f1-6 > $tmpd/a
	awk -v OFS="\t" '{if($6=="+"){ print $1,$2,$3,$5;}}' $tmpd/a >  $tmpd/a_fwd
	awk -v OFS="\t" '{if($6=="-"){ print $1,$2,$3,$5;}}' $tmpd/a > $tmpd/a_bwd
	if [  `cat $tmpd/a_fwd | wc -l` > 0 ];then
		echo ${@:2} | perl -ne '$_=~ s/name=([\w|\"]+)/name=${1}_fwd/;print';
		cat $tmpd/a_fwd;
	fi
	if [  `cat $tmpd/a_bwd | wc -l` > 0 ];then
		echo ${@:2} | perl -ne '$_=~ s/name=([\w|\"]+)/name=${1}_bwd/;print';
		cat $tmpd/a_bwd;
	fi
	rm -rf $tmpd
}
bedgraph_to_png(){
usage="
$FUNCNAME <bed> <interval> [<Rplot_options>]
 <Rplot_options>: 
	default) width=480 height=200 name=\"noname\" color=black ylab=name
"
if [ $# -lt 2 ];then echo "$usage"; return; fi
	bed=$1; interval=$2; out=$3; ### other prameters 
cmd='
	width=480; height=200; name=\"noname\"; color=black; ylab=name;
	color="black";
	width=480; height=200;
	PARAMS ## override parameters here
	main="INTERVAL";
	a=strsplit("INTERVAL",":")[[1]]; chrom=a[1]; grange=as.integer(strsplit(a[2],"-")[[1]]);
	out="OUTPUT";

	tt=read.table("INPUT",header=F); 
	s=tt[,2];e=tt[,3]-1;c=tt[,4];n=length(c);
	x=c(s-0.1,s,e,e+0.1);
	y=c(rep(0,n),c,c,rep(0,n));
	ix=sort(x,index=T)$ix; x=x[ix];y=y[ix];

	png(out, width=width,height=height);
	plot(x,y,type="l",xlim=grange,main=main,col=color)
	dev.off();
'
	tmpd=`make_tempdir`
	mycat $bed > $tmpd/a;
	mycat $bed | awk '$2 ~ /^[0-9]+$/' > $tmpd/a; ## hmm this is smarter than tail -n+2
	if [ `cat $tmpd/a | wc -l` = "0" ];then echo "$bed is empty" >&2; return; fi

	params=${@:4};
	a=`echo $interval | tr ":-" "\t"`;
	cmd=${cmd/CHROM/${a[0]}};
	cmd=${cmd/CSTART/${a[1]}};
	cmd=${cmd/CEND/${a[2]}};
	cmd=${cmd/INPUT/$tmpd/a};
	cmd=${cmd/OUTPUT/$out};
	cmd=${cmd/PARAMS/${@:4}};
	echo "$cmd";
#	echo "$cmd" | R --no-save -q
	rm -rf $tmpd
}

split_tracks(){
cmd='use strict; 
	my %fh=(); my $name="out";
	my $odir="ODIR/"; 
	while (<STDIN>){chomp;
		if($_=~/track/){
			my $type="bed";my $name="unknown";
			if($_=~/type=(\w+)/){ $type=$1;}
			if($_=~/name=(\w+)/){ $name=$1;}
			close F unless defined fileno F;
			open(F,">",$odir.$name.".".$type) or die "$!";
		}elsif($_=~/browser position:(\w+:\d+-\d+)/){
			open(F,">",$odir."interval.txt") or die "$!";
			print F $1;
			close(F);
		}
		print F $_,"\n";
	}
	close F unless defined fileno F;
'
	cmd=${cmd/ODIR/$2};
	cat $1 | perl -e "$cmd";
}





track_to_png(){
cmd='
	getrgb=function(x){ 
		if(is.null(x)){
			return(rgb(0,0,0));
		}
		x=as.integer(strsplit(as.character(x),",")[[1]]); 
		y=x/255; 
		return(rgb(y[1],y[2],y[3])); 
	}
	draw_genetrack=function(data,chrom,cstart,cend,main){
		gene_track.n = nrow(data);
		gene_track.top = gene_track.n;
		gene_track.bot = 0;
		gene_track.h = ((gene_track.top-gene_track.bot)/gene_track.n);
		ylim=c(gene_track.bot, gene_track.top);
		plot(NULL,main=paste("GRANGE",name,sep=":"),ylim=ylim,xlim=c(cstart,cend));
		for( i in 1:gene_track.n){
			s=data[i,2]; e=data[i,3]; st=data[i,6]; color=getrgb(data[i,9]);
			ybot=gene_track.top - i*gene_track.h; ytop=ybot+gene_track.h/2;
			lengths=as.integer(strsplit(as.character(data[i,11]),",")[[1]])
			starts=as.integer(strsplit(as.character(data[i,12]),",")[[1]])
			for (j in 1:data[i,10] ){
				start=s+starts[j];	
				end=start+lengths[j];
				rect(start,ybot,end,ytop,col=color);
				if (j < data[i,10]){
					if ( st == "+"){
						lines(c(end, starts[j+1] +s),c(ybot,ybot),col=color);
					}else{
						lines(c(end, starts[j+1] +s),c(ytop,ytop),col=color);
					}
				}
			}
		}
	}
	draw_bedgraph=function(data,chrom,cstart,cend,color,main){
		s=data[,2];e=data[,3]-1;c=data[,4];n=length(c);
		x=c(s-0.1,s,e,e+0.1);
		y=c(rep(0,n),c,c,rep(0,n));
		ix=sort(x,index=T)$ix; x=x[ix];y=y[ix];
		plot(x,y,type="l",xlim=c(cstart,cend),col=color,main=main)
	}


	getv=function(k,x){ ## return value from key=value pair
		i=grep(k,x);
		if(length(i)==1){ return(gsub(k,"",x[i])); }
		return(NULL);
	}
	data=list();
	con=file("INPUT","r");
	line=readLines(con,n=1);
	chrom=""; cstart=-1; cend=-1;
	while(length(line) > 0){
		if(line == ""){ line=readLines(con,n=1);next ;}
		tmp=strsplit(line," ")[[1]]
		## handle track
		if( tmp[1] == "track"){
			name=getv("name=",tmp); 
			data[[ name ]] = list(head=line,bed=""); 
		}else if( tmp[1] == "browser"){
			a=strsplit(tmp[3],":")[[1]];
			chrom=a[1];
			b=strsplit(a[2],"-")[[1]];
			cstart=as.integer(b[1]);
			cend=as.integer(b[2]);
		}else{
			data[[ name ]]$bed = paste(data[[ name ]]$bed,line,sep="\n");
		}
		line=readLines(con,n=1);
	}
	close(con);
	png("OUT");
	par(mfrow=c(length(data),1));
	for( n in names(data)){
		tmp=strsplit(data[[n]]$head," ")[[1]];
		type=getv("type=",tmp);
		color=getrgb(getv("color=",tmp));
		if(type=="bed"){
			draw_genetrack(read.table(text=data[[n]]$bed),chrom,cstart,cend,n);
		}else if(type=="bedGraph"){
			draw_bedgraph(read.table(text=data[[n]]$bed),chrom,cstart,cend,color,n);
		}
	}
	dev.off();

'
	tmpd=`make_tempdir`
	mycat $1 > $tmpd/a
	cmd=${cmd/INPUT/$tmpd/a};
	#a=( `echo $2 | tr ":-" "\t"` );
	#cmd=${cmd/CHROM/${a[0]}}; cmd=${cmd/CSTART/${a[1]}}; cmd=${cmd/CEND/${a[2]}};
	cmd=${cmd/OUT/$2};
	echo "$cmd" | R --no-save 
	rm -rf $tmpd
}

test(){

#echo \
#"chr1	1	2	p1	1	+
#chr1	2	3	p2	2	+
#chr1	5	6	p3	7	+
#chr1	7	8	p4	1	+
#chr1	5	6	n1	7	-
#chr1	7	8	n2	1	-" \
#| polya_to_track - 'track name=mypolya1 type=bedGraph description="myPolyA 1" color=0,0,255,' > obs
echo \
'
browser position chr1:0-100
track name=mypolya1_fwd type=bedGraph description="myPolyA 1" color=0,0,255,
chr1	1	2	1
chr1	2	3	2
chr1	5	6	7
chr1	7	8	1
track name=mypolya1_bwd type=bedGraph description="myPolyA 1" color=0,0,255,
chr1	5	6	7
chr1	7	8	1
track name=gene1 type=bed description="gene 1" color=0,0,255,
chr1	0	100	n1	0	+	10	90	0,0,0	3	10,20,30	0,20,70
chr1	0	100	n2	0	-	10	90	0,0,255	3	10,20,30	0,20,70' \
| track_to_png - out.png

#track_to exp chr1:0-10 out
#
#check exp obs
#track_to_png "chr:0-10" obs
#rm -rf exp obs

echo \


}


