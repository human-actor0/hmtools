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

genetrack(){
cmd='
	genes=read.table("INPUT",header=F);
	a=strsplit("GRANGE",":")[[1]];
	chrom=a[1]; grange=as.integer(strsplit(a[2],"-")[[1]]);

	gene_track.n = nrow(genes);
	gene_track.top = gene_track.n;
	gene_track.bot = 0;
	gene_track.h = ((gene_track.top-gene_track.bot)/gene_track.n);
	
	ylim=c(gene_track.bot, gene_track.top);
	xlim=grange;
	png("OUTPUT",width =480, height=200);
	plot(NULL,main="GRANGE",ylim=ylim,xlim=xlim);
	for( i in 1:gene_track.n){
		s=genes[i,2];
		e=genes[i,3];
		st=genes[i,6];
		ybot=gene_track.top - i*gene_track.h;
		ytop=ybot+gene_track.h/2;
		lengths=as.integer(strsplit(as.character(genes[i,11]),",")[[1]])
		starts=as.integer(strsplit(as.character(genes[i,12]),",")[[1]])
		for (j in 1:genes[i,10] ){
			start=s+starts[j];	
			end=start+lengths[j];
			rect(start,ybot,end,ytop);
			if (j < genes[i,10]){
				if ( st == "+"){
					lines(c(end, starts[j+1] +s),c(ybot,ybot));
				}else{
					lines(c(end, starts[j+1] +s),c(ytop,ytop));
				}
			}
		}
	}
	dev.off();
'
	tmpd=`make_tempdir`;
	mycat $2 > $tmpd/a
	cmd=${cmd//GRANGE/$1};
	cmd=${cmd//INPUT/$tmpd/a};
	cmd=${cmd//OUTPUT/$3};
	echo "$cmd";
	echo "$cmd" | R --no-save -q
	rm -rf $tmpd
}
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
	echo "$@" >&2
	interval=$1; bed=$2; out=$3;
	n=`cat $bed | wc -l`;
	tmpd=`make_tempdir`
cmd='
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
	mycat $bed | awk '$2 ~ /^[0-9]+$/' > $tmpd/a;
	if [ `cat $tmpd/a | wc -l` = "0" ];then echo "$bed is empty" >&2; return; fi
	cmd=${cmd//INTERVAL/$interval};
	cmd=${cmd/INPUT/$tmpd/a};
	cmd=${cmd/OUTPUT/$out};
	cmd=${cmd/PARAMS/${@:4}};
	echo "$cmd";
	echo "$cmd" | R --no-save -q
	rm -rf $tmpd
}

_prep_track_to_png(){
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
	interval=$1; tracks=$2;
cmd='my $interval="'$interval'";
	my ($CHROM,$START,$END);
	my ($name, $type, $color); 
	my @names=();

	if($interval=~/(\w+):(\d+)-(\d+)/){
		$CHROM=$1;$START=$2;$END=$3;
	}else{
		exit -1;	
	}
	my %data = ();
	while(<STDIN>){chomp;
		if($_=~ /^track/){
			$name="noname"; $type="bed"; $color="0,0,0";
			if($_=~/name=(\w+)/){ $name=$1; }
			if($_=~/color=(\d+),(\d+),(\d+)/){ 
				$color="rgb(".($1/255).",".($2/255).",".($3/255).")";
			}	
			if($_=~/type=(\w+)/){ $type=$1; }
	
			$data{$name}{type}=$type;
			$data{$name}{color}=$color;
			$data{$name}{starts}=();
			$data{$name}{ends}=();
			push @names,$name; ## keep the order
			next;
		}
		my @a=split /\t/,$_;
		if($type eq "bedGraph"){
			push @{$data{$name}{starts}}, $a[1];
			push @{$data{$name}{ends}}, $a[2];
			push @{$data{$name}{scores}}, $a[3];
		}elsif($type eq "bed12"){
			$data{$name}{start}=$a[1];
			$data{$name}{end}=$a[2];
			$data{$name}{name}=$a[3];
			$data{$name}{strand}=$a[5];
			$data{$name}{size}=$a[9];
			$data{$name}{starts}=$a[10];
			$data{$name}{ends}=$a[10];
		}
	}

	my $width=480; my $height=300;
	my $ntracks=scalar @names;
	print "png(\"OUT.png\",width=",$width,",height=",$height*$ntracks,")\n";
	print "par(mfrow=c(",$ntracks,",1))","\n";	

	foreach my $name (@names){
		my $color=$data{$name}{color};
		if( $data{$name}{type} eq "bedGraph"){
			print "s=c(",join( ",",@{$data{$name}{starts}}),");\n";
			print "e=c(",join( ",",@{$data{$name}{ends}}),");\n";
			print "e=e-1;"; 
			print "c=c(",join( ",",@{$data{$name}{scores}}),");\n";
			print "n=length(c); x=c(s-0.1,s,e,e+0.1); y=c(rep(0,n),c,c,rep(0,n)); ix=sort(x,index=T)\$ix; x=x[ix];y=y[ix];","\n";
			print "plot(x,y,type=\"l\",main=\"",$name,"\",xlim=c(".$START.",".$END."),col=",$color,");\n";
		}else{
		}
	}
	print "dev.off();","\n";
'
	cat $tracks | perl -e "$cmd"
#	mkdir -p $2;
#	_prep_track_to_png $1 $2
#	for f in $2/*.bedGraph;do
#		interal=`cat $2/interval.txt`
#		color=`head -n 1 $f | perl -ne 'if($_=~/color=(\d+),(\d+),(\d+)/){ print $1/255,",",$2/255,",",$3/255;}'`
#		if [ -x $color ];then color="0,0,0"; fi
#		bedgraph_to_png $interval $f ${f%.bedGraph}.png "color=rgb($color)"
#	done
}


test(){

echo \
"chr1	1	2	p1	1	+
chr1	2	3	p2	2	+
chr1	5	6	p3	7	+
chr1	7	8	p4	1	+
chr1	5	6	n1	7	-
chr1	7	8	n2	1	-" \
| polya_to_track - 'track name=mypolya1 type=bedGraph description="myPolyA 1" color=0,0,255,' > obs
echo \
'track name=mypolya1_fwd type=bedGraph description="myPolyA 1" color=0,0,255,
chr1	1	2	1
chr1	2	3	2
chr1	5	6	7
chr1	7	8	1
track name=mypolya1_bwd type=bedGraph description="myPolyA 1" color=0,0,255,
chr1	5	6	7
chr1	7	8	1' >  exp

check exp obs
track_to_png "chr:0-10" obs
rm -rf exp obs

echo \
'track name=mypolya1_fwd type=bedGraph description="myPolyA 1" color=0,0,255,
chr1	0	100	n1	0	+	10	90	0,0,0	3	10,20,30	0,20,70' \
| track_to_png "chr:0-100" -

}


