#!/bin/bash
. $HMHOME/src/root.sh
. $HMHOME/src/bed.sh

ucsc.chrom_size(){
if [ $1 = "hg19" ];then
echo \
"chr1	249250621
chr2	243199373
chr3	198022430
chr4	191154276
chr5	180915260
chr6	171115067
chr7	159138663
chr8	146364022
chr9	141213431
chr10	135534747
chr11	135006516
chr12	133851895
chr13	115169878
chr14	107349540
chr15	102531392
chr16	90354753
chr17	81195210
chr18	78077248
chr19	59128983
chr20	63025520
chr21	48129895
chr22	51304566
chrX	155270560
chrY	59373566
chrM	16571"
else
	echo "no such $1 genome" >&2
fi
}


ucsc.ens2genename(){
usage=" $FUNCNAME <file>
"
if [ $# -ne 1 ];then echo "$usage";return; fi

        cat $1 | perl -e 'use strict;
                my %eg=();
                my $file=$ARGV[0];
                open(F,$file) or die;
                while(<F>){ my ($k,$v)=split /\s/,$_; $eg{$k}=$v; }
                close(F);
                while(<STDIN>){ chomp;
                        if($_=~/(ENSG\d+)/){
                                my $tmp=$eg{$1};
                                if(defined $tmp){
                                        $_=~s/ENSG\d+/$tmp/g;
                                }
                        }
                        print $_,"\n";
                }
        ' $HMHOME/data/hg19/ensgToGenename.txt
}
ucsc.bg2bw(){
usage="
USAGE: $FUNCNAME <chrom.size> <bedGraph> [<bedGraph>]
"
if [ $# -lt 2 ];then echo "$usage";return; fi
	local tmpd=`mymktempd`;
	if [ -f $1 ];then
		chromsize=$1;
	else
		ucsc.chrom_size $1 > $tmpd/a
		chromsize=$tmpd/a	
	fi
	local cmd='use strict; my %chrom=();
		open(F,"<","'$chromsize'") or die "cannot open $!";
		while(<F>){ chomp; my ($c,$l)=split /\t/,$_; $chrom{$c}=$l; }
		close(F);
		while(<STDIN>){ chomp; my @a=split/\t/,$_;
			my $l= (defined $chrom{$a[0]}) ? $chrom{$a[0]}: "";
			next if ($l eq "");	
			if($a[1] < 0){ $a[1]=0;}
			if($a[2] < 1){ $a[2]=2;}
			if($a[1] > $l-1){ $a[1]=$l-1;}
			if($a[2] > $l){ $a[2]=$l;}
			print join("\t",@a),"\n";
		}	
	'
		
	for f in ${@:2};do
		o=${f/.bedGraph/.bw};
		bed.sort $f | perl -e "$cmd" > $tmpd/b
		bedGraphToBigWig $tmpd/b $tmpd/a $o
	done
	rm -rf $tmpd
}

ucsc.bg2bw.test(){
echo \
"chr1	1	2	1
chr1	2	3	4
chr111	2	3	4
chr1	-2	1	4
chr1	249250621	249250629	100
chr2	3	4	7" > tmp.bedGraph
ucsc.bg2bw hg19 tmp.bedGraph 
rm -rf tmp.bedGraph tmp.bw
}

ucsc.url(){
	$HMHOME/bin/ucsc_url.pl $@
}
ucsc.bed12(){
	mycat $1 | perl -ne '
	    chomp;
	    my @aa = split /\t/,$_;
	    my ($bin,$name,$chr,$strand,$start,$end,$thickStart,$thickEnd,$blockCount,$blockStarts,$blockEnds,$id,$name2) = split /\t/, $_;
	    my $itemRgb = "255,0,0";
	    my $score = 0;
	    if(defined $name2){
	        $name = $name."|".$name2;
	    }
	    print $chr,"\t",$start,"\t",$end,"\t",$name,"\t",$score,"\t",$strand,"\t",$thickStart,"\t",$thickEnd,"\t",$itemRgb,"\t",$blockCount,"\t";
	    my @ss = split /,/,$blockStarts;
	    my @ee = split /,/,$blockEnds;
	    for(my $i=0;$i<$blockCount;$i++){
	        my $length = $ee[$i]-$ss[$i];
	        print $length,",";
	    }
	    print "\t";
	    for(my $i=0;$i<$blockCount;$i++){
	        my $relstart = $ss[$i]-$start;
	        print $relstart,",";
	    }
	    print "\n";
	'
}

###
#  Repeat handling in ucsc
#  ref: https://www.biostars.org/p/12705/


ucsc.bed_nestedrep(){
#0  `bin` smallint(6) NOT NULL default '0',
#  `chrom` varchar(255) NOT NULL default '',
#  `chromStart` int(10) unsigned NOT NULL default '0',
#  `chromEnd` int(10) unsigned NOT NULL default '0',
#  `name` varchar(255) NOT NULL default '',
#  `score` int(10) unsigned NOT NULL default '0',
#  `strand` char(1) NOT NULL default '',
#7  `thickStart` int(10) unsigned NOT NULL default '0',
#8  `thickEnd` int(10) unsigned NOT NULL default '0',
#  `reserved` int(10) unsigned NOT NULL default '0',
#10  `blockCount` int(11) NOT NULL default '0',
#11  `blockSizes` longblob NOT NULL,
#12  `chromStarts` longblob NOT NULL,
#13  `blockStrands` longblob NOT NULL,
#  `id` int(10) unsigned NOT NULL default '0',
#  `repClass` varchar(255) NOT NULL default '',
#  `repFamily` varchar(255) NOT NULL default '',
	mycat $1 | perl -ne 'chomp; my @a=split/\t/,$_;
		my $chr=$a[1];
		my $start=$a[2];
		my $end=$a[3];
		my $name=$a[4];
		my @sizes=split/,/,$a[11];
		my @starts=split/,/,$a[12];
		my @strands=split/,/,$a[13];
		my ($id,$cl,$fm) = @a[14..16];
		for(my $i=0; $i<$a[10];$i++){
			print $chr,"\t",$start + $starts[$i];
			print "\t",$start + $starts[$i] + $sizes[$i];
			print "\t$name|$cl|$fm\t0\t",$strands[$i],"\n";
		}
	'	
}
ucsc.bed_nestedrep.test(){
echo "585	chr1	23803	24448	L2b	49	+	23803	24448	0	2	235,194,	0,451,	+,+,	12	LINE	L2" \
| ucsc.bed_nestedrep -
}

