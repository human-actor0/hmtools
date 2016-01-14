#!/bin/bash

. $HMHOME/src/root.sh

ucsc.bed12(){
	mycat $1 perl -ne '
	    chomp;
	    my @aa = split /\t/,$_;
	    my ($bin,$name,$chr,$strand,$start,$end,$thickStart,$thickEnd,$blockCount,$blockStarts,$blockEnds,$id,$name2) = split /\t/, $_;
	    my $itemRgb = "255,0,0";
	    my $score = 0;
		print $blackCount,"\n"; exit;

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

