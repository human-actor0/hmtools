. $HMHOME/src/bed.sh

bruseq.bedgraph(){
usage="
FUNCT:  count the number of intervals overlapping with fixed bins
USAGE: $FUNCNAME <bed> <bin_size> <out_name>
	will produce <outdir>/fwd.bedGraph and <outdir>/bwd.bedGraph
"
if [ $# -lt 3 ]; then echo "$usage"; return; fi
	local tmpd=`mymktempd`;
	mycat $1 | awk -v OFS="\t" -v B=$2 -v O=${tmpd} '{
		fout=O"/"$1","$6;
		s=int($2/B);e=int(($3-1)/B);	
		for(i=s; i<=e; i++){
			print i >> fout;
		}
	}'	

	rm -rf ${3}_fwd.bedGraph ${3}_bwd.bedGraph
	for f in $tmpd/*;do
		n=( `echo ${f##*/} | tr "," " "` );
		sort -n $f | uniq -c \
		| awk -v C=${n[0]} -v S=${n[1]} -v OFS="\t" -v B=$2 -v O=$3 '{ 
			if(S == "+"){ fout=O"_fwd.bedGraph";
			}else{ fout=O"_bwd.bedGraph"; }
			print C,$2*B,($2+1)*B,$1 >> fout;
		}' 
	done
	rm -rf $tmpd;
}
bruseq.bedgraph.test(){
echo \
"chr1	100	200	.	255	+
chr1	110	220	.	255	+
chr1	100	200	.	255	-
chr1	110	220	.	255	-
chr2	100	200	.	255	-
chr2	110	220	.	255	-" | bruseq.bedgraph - 50 obs
cat obs_* > obs 
echo \
"chr1	100	150	2
chr1	150	200	2
chr1	200	250	1
chr2	100	150	2
chr2	150	200	2
chr2	200	250	1
chr1	100	150	2
chr1	150	200	2
chr1	200	250	1" > exp
check obs exp
rm -rf obs exp

}

bruseq.bedgraph_smart(){
usage="
FUNCT: bed6 to 1bp resolution bedGraph
USAGE: $FUNCNAME <bed> <out_name>
	will produce <outdir>/fwd.bedGraph and <outdir>/bwd.bedGraph
EXAMP:
	INPUT
	[  1  )   [  2  )
	  [3 )  [  4  )

	OUTPUT
	[1) [1) 
          [ )
"
if [ $# -lt 2 ]; then echo "$usage"; return; fi
	rm -rf ${2}_fwd.bedGraph ${2}_bwd.bedGraph;
	local tmpd=`mymktempd`;
	mycat $1 > $tmpd/a
	bed.flat $tmpd/a > $tmpd/b
	bed.count $tmpd/b $tmpd/a -s \
	| awk -v OFS="\t" -v O=$2 '{
		if($6=="+"){ fout=O"_fwd.bedGraph";	
		}else{ fout=O"_bwd.bedGraph"; }
		print $1,$2,$3,$7 >> fout;
	}'
	rm -rf $tmpd;
}

bruseq.bedgraph_smart.test(){
echo \
"chr1	100	200	.	1	+
chr1	150	250	.	2	+
chr1	100	200	.	1	-
chr2	110	220	.	1	-" | bruseq.bedgraph_smart - tmp
head tmp_fwd.bedGraph tmp_bwd.bedGraph > tmp.obs
echo \
"==> tmp_fwd.bedGraph <==
chr1	100	150	.	1	+
chr1	200	250	.	2	+
chr1	150	200	.	3	+

==> tmp_bwd.bedGraph <==
chr1	100	200	.	1	-
chr2	110	220	.	1	-" > tmp.exp
check tmp.obs tmp.exp

rm tmp_*.bedGraph tmp.obs tmp.exp
}

bruseq.track(){
usage="$FUNCNAME <bed12> <track_name> <color>
	<color>: RGB (ex. black = 0,0,0 )
"
if [ $# -ne 3 ];then echo "$usage"; return; fi
## separate by strand and splicing or unsplicing	
	local tmpd=`mymktempd`;
	local track_name=$2;
	mycat $1 > $tmpd/a; 
	bed.exon $tmpd/a | bruseq.bedgraph_smart - $tmpd/b

	if [ -f $tmpd/b_fwd.bedGraph ];then
		echo "track name=${track_name}_fwd type=bedGraph color=0,0,0 parent=${track_name}"
		cat $tmpd/b_fwd.bedGraph 
	fi
	if [ -f $tmd/b_bwd.bedGraph ];then
		echo "track name=${track_name}_bwd type=bedGraph color=0,0,0 parent=${track_name}"
		cat $tmpd/b_bwd.bedGraph
	fi
	rm -rf $tmpd
}
bruseq.track.test(){
echo \
"chr1	1000	2000	g1	0	+	1000	2000	0,0,0	3	100,100,100	0,500,900"> tmp.gene
echo \
"chr1   1050    1550    r1      1       +       1050    1550    255,0,0 2       50,50   0,450
chr1   1050    1950    r1      1       +       1050    1950    255,0,0 2       50,50   0,850"> tmp.read
bruseq.track tmp.read track1 0,0,0
rm tmp.read tmp.gene
}


