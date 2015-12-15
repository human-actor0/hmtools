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
