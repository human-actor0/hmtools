. $HMHOME/src/root.sh;


fq.flat(){
	awk '{printf("%s\t",$0); if(NR%4==0){ printf("\n");}}' $1;
}
fq.filter5(){
usage="
USAGE: $FUNCNAME <fq> <barcode>
"
if [ $# -ne 2 ]; then echo "$usage"; return; fi
	mycat $1 | fq.flat -  \
	| perl -ne 'my $bc="'$2'";
		chomp; my @a=split/\t/,$_;
		if( $a[1]=~/^$bc/){
			print join("\n",@a),"\n";
		}
	'
}

fq.filter5.test(){
echo \
"@a
GATTACA
+
IIIIIII
@a
TTGATTA
+
IIIIIII" | fq.filter5 - TT
}
fq.cut3(){
if [ $# -lt 2 ];then echo "
$FUNCNAME <fastq> <3adapter> <mis>
";fi
	mycat $1 | fq.flat - | perl -ne 'BEGIN{ $T=0;$T1=0;} chomp; my $a="'$2'"; my $S=1; my $M='${3:-0}';
		my ($id,$s,$o,$q)=split/\t/,$_;
		#print $_," ",$a," ",length($s)," ",length($a),"\n";
		$T++;
		for(my $i=$S; $i<length($s)-length($a)+1+$M;$i++){
			my $m=0;
			for(my $j=0; $j < length($a); $j++){
				if( $i+$j > length($s)-1){ $m+= length($a)-$j; last}
				my $sj=substr($s,$i+$j,1);
				my $aj=substr($a,$j,1);
				if( $sj ne $aj){ $m++;}
			}
			if($m <= $M){
				$T1++;
				print $id,"\n";
				print substr($s,0,$i),"\n";
				print $o,"\n";
				print substr($q,0,$i),"\n";
				last;
			}
		}
		END{ print {*STDERR} print "$T1/$T passed\n";}
	'

}

fq.cut3.test(){
echo \
"@a
CGGGCTTGAACACG
+
IIIIIIIIIIIIII" > tmp.a
cat tmp.a

#fq.cut3 tmp.a GAAC 
#fq.cut3 tmp.a GAACACG 
#fq.cut3 tmp.a AGAACACG 1 
fq.cut3 tmp.a GAACACGT 0 
fq.cut3 tmp.a GAACACGT 1 
rm tmp.a
}


