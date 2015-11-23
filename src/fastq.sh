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
