. $HMHOME/src/root.sh;

#eval type mymktempd

fastq.flat(){
	awk '{printf("%s",$0); if(NR%4==0){ printf("\n");}else{ printf("\t");}}' $1;
}

fastq.cut5p_into_head(){
usage="
FUNCTION: Cut 5p barcode and add it to the header (at the end of first non-space tocken)
USAGE: $FUNCNAME <fastq> <len>
"
if [ $# -ne 2];then echo "$usage";return; fi
	fastq.flat $1 \
        | perl -ne 'chomp; my @a=split/\t/,$_; my $L="'$2'";
                my @b=split /\s+/,$a[0];
                $a[0] = $b[0].":".substr($a[1],0,$L);
                $a[1] = substr($a[1],$L);
                $a[3] = substr($a[1],$L);
                print join("\n",@a),"\n";
        '

}

fq.cut3p_go(){
usage="$FUNCNAME <input> <3adaptor> [<num_mis_match>]";
local cmd='
package main
import (
    "bufio"
    "fmt"
	//"strings"
	"log"
    "os"
)
func main() {
	var adapt string="'$2'";
	var bufs []string 
	//file, err := os.Open(os.Stdin)
	//if err != nil { log.Fatal(err) }
	//defer file.Close()
	var misMatch int = '${3:-0}'; //fmt.Println(misMatch);
	var minLen int=0; _ = minLen;
	var minMatch int=len(adapt) - misMatch; _=minMatch;
	var line_num int = 0;
	var i,j,m int; _,_,_ = i,j,m;

    	//scanner := bufio.NewScanner(file)
    	scanner := bufio.NewScanner(os.Stdin)
    for scanner.Scan() {
	var line string=scanner.Text();
	bufs = append(bufs,line);
	if len(bufs) == 4 {
		for i := len(bufs[1])-1; i >= 0; i--{
			m = 0;
			for j :=0; j< len(adapt) && i+j < len(bufs[1]); j++{
				if bufs[1][i+j] == adapt[j] { m ++;}
				//fmt.Printf("comp %s %s %d\n",bufs[1][(i+j):], adapt[j:],m);
			}
			if m >= minMatch {
				fmt.Printf("%s\n",bufs[0])
				fmt.Printf("%s\n",bufs[1][:i])
				fmt.Printf("%s\n",bufs[2])
				fmt.Printf("%s\n",bufs[3][:i])
				//fmt.Printf("%s\n",strings.Join(bufs,"\n"))
				break; }
		}
		bufs= nil;
	}

	line_num ++;
    }

    if err := scanner.Err(); err != nil {
        log.Fatal(err)
    }
}
'
if [ $# -lt 2 ];then
	echo "$usage"; return;
fi
local tmp=`mymktempd`; echo "$cmd" > $tmp/script.go;
cat $1 | go run $tmp/script.go

}
fq.cut3p_go.test(){
echo \
"@a
CGGGCTTGAACACG
+
IIIIIIIIIIIIII
@b
NGGAANAAAGAAAGGCAGACTGCCAGAACCAGCGCCTCATTTAATCTCGTATGCCGTCTTCTGCTTGAAAAAGAAT
+
#AA<A#FFFFFFFFFFFF)FFFFF.FFFFF<AFFFFFFAF<FF.<AAF<F.AAAFAFFAAF<AFFF7.FFFF7A..
@c
NAGAACAAGCCTACAGCACCCGGTATCTAGTATGCCGTCTTCTGCTTGAAAAAAAAAAGGGGGGGGGGGGGGGGGG
+
#AAAA#FFFFFFFFFFF7FFFFFFFAF<)F<<FF)FFF7.FF.7F)FAFFFFFFFFF77FFFFFFFFFFFFFFFFF
@d
CGGGCTTGAACACG
+
IIIIIIIIIIIIII"  > tmp.a

fq.cut3p_go tmp.a GAAC  
#fq.cut3 tmp.a GAAC 
#fq.cut3 tmp.a GAACACG 
#fq.cut3 tmp.a AGAACACG 1 
#fq.cut3 tmp.a GAACACGT 0 
#fq.cut3 tmp.a GGGGGGG 0 
#rm tmp.a
}



fastq.cut5(){ usage="
$FUNCNAME <fastq> <len> 
"
	if [ $# -ne 2 ]; then echo "$usage"; return; fi
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
	mycat $1 | fq.flat - \
	| perl -e 'use strict; 
	my $a="'$2'"; my $S=1; my $M='${3:-0}';
	my $T=0; my $T1=0; 
	while(<STDIN>){chomp; my ($id,$s,$o,$q)=split/\t/,$_;
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
	}
	print {*STDERR} $T1,"/",$T," passed\n";
	'

}



