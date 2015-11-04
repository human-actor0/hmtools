. $HMHOME/src/root.sh 
. $HMHOME/src/sam.sh

point(){
usage="
FUNCTION: extract 1 bp ustream position of 3'end and it's nucleotide from the reads 
USAGE: $FUNCNAME <bed> <relpos> <size>
"
if [ $# -ne 3 ];then echo "$usage" >&2; return; fi
	cat $1 | perl -e 'use strict; my %res=();
		my $relpos='$2'; my $size='$3';
		while(<STDIN>){ chomp; my @a=split /\t/,$_; 
			## cacluate 1bp upstream of the 3 prime end of reads
			## read           <_*__________/
			## mRNA ------------------------->>
			my $chrom=$a[0];
			my $len=$a[2]-$a[1];
			my $seq=$a[3];
			my $strand=$a[5];
			my $start=-1; 
			my $end=-1; 
			if($strand eq "+"){
				$start = $relpos < 0 ? $len + $relpos : $relpos;
				$end = $start + $size;
			}else{
				$end= $relpos < 0 ? -$relpos : $len - $relpos ;
				$start = $end - $size; 

			}
			my $subseq=substr($seq,$start, $end - $start);
			if($strand eq "-"){
				$subseq=reverse($subseq); $subseq=~tr/ACGTacgt/TGCAtgca/;
			}

			$start += $a[1]; $end += $a[1];
			my $id=join("\t",($chrom,$start,$end,$strand));
			$res{ $id }{$subseq} ++;
		}
		foreach my $k (keys %res){
			my @b=split/\t/,$k;
			my @name=();
			my $s=0; foreach my $k2 (keys %{$res{$k}}){ 
				my $v = $res{$k}{$k2};
				push @name,$k2.":".$v;
				$s += $v;
			}
			print join("\t",@b[0..2]),"\t",join(",",@name),"\t$s\t",$b[3],"\n";
		}
	' 
}
test__point(){
echo \
"chr1	100	110	ACGTgacgat	0	+
chr1	100	110	ACGTgacgat	0	-" > inp
echo \
"chr1	9	10	t:1	1	+
chr1	0	1	T:1	1	-
chr1	0	2	GT:1	1	-
chr1	8	10	at:1	1	+
chr1	7	10	atc:1	1	-
chr1	0	3	ACG:1	1	+" > exp
point inp -1 1 > obs
point inp -2 2 >> obs
point inp 0 3 >> obs
check obs exp
rm -rf obs exp inp
}

