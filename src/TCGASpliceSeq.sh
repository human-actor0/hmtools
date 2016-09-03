TCGASpliceSeq.parse_table(){
cat PSI_download_BLCA.txt | awk '$4==10'\
 | perl -ne 'chomp; my @a=split/\t/,$_;
my $sum=0; my $n=0;
foreach my $v (@a[10..$#a]){
    if($v=~/[\d|\.]+/){ $sum += $v; $n++;}
}
print $sum," ",$n," ",$sum/$n,"\n";'
}

