. $HMHOME/src/root.sh

mhits.clique(){
	cat $1 | perl -e 'use strict;
		my %res=();
		while(<STDIN>){ chomp; 
			my ($t,$r) = split/\t/,$_;
			$res{$r}{$t}=1;
		}
		foreach my @$r (keys %res){
			if($#r > 0){
			}
		}
	';
}

mhits.clique.test(){
echo "t1	r1
t1	r2
t2	r1
t2	r2
t3	r1" | mhits.clique -
}
