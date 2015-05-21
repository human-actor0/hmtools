#!/bin/bash
. $HMHOME/src/root.sh

test_lineartrend(){
cmd='
from scipy.stats.stats import pearsonr
from scipy.stats import chi2
import sys

def test_lineartrend(x1,y1,x2,y2):
	if len(x1) < 2 or len(x2) < 2 or sum(y1) == 0 or sum(y2) == 0:
		raise ValueError
	sx = [];
	sy = [];
	for i in range(len(x1)):
			for j in range(y1[i]):
				sx.append(x1[i]);
				sy.append(1); 
	for i in range(len(x2)):
			for j in range(y2[i]):
				sx.append(x2[i]);
				sy.append(2);
	if len(set(sx)) == 1 :
		return 0,1;
	r,pval = pearsonr(sx,sy);
	s = r*r*(len(sx)-1);
	pval = 1-chi2.cdf(s,1)
	return r,pval

def s2a(x):
	return [ int(float(x)) for x in x.split(",")];

for line in sys.stdin:
	if line[0] == "#":
		print line.rstrip();
		continue;
	a = line.rstrip().split("\t");
	x1,y1,x2,y2 = map(s2a, a[1:]);
	try:
		r,p = test_lineartrend(x1,y1,x2,y2)
		print "\t".join(map(str, a + [r,p]))
	except ValueError:
		1
'
	tmpd=`make_tempdir`;	
	echo "$cmd" > $tmpd/cmd
	cat $1 | python $tmpd/cmd

}

test(){
	echo "testing test_lineartrend .. "
	echo \
	"id1	1,2,3	10,20,30	1,2,3	300,200,100" > inp
	echo \
	"id1	1,2,3	10,20,30	1,2,3	300,200,100	-0.249029122546	1.62848179386e-10" > exp
	test_lineartrend inp > obs	
	check exp obs	
	rm inp exp obs
}
