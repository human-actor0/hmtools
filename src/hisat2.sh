. $HMHOME/src/root.sh

hisat2.fixchr(){
	cat $1 | perl -ne 'chomp; 
	if($_=~/GL/){
	}else{
		$_=~s/(\d+|MT|X|Y)/chr$1/g; 
		$_=~s/MT/M/g; 
	}
	print $_,"\n";'
}

hisat2.fixchr.test(){
echo \
"1
Y
X
MT" | hisat2.fixchr -
}

