#!/bin/bash

IFS=$' \n\t';

## check executive files
bedGraphToBigWig(){
	$HMHOME/bin/bedGraphToBigWig $@
}


escape(){
cat $1 | perl -ne 'chomp;
	$_=~ s/([^a-zA-Z0-9])/\\$1/g;
	print $_,"\n";	
'
}

mycat(){
	if [[ -d $1 ]];then
		cat ${1%/}/*;
        elif [[ ${1##*.} = "gz" ]];then
                gunzip -dc $1;
        elif [[ ${1##*.} = "bz2" ]];then
                bunzip2 -dc $1;
        elif [[ ${1##*.} = "zip" ]];then
                unzip -p $1;
        else
                cat $1;
        fi
}

swapcol(){
	awk -v I=$2 -v J=$3 -v OFS="\t" '{
		if(J < 0){ J=NF+J+1;}
		a=$(I); $(I)=$(J);$(J)=a;	
		print $0;
	}' $1;
}
swapcol.test(){
echo \
"a	b	c
d	e	f" | swap_columns - 2 3
}

check(){ ## obtained from bamtools test
	if diff $1 $2; then
		echo ok
		return 1
	else
		echo fail
		return 0
	fi
}

rmempty(){
	for f in $@;do
		(( $(stat -c %s "$f") )) || rm "$f"
	done
}

mymktempd(){
	mktemp -d 2>/dev/null || mktemp -d -t 'hmtmpdir'
}

run_R(){
usage="
usage: $FUNCNAME <Rscript> [flag]
	[flag] := log: print log
"
        cmd=$1;
	local tmpd=`mymktempd`;
        cmd=${cmd//stdout/$tmpd/out}

        echo "$cmd" > $tmpd/cmd
        R --no-save -q -f $tmpd/cmd &> $tmpd/log;
	if [ -f $tmpd/out ];then
		cat $tmpd/out
	else
		cat $tmpd/log
	fi
	local flag=${2:-""};

	if [ "$flag" = "log" ];then 
		cat $tmpd/log
	fi
	rm -rf $tmpd;
}
unset -f quote
quote(){
	#perl -ne 'chomp; print join("'$del'", map{ "\"$_\"" } split( /\s+/,$_ )),"\n"';
	local a=( ${@/,/" "} ); a=( "${a[@]/#/\"}" ); a=( "${a[@]/%/\"}" ); echo "${a[*]}"
}

# array=( "${array[@]/%/_content}" )
# array=( "${array[@]/#/prefix_}" )

# array=( "${array[@]/%/_content}" )
# array=( "${array[@]/#/prefix_}" )


###########################################################
# test 
###########################################################
test(){
	echo "too clean to test"
}


