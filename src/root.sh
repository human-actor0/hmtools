#!/bin/bash

mycat(){
        if [[ ${1##*.} = "gz" ]];then
                gunzip -dc $1;
        elif [[ ${1##*.} = "bz2" ]];then
                bunzip2 -dc $1;
        elif [[ ${1##*.} = "zip" ]];then
                unzip -p $1;
        else
                cat $1;
        fi
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

make_tempdir(){
	mktemp -d 
}

run_R(){
usage="
usage: $FUNCNAME <Rscript> [<output_dir>]
"
	if [ $# -gt 1 ];then
		mkdir -p $2; tmpd=$2;
	elif [ $# -gt 0 ]; then
		tmpd=`make_tempdir`
	else
		echo "$usage"; return;
	fi
        cmd=$1;
        tmpout=$tmpd/out;
        cmd=${cmd/stdout/$tmpout}
        echo "$cmd" > $tmpd/cmd
        R --no-save -q -f $tmpd/cmd &> $tmpd/log;
	if [ -f $tmpout ];then
        	cat $tmpout
	fi
	if [ $# -lt 2 ];then
		rm -rf $tmpd
	fi
}


###########################################################
# test 
###########################################################
test(){
	echo "too clean to test"
}


