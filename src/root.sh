#!/bin/bash


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
usage: $FUNCNAME <Rscript> 
"
        cmd=$1;
	local tmpd=`mktemp -d`;
        cmd=${cmd/stdout/$tmpd/out}
        echo "$cmd" > $tmpd/cmd
        R --no-save -q -f $tmpd/cmd &> $tmpd/log;
	if [ -f $tmpd/out ];then
		cat $tmpd/out
	else
		cat $tmpd/log
	fi
	rm -rf $tmpd;
}


# array=( "${array[@]/%/_content}" )
# array=( "${array[@]/#/prefix_}" )


###########################################################
# test 
###########################################################
test(){
	echo "too clean to test"
}


