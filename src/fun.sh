#!/bin/bash
## bash functional programming

## obtained from https://github.com/abesto/fun.sh/blob/master/fbash.sh



list() {
	for i in "$@"; do
		echo "$i"
	done
}

each_line(){
        fun="$@";
        awk 'NF' | while read -r; do
                list $REPLY | $fun;
        done
}


fold() {
	funct="$@"
	read acc
	while read elem; do
		acc="$(printf "$acc\n$elem" | $funct)"
	done
	echo $acc
}


lambda() {
    lam() {
        unset last
        for last; do
            shift
            if [[ $last = ":" ]]; then
                echo "$@"
                return
            else
                echo "read $last;"
            fi
        done
    }
    y="stdin"
    for i in "$@"; do
        if [[ $i = ":" ]]; then
            y="args"
        fi
    done
    if [[ "$y" = "stdin" ]]; then
        read funct
        eval $(lam "$@ : $funct")
    else
        eval $(lam "$@")
    fi
    unset y
    unset i
    unset funct
}

clean_line(){
	perl -ne 'chomp; $_=~s/^\s+//g; $_=~ s/#.+//g; $_=~ s/\s\s+/ /g; print $_,"\n";' | awk 'NF'
}
each_line(){
        fun="$@";
        clean_line | while read -r; do
                list $REPLY | $fun;
        done
}
extend_lines(){
usage="$FUNCNAME <lines1> <lines2>
output: all combinations of two lines
"
	oIFS=$IFS; IFS=$'\n';
	for line1 in `echo "$1" | clean_line`;do
		for line2 in `echo "$2" | clean_line`;do
			echo "$line1 $line2";
		done
	done
	IFS=$oIFS;
        #local cr="$2";
        #echo "$1" | each_line lambda a b : 'echo "$cr" | each_line lambda c : echo \$a \$b \$c\';
	#IFS=$OIFS;
}
test__extend_lines(){
a="
## do not print this
 a b c
d e f";
b="1 #not me
2";
extend_lines "$a" "$b"
}
#test__extend_lines
