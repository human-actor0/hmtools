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

extend_lines(){
usage="$FUNCNAME <lines1> <lines2>
output: all combinations of two lines
"
        local cr="$2";
        echo "$1" | each_line lambda a b : 'echo "$cr" | each_line lambda c : echo \$a \$b \$c\';
}
test__extend_lines(){
a=`echo -e "a b\nc d"`;
b=`echo -e "1\n2"`
extend_lines "$a" "$b"
rm a b
}
