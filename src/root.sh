#!/bin/bash

IFS=$' \n\t';


crop_prefix(){
# obtained from http://www.alecjacobson.com/weblog/?p=141
# and modified here
	if [ "$#" -eq 0 ] || [ "$1" = "-" ];then
		local input=`cat $1 | tr " \t" "\n" | awk 'NF'`;
	else
		local input=`echo $@ | tr " \t" "\n" | awk 'NF'`;
	fi
	local longest=$(echo "$input" | wc -L)
	local count=$(echo "$input" | wc -l)
	local prefix=1
	if [ $count -gt 1 ]; then
		for i in $(seq $longest); do
			[ $(echo "$input" | uniq -w$i | wc -l) -eq 1 ] || break
			prefix=$(($prefix+1))
		done
	fi
	echo "$input" | cut -c$prefix-
}

crop_prefix.test(){
	crop_prefix abc abd abcef   
	echo "abc abd abcef" | crop_prefix -
}
myjoin(){
usage="
usage: $FUNCNAME [options] <file> [<file>..] 
 [options]:
	-d <str> : delimitor (default @)
	-e <str> : replace empty input with <str> (default NA)
";

	local empty="NA";
	local sep="@";
	local OPTIND
	while getopts ":d:e:" opt; do
	case $opt in
		d) sep=$OPTARG;;
		e) empty=$OPTARG;;
		\?) echo "Invalid option -$OPTARG"; return;;
	esac
	done
	shift $(( $OPTIND -1 ));
	if [ $# -lt 1 ] || [ $1 = "-" ] ;then 
		echo "$usage"; return;
	fi
	local names=`crop_prefix $@ | tr "\n" " "`;
        echo "$@" | perl -e '
                use strict;
		my $sep="'$sep'";
                my @files=split /\s+/,<STDIN>;
                my @names=split /\s+/,"'"$names"'";
                my $i=0;
                my $n=scalar @files;
                my %data=();
                my %nc = ();
                foreach my $f (@files){
                        open my $in, "<", $f or die "$f not exists";
                        while(<$in>){ chomp; my @a=split/\t/,$_;
                                my $k= $a[0]; 
                                my $v= join($sep,@a[1..$#a]);
                                $data{$k}{$i}=$v;
                                if(!defined $nc{$i}){
                                        $nc{$i} = $#a;
                                }
                        }
                        close($in);
                        $i++;
                }
		print "id";
		for(my $i=0; $i < $n; $i++){
		for(my $j=0; $j < $nc{$i}; $j++){
			print "\t",$names[$i],".",$j;
		}}
		print "\n";
                foreach my $k (keys %data){
                        print $k;
                        for(my $i=0; $i < $n; $i++){
                                print "\t", defined $data{$k}{$i}? $data{$k}{$i} : join("\t",( "'$empty'" ) x $nc{$i}) ;
                        }
                        print "\n";
                }
        ' | tr "$sep" "\t"
	OPTIND=$optind
}

myjoin.test(){
echo \
"a	1	2
c	5	6
b	3	4" > tmp.a
echo \
"a	10	20	30
b	30	40	50" > tmp.b
echo \
"a	1
b	2" > tmp.c
	myjoin -d "," tmp.*
	rm tmp.*
}

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


