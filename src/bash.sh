## here I collect some tips
local ifs="$IFS";
local a=a,b
IFS=,
echo $a; ## apply separator
echo "$a"; ## keep the original value
b=( $a ); ## define an array b
echo "${b[*]}" ## apply seperator 
echo ${b[*]} ## use default separator
IFS=":"; 
echo "${b[*]}" ## use default space separator
echo ${b[*]} ## use default array separator

IFS="$ifs"; 
