#!/bin/bash

dig_sub(){
	for i in `seq 1 255`;do
		echo $1.$i `dig +short -x $1.$i`
	done
}
