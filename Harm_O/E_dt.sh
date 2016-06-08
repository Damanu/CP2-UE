#!/bin/bash

dt=0.0001
x=1
v=0
job=10000
pt=1
c=0

while [ $c -lt 10 ] ; do
	let c=c+1
	echo $c
	echo $dt $x $v $job $pt | ./init_holf 
	./holf
done
