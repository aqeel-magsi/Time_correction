#!/bin/bash


i=1 
while read line;
do 
	echo $line; 

	name=$(echo $line |awk -F "[_]" '{print $1}') ; 

	cp $line $name"_$i"".mat";

	i=$(echo "$i+1" |bc); 

done < lst

