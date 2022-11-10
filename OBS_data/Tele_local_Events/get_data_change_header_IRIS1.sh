#!/bin/bash
#this code will help us to get events information from file, downloaded from IRIS.
#help us to change the header based on event information

while read p; do
  echo "$p" > info.lst
#more info.lst
all=$(awk '{print $1}' info.lst)
location=$(awk '{print $6 $7 $8}' info.lst)
year=$(awk '{print substr($1,1,4)}' info.lst)
date=$(awk '{print substr($1,1,10)}' info.lst)
echo $date
day=$(date -d $date +%j)
hour=$(awk '{print substr($1,12,2)}' info.lst)
min=$(awk '{print substr($1,15,2)}' info.lst)
sec=$(awk '{print substr($1,18,2)}' info.lst)
evdp=$(awk '{print substr($4,1,6)}' info.lst)
evla=$(awk '{print $2}' info.lst)
evlo=$(awk '{print $3}' info.lst)
mag=$(awk '{print $5}' info.lst)
echo $location
echo $year $day $hour $min $sec 
echo $all
echo $evdp $evlo $evla $day "day"
#make directory for each event 
#''''''''''''''''''''''''''''''''''''''''
mkdir $day$hour$min"_"$year

#''''''''''''''''''''''''''''''''''''''''
#copy the data files for each events's day 
#get data for each day
dir=/NAS2/Abbas/TS14_OBS
for stnm in $dir/GUMO $dir/H?? $dir/K??
do
	sta=$(echo $stnm |awk -F "/" '{print $5}')
	echo $stnm "stttnm"
	echo $sta "staaaa"
        echo $year
   	for file in $stnm/$year
  	do
	echo $file
    	cp $file/$day.shz ./$day$hour$min"_"$year/$day"_"$sta.shz
	#cp $file/$day.shy ./$hour"_"$day"_"$year/$day"_"$sta.shy
    	#cp $file/$day.shx ./$hour"_"$day"_"$year/$day"_"$sta.shx
    	cp $file/$day.hyd ./$day$hour$min"_"$year/$day"_"$sta.hyd
  	done
# after extracting all info from events file, Lets change the header 
for sac in ./$day$hour$min"_"$year/$day"_"$sta.sh* ./$day$hour$min"_"$year/$day"_"$sta.hyd
 do
echo $sac
sac << EOF
r $sac
qdp off 
rmean
synchronize 
ch evdp $evdp evla $evla evlo $evlo
ch o gmt $year $day $hour $min $sec
ch allt (0 - &1,o&) iztype IO
ch user0 $mag
wh 
cut 0 800
r $sac
w $sac
q
EOF
done
done
done <event_lsst.txt




