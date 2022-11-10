#!/bin/sh
#Cut data into 4 hours segment

export SAC_DISPLAY_COPYRIGHT=1
output_dir=~/Aqeel_practice/Noise_Correlation/data_hyd_6h
input_dir=/NAS2/Abbas/TS14_OBS
secc=21600	
for stnm in H36
do
	cd $stnm
	#echo $stnm
   	   for file in $input_dir/$stnm/20*/???.hyd
           do
		 echo $file
       		 #year=`saclst  kzdate f $file | awk '{print substr($2,1,4)}'`
        	 dat=`saclst  kzdate f $file| awk '{print substr($2,1,10)}'`
       		 jday=$(date -d $dat +%j)
        	 #jday=`saclst  kztime f $sacfile | awk '{print $3}'`
        	 hour=`saclst  kztime f $file | awk '{print substr($2,1,2)}'`
        	 min=`saclst  kztime f $file | awk '{print substr($2,4,2)}'`
       		 sec=`saclst  kztime f $file | awk '{print substr($2,7,9)}'`
                 name=$(echo $file | awk -F "/" '{print substr($7,1,3)}')
		 cmp=$(echo $file | awk -F "/" '{print substr($7,5,7)}')
		 year=$(echo $file | awk -F "/" '{print substr($6,1,4)}')
		 #name=$jday
		 echo $name  $year  "####################"
		 echo "jday" $name "hour" $hour "min" $min "sec" $sec $name $cmp "cmp"
#cd H34
#now we need to cut the data in 4 hours segments 
                i=1 
		start=0
	        while [ $i -lt 5 ]
		do 
		   name1=$name'_'$i'_'$cmp'.sac'
rm $name1
echo $name1  "####################"
sac << EOF
r $file
lh o
ch o gmt $year $name 00 00 00		
ch allt (0 - &1,o&) iztype IO
lh o
w $name1
cut o $start $secc
r $name1
w over
q
EOF
		  secc=$(echo "$secc+21600" |bc)
		  start=$(echo "$start+21600" |bc)
	          i=$(echo "$i+1" |bc)
name2=$name'_hyd.sac'
rm $output_dir/$stnm/$year/$name2
mv $name1 $output_dir/$stnm/$year
done
secc=21600
	done
done
