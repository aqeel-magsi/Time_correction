#!/bin/sh

#extract the event information for all events including lon,lat,day and date 

path=/home/aqeel/Aqeel_practice/OBS_Correction/Waveforms_Events/Tele_Events
sta=H36
cmpn=shz
FILE1=time_picks.txt
FILE2=phase.arr
output1=$sta'_pred.txt'
output2=$sta'_pick.txt'
output3='local_'$cmpn'_'$sta'.dat'
#rm $path/$output1 $path/$output2 $path/$output3
	cd $path || exit
	rm -f $FILE
	rm $output1 $output2
	rm /home/aqeel/Aqeel_practice/OBS_Correction/Waveforms_Events/Local_Events/diff_results/$output3 
	for dir in *_2018 *_2019; 
	do
    		#cd $dir || continue 
		cd $dir
    			if [ -f $FILE1 ] && [ -f $FILE2 ];
       			then
			echo $dir "required file exist"
			day=$(echo $dir | awk '{print substr($1,1,3)}')
			echo $day "day"
			data=$(echo $day $sta $cmpn | awk '{print $1 "_" $2 "." $3}')
			echo $data "data"
			evlat=$(saclst evla f $data | awk '{print $2}')
   			evlon=$(saclst evlo f $data | awk '{print $2}')
			date=$(saclst kzdate f $data | awk '{print $2}')
			echo  $date $evlon $evlat $day >> $path/H36_tele.dat
			#paste $path/$output1 $path/$output2 |awk '{print substr($4,1,7),$5,$7, $6-$3}' >> $path/$output3
                	#cat tempFile >> ../34_pred.txt
                	#> tempFile
			else
			echo $dir "file not exist"
			fi
    		cd ..
	done
#rm $output1 $output2
#mv $path/$output3 /home/aqeel/Aqeel_practice/OBS_Correction/Waveforms_Events/Local_Events/diff_results

#289_2018 335_2018 350_2018 354_2018 363_2018 006_2019 043_2019 048_2019 083_2019 089_2019 113_2019 108_2019 134_2019 151_2019 188_2019 179_2019 192_2019 195_201
