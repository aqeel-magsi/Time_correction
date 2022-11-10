#!/bin/sh
# calculate the predicted arrivals and plot event record profile
# required tools: GMT 4, pssac and taup

# Usage
# ./plt_eqprf.sh 201*
 
gmt set FRAME_WIDTH = 0.03
gmt set HEADER_FONT_SIZE 12p
cmp=hyd
#the component, e.g. z r t
amp=0.7
for eve
do
   cd $eve
   mapname=wfprf_$eve.ps
   gmt set HEADER_FONT_SIZE 16p
   gmt set HEADER_OFFSET 0
   dmin=`saclst dist f *.$cmp | sort -g -k 2 | head -1 | awk '{print $2-20}'`
   dmax=`saclst dist f *.$cmp | sort -g -k 2 | tail -1 | awk '{print $2+20}'`
   #T1=$(sort -g -k 7 phase.arr | head -1  | awk '{print $7}')
   #$T2=$(sort -g -k 7 phase.arr | tail -1 | awk '{print $7+20}')
   echo "time window" $T1 $T2

############ plot waveform using PSSAC  ##############################################
## -E option determines could refer to pssac usage, the pssac could use to plot waveform profile aligned by distance (k), azimuth (z), back-azimuth (b)  and traces (n) by changing the -E parameter.
# tn controls the align up time, n= -5(b), -3(o), -2(a)
   for com in $cmp
   do 
      ls *.$com
      ls *.$com |gmt pssac -JX17/10 -R450/510/$dmin/$dmax -C -Ekt-3 -M$amp -Fr -Ba10f5:"Time (s)":/a50f10:"Distance (km)"::."$eve/$cmp":WSne -K > $mapname

      nawk -v T1=450 '{printf "%8.2f %8.2f 12 0 0 1 %s\n",T1+1,$6,substr($1,1,4)}' phase.arr |gmt pstext -JX -R -O -K >> $mapname
      grep -e $cmp time_picks.txt | nawk '{print $3,$4}' |gmt psxy -JX -R -Sy0.7 -O -K -W5,blue  >> $mapname
      nawk '{print $7,$6}' phase.arr |gmt psxy -JX -R -Sy0.7 -O -K -W5,red  >> $mapname
      
      ps2pdf $mapname 
      rm $mapname
      cp *.pdf ..
   done
  cd ..
done
