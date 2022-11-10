#!/bin/sh
# calculate the predicted arrivals
# required tools: taup

# Usage
# ./caltime-taup.sh 202*

cmp=shz
#the component, e.g. z r t 
for eve
do
   cd $eve
   rm phase.arr
   evlat=$(saclst evla f `ls *.$cmp | head -1` | awk '{print $2}')
   evlon=$(saclst evlo f `ls *.$cmp | head -1` | awk '{print $2}')
   dep=$(saclst evdp f `ls *.$cmp | head -1` | awk '{print $2}')
   dmin=`saclst dist f *.$cmp | sort -g -k 2 | head -1 | awk '{print $2-200}'`
   dmax=`saclst dist f *.$cmp | sort -g -k 2 | tail -1 | awk '{print $2+200}'`
   saclst kstnm stla stlo dist f *.$cmp | awk '{print $2,$3,$4,$5}' > junk1
   for sta in $(cat junk1 | awk '{print $1}')
   do 
      stla=$(grep $sta junk1 | awk '{print $2}')
      stlo=$(grep $sta junk1 | awk '{print $3}')
      dist=$(grep $sta junk1 | awk '{print $4}')
echo $stla $stlo $dist

############  calculate the P and S arrivals #######################################
      P=`taup_time -h $dep -sta $stla $stlo -evt $evlat $evlon -ph P  --time | head -1 | awk '{print $1}'`
      if test -z "$P"
      then
      P=`taup_time -h $dep -sta $stla $stlo -evt $evlat $evlon -ph p --time | head -1 | awk '{print $1}'`
      else
      echo "P exists!"
      fi
      echo $sta $stla $stlo $evlat $evlon  $dist $P $S >> phase.arr
      echo $sta $stla $stlo $evlat $evlon  $dist $P $S 
   done
   rm junk*
  cd ..
done
