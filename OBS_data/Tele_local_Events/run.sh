#!/bin/bash

for file in *

do 
echo $file "fileeeeeeeeeeeeeeee"
cd $file 
sac << EOF
r *hyd *GUMO*
qdp off 
rmean
bp c 2 5
sort dist
ppk
saveimg ../$file.ps
q
EOF
cd ..
done
