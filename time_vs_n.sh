#!/bin/bash

nvals=(10 20 30)
#nvals=(10 20 30 40 50 60 70 80 90 100)

solver=("lu" "cg" "cg_pd" "cg_pgs")
#solver=("cg" "cg_pd" "cg_pgs")
for i in ${solver[@]}; do if [ -f "$i.dat" ]; then rm "$i.dat"; fi; done

for s in ${!solver[@]}; do
 echo "${solver[$s]}"
 for n in ${!nvals[@]}; do
 
  printf "running for n = %3d  " ${nvals[$n]}
  echo "octave run_lin.m "-$s" -ustrain -nx ${nvals[$n]} -ny ${nvals[$n]} > output.dat"
 
  time[$n]=$(awk '/time/{print $12;exit;}' output.dat)
  its[$n]=$(awk '/time/{print $9;exit;}' output.dat)
  printf "time = %f its = %3d\n" ${time[$n]} ${its[$n]}

  printf "%3d %f %3d\n" $((nvals[$n] * nvals[$n])) ${time[$n]} ${its[$n]} >> "ustrain_${solver[$s]}.dat"
 
 done
done

echo ${time[@]}
