#!/bin/bash

nvals=(10 20 30 40 50 60 70 80)
#nvals=(10 20 30 40 50 60 70 80 90 100)

#------------------------------------------------------------
# ustrain
solver=("lu" "cg" "cg_pd" "cg_pgs")

for i in ${solver[@]}; do if [ -f "ustrain_$i.dat" ]; then rm "ustrain_$i.dat"; fi; done

for s in ${!solver[@]}; do
 echo "${solver[$s]}"
 for n in ${!nvals[@]}; do
 
  printf "running for n = %3d  " ${nvals[$n]}
  octave run_lin.m "-${solver[$s]}" -ustrain -nx ${nvals[$n]} -ny ${nvals[$n]} -nexp 1 > output.dat
 
  time[$n]=$(awk '/time/{print $12;exit;}' output.dat)
  its[$n]=$(awk '/time/{print $9;exit;}' output.dat)
  printf "time = %f its = %3d\n" ${time[$n]} ${its[$n]}
  printf "%3d %f %3d\n" $((nvals[$n] * nvals[$n])) ${time[$n]} ${its[$n]} >> "ustrain_${solver[$s]}.dat"
 
 done
done

#------------------------------------------------------------
# ustress
solver=("lu")

for i in ${solver[@]}; do if [ -f "ustress_$i.dat" ]; then rm "ustress_$i.dat"; fi; done

for s in ${!solver[@]}; do
 echo "${solver[$s]}"
 for n in ${!nvals[@]}; do
 
  printf "running for n = %3d  " ${nvals[$n]}
  octave run_lin.m "-${solver[$s]}" -ustress -nx ${nvals[$n]} -ny ${nvals[$n]} -nexp 1 > output.dat
 
  time[$n]=$(awk '/time/{print $12;exit;}' output.dat)
  its[$n]=$(awk '/time/{print $9;exit;}' output.dat)
  printf "time = %f its = %3d\n" ${time[$n]} ${its[$n]}
  printf "%3d %f %3d\n" $((nvals[$n] * nvals[$n])) ${time[$n]} ${its[$n]} >> "ustress_${solver[$s]}.dat"

 done
done

#------------------------------------------------------------
# per_ms
solver=("lu" "cg" "cg_pd" "cg_pgs")

for i in ${solver[@]}; do if [ -f "per_ms_$i.dat" ]; then rm "per_ms_$i.dat"; fi; done

for s in ${!solver[@]}; do
 echo "${solver[$s]}"
 for n in ${!nvals[@]}; do
 
  printf "running for n = %3d  " ${nvals[$n]}
  octave run_lin.m "-${solver[$s]}" -per_ms -nx ${nvals[$n]} -ny ${nvals[$n]} -nexp 1 > output.dat
 
  time[$n]=$(awk '/time/{print $12;exit;}' output.dat)
  its[$n]=$(awk '/time/{print $9;exit;}' output.dat)
  printf "time = %f its = %3d\n" ${time[$n]} ${its[$n]}
  printf "%3d %f %3d\n" $((nvals[$n] * nvals[$n])) ${time[$n]} ${its[$n]} >> "per_ms_${solver[$s]}.dat"
 
 done
done

#------------------------------------------------------------
# per_lm
solver=("lu")

for i in ${solver[@]}; do if [ -f "per_lm_$i.dat" ]; then rm "per_lm_$i.dat"; fi; done

for s in ${!solver[@]}; do
 echo "${solver[$s]}"
 for n in ${!nvals[@]}; do
 
  printf "running for n = %3d  " ${nvals[$n]}
  octave run_lin.m "-${solver[$s]}" -per_lm -nx ${nvals[$n]} -ny ${nvals[$n]} -nexp 1 > output.dat
 
  time[$n]=$(awk '/time/{print $12;exit;}' output.dat)
  its[$n]=$(awk '/time/{print $9;exit;}' output.dat)
  printf "time = %f its = %3d\n" ${time[$n]} ${its[$n]}
  printf "%3d %f %3d\n" $((nvals[$n] * nvals[$n])) ${time[$n]} ${its[$n]} >> "per_lm_${solver[$s]}.dat"
 
 done
done
