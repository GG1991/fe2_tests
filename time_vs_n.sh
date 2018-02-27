#!/bin/bash

program=run_periodic_ms.m
line_nx=18
line_ny=19

#solver=("cg_pd" "lu" "cg")
solver=("cg" "cg_pd" "cg_pgs")
for i in ${solver[@]}; do if [ -f "$i.dat" ]; then rm "$i.dat"; fi; done

nvals=(10 20 30 40 50 60 70)
#nvals=(10 20 30 40 50 60 70 80 90 100)

for s in ${!solver[@]}; do
 echo "${solver[$s]}"
 for n in ${!nvals[@]}; do
 
  printf "running for n = %3d  " ${nvals[$n]}
  sed -e "${line_nx}s/.*/global nx = ${nvals[$n]};/" -e "${line_ny}s/.*/global ny = ${nvals[$n]};/" $program > program_proto.m
  # default is cg_pd
  if   [ ${solver[$s]} == "lu" ]; then
    sed -i '/b;/ s/%/ /'       program_proto.m
    sed -i '/cg_pd/  s/^\s/%/' program_proto.m
  elif [ ${solver[$s]} == "cg" ]; then
    sed -i '/cg(/ s/%/ /'      program_proto.m
    sed -i '/cg_pd/  s/^\s/%/' program_proto.m
  elif [ ${solver[$s]} == "cg_pgs" ]; then
    sed -i '/cg_pgs/ s/%/ /'   program_proto.m
    sed -i '/cg_pd/  s/^\s/%/' program_proto.m
  fi
  sed -i '49s/.*/for i = 1 : 1/' program_proto.m
  octave program_proto.m > out_proto.dat
 
  time[$n]=$(awk '/time/{print $12;exit;}' out_proto.dat)
  its[$n]=$(awk '/time/{print $9;exit;}' out_proto.dat)
  printf "time = %f its = %3d\n" ${time[$n]} ${its[$n]}

  printf "%3d %f %3d\n" ${nvals[$n]} ${time[$n]} ${its[$n]} >> "${solver[$s]}.dat"
 
 done
done

echo ${time[@]}


