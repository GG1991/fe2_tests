#!/bin/bash

program=run_periodic_ms.m
line_nx=18
line_ny=19

#nvals=(5 10 15 20)
nvals=(10 20 30 40 50 60 70 80 90 100)

for n in ${!nvals[@]}; do

 printf "running for n = %3d  " ${nvals[$n]}
 sed -e "${line_nx}s/.*/global nx = ${nvals[$n]};/" -e "${line_ny}s/.*/global ny = ${nvals[$n]};/" $program > program_proto.m
 octave program_proto.m > out_proto.dat

 time[$n]=$(awk '/time/{print $12;exit;}' out_proto.dat)
 printf "time = %f\n" ${time[$n]}

done

echo ${time[@]}


