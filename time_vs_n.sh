#!/bin/bash

program=run_periodic_ms.m
line_nx=18
line_ny=19

nvals=(5 10 15 20)

sed -e "${line_nx}s/.*/global nx = ${nvals[0]};/" -e "${line_ny}s/.*/global ny = ${nvals[0]};/" $program > program_proto.m
octave program_proto.m > out_proto.dat

