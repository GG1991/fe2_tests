#!/bin/bash
#SBATCH --job-name="job_plast"
#SBATCH --time=00-00:30:00
#SBATCH --output=job_plast.out
#SBATCH --error=job_plast.err
#SBATCH --nodes=1

octave run_plast.m -nx 20 -ny 20 -Sy_f 1.0e11 -Sy_m 1.0e17
