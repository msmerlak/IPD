#!/bin/bash

#SBATCH --ntasks=10
#SBATCH --partition=epyc
#SBATCH --mail-user=smerlak@mis.mpg.de
#SBATCH--mail-type=ALL

julia-1.6 /usr/people/smerlak/projects/iterated-prisoners-dilemma/IPD/_research/test.jl