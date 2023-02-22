#!/bin/bash
#PBS -q small
##PBS -m abe
##PBS -M (メールアドレス)
#PBS -l nodes=1:ppn=4
#PBS -l walltime=00:30:00

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=4

logfile=$(date +log.%Y-%m-%d-%H%M%S)

./a.out >& $logfile
