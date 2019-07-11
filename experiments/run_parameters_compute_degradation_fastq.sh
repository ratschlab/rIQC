#!/bin/bash

set -e

for R in 10000 25000 50000 0.1 0.3 0.5
do
    bsub -J "degArr[1-199]" -R "rusage[mem=10000]" -W 4:00 -n 1 $(pwd)/compute_degradation_fastq_cluster.sh $R
done
