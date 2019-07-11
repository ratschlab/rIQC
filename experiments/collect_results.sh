#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/m53

cat ${basedir}/bam_based/*.tsv | sort > ${basedir}/bam_based.tsv

for K in 0.1 0.3 0.5 10000 25000 50000
do
    cat ${basedir}/fastq_based/*_${K}_sample_a_ratio_uq.tsv > ${basedir}/fastq_based_${K}.tsv
done
