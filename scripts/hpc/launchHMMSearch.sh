#!/bin/bash

cd /mnt/research/germs/shane/transActRNA/data

batchn=1
for fname in /mnt/research/germs/shane/transActRNA/scripts/hpc/hmmSearch/*.sb;
do
    ((batchn++))
    #baseFileName=$(basename "$fname" .sb)
    echo $batchn $fname
    sbatch --job-name=HMM$batchn -A shadeash-colej --cpus-per-task=1 --mem=2G --output=/mnt/research/germs/shane/transActRNA/data/logs/HMM_5_02_run_$batchn.log --ntasks=1 --time=1:00:00 $fname
done