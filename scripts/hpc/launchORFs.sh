#!/bin/bash

cd /mnt/research/germs/shane/transActRNA/data

batchn=1
for fname in /mnt/research/germs/shane/transActRNA/scripts/hpc/orfs/*.sb;
do
    ((batchn++))
    #baseFileName=$(basename "$fname" .sb)
    echo $batchn $fname
    sbatch --job-name=ORFs$batchn -A shadeash-colej --cpus-per-task=1 --mem=10G --output=/mnt/research/germs/shane/transActRNA/data/logs/ORFs_5_3_run_$batchn.log --ntasks=1 --time=0:45:00 $fname
done