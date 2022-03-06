#!/bin/sh

pathlist=$1 
linecount=$(wc -l $pathlist | awk '{print $1}')
echo $linecount

log=logs/slurm_%A.log 
sbatch -p hpc,vgl,vgl_bigmem -c 1 --error=$log --output=$log --array=1-$linecount gfastats_stats.sh $pathlist


