#!/bin/bash

while getopts i:l:n:p: option
do
  case "${option}"
    in
      i) infile=${OPTARG};;
      l) infilepath=${OPTARG};;
      n) nodeuse=${OPTARG};;
      p) ploc=${OPTARG};;
    esac
done

#nodeuse is the number of simultaneous processes to launch. Decrease if memory errors occur

#Total number of job batches = Number of commands (here 22) / number of commands used per job, rounded up 
#Small error: this seems lead to processing chr 23 and d4 as well
 totjobs=$(( (22 + $nodeuse - 1 ) / $nodeuse ))

if [ ! -d "grm" ]; then
  mkdir grm
fi

#Make a GRM for all chromosomes
for job in $(seq 1 1 $totjobs)
do
 jstart=$((($job-1)*$nodeuse +1))
 jstop=$(($job*$nodeuse))
 for chr in  $(seq $jstart 1 $jstop)
 do
  chr="$PBS_ARRAYID"
  echo $PBS_ARRAYID
  #awk -v chr=$chr '{if ($1 != chr) print $2}' "$infilepath"/"$infile".bim > grm/"$infile"_nochr"$chr".snplist
  #--extract grm/"$infile"_nochr"$chr".snplist 
  $ploc --bfile "$infilepath"/"$infile" --autosome --make-rel square --out grm/"$infile"_nochr"$chr" --memory 2000 &
 done
 wait
done
wait

#Merge all grm





