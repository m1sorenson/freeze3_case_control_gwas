#!/bin/bash
nodeuse=15

cd /home/pgca1pts/freeze3_gwas

counter=1
for file in $(awk 'BEGIN{FS="/"}{print $2}' $metafile)
do
  echo $file
  /home/maihofer/trauma_gwas/METAL-master/build/bin/metal metal_scripts/$file > errandout/"$file"_errorlogs 2>&1 &
  counter=$((counter+1))
  if [ $counter  -eq 8 ]
  then
    wait
    counter=1
  fi

done

wait
