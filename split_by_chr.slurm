#!/bin/bash
# module load plink2


cp ${WORKING_DIR}/sumstats/${infile} $TMPDIR/.


echo "working on $study"

for chr in {1..23} X
do
  zcat "$TMPDIR"/$infile | awk -v chr=$chr -v chrcol=$chrcol '{if (NR == 1 || $chrcol == chr ) print}'  | gzip > /home/pgca1pts/freeze3_gwas/sumstats/bychr/${infile%.txt.gz}.chr${chr}.txt.gz
done
