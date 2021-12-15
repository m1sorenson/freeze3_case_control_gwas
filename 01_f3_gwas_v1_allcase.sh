#!/bin/bash

WORKING_DIR=/home/pgca1pts/freeze3_gwas
cd $WORKING_DIR
mkdir -p ${WORKING_DIR}/metal_results
mkdir -p ${WORKING_DIR}/metal_scripts
mkdir -p ${WORKING_DIR}/results_cat
mkdir -p ${WORKING_DIR}/results_filtered
mkdir -p ${WORKING_DIR}/errandout

# RUNTYPE options -
# males   - all male only studies
# females - all female only studies
# broad   - all studies
# narrow  - all EHR narrow definition studies
RUNTYPE=broad

if [ $RUNTYPE == "males" ]; then
  sex=males
elif [ $RUNTYPE == "females" ]; then
  sex=females
elif [ $RUNTYPE == "broad" ]; then
  sex=all
elif [ $RUNTYPE == "narrow" ]; then
  sex=all
else
  echo "RUNTYPE not defined/valid, options are: males, females, broad, narrow"
  exit 1
fi

### 1)Study Level Analysis steps:

##Group 1 and 2 GWAS only: GWAS

cov=pcs

#Make sure time codes are correct
IFS=$'\n'
for line in $(tail -n+2 config/dosage_locations_f3.csv)
do
  study_1=$(echo $line | awk 'BEGIN{FS=","}  {print $1}')
  study_2=$(echo $line | awk 'BEGIN{FS=","}  {print $2}')

  ancgroup=$(echo $line | awk 'BEGIN{FS=","} {print $3}')
  timecode=$(echo $line | awk 'BEGIN{FS=","} {print $4}')
  exclude=$(echo $line | awk 'BEGIN{FS=","}  {print $5}')

  if [ $exclude != "1" ]
  then
    echo gwas for $study_1 $study_2 $timecode
    sbatch --time=$timecode --error errandout/"$study_1"_"$study_2"_"$cov"_"$ancgroup"_"$sex".e --output errandout/"$study_1"_"$study_2"_"$cov"_"$ancgroup"_"$sex".o  \
    --export=ALL,study="$study_1",cov="$cov",study_2="$study_2",ancgroup="$ancgroup",sex="$sex" run_trauma_gwas_v2_freeze3_allcase.slurm -D $working_dir
  fi

done


### 2) Summary stat splitting
#Group 3 and 4 and 5 only: Summary stat conversion from genome-wide to split by chromosome
#Caution: these data may need to be reformatted first for meta-analysis. Do this PRIOR to running this script

#To reformat by chromosome, supply a list .gz files. Each row should be a file name followed by the column number for the chromosome
mkdir -p ${WORKING_DIR}/sumstats/bychr
while IFS= read -r line; do
  study=$(echo $line | awk '{print $1}')
  include_male=$(echo $line | awk '{print $2}')
  include_female=$(echo $line | awk '{print $3}')
  include_broad=$(echo $line | awk '{print $4}')
  chrcol=$(echo $line | awk '{print $5}')
  # study != study to skip first line
  if [[ $RUNTYPE == "males" ]] && [[ $study != "study" ]] && [ $include_male -eq 1 ]; then
    if [[ $study == "UKBB" ]] || [[ $study == "VETS" ]]; then
      infile=${study}_males_eur_allcase.txt.gz
    else
      infile=${study}_males_eur.txt.gz
    fi
    sbatch --time=00:15:00 --error errandout/"$infile"_split_by_chr.e --output errandout/"$infile"_split_by_chr.o  \
    --export=ALL,chrcol=${chrcol},infile=${infile},study=${study},WORKING_DIR=${WORKING_DIR} split_by_chr.slurm -D $WORKING_DIR
  elif [[ $RUNTYPE == "females" ]] && [[ $study != "study" ]] && [ $include_female -eq 1 ]; then
    if [[ $study == "UKBB" ]] || [[ $study == "VETS" ]]; then
      infile=${study}_females_eur_allcase.txt.gz
    else
      infile=${study}_females_eur.txt.gz
    fi
    sbatch --time=00:15:00 --error errandout/"$infile"_split_by_chr.e --output errandout/"$infile"_split_by_chr.o  \
    --export=ALL,chrcol=${chrcol},infile=${infile},study=${study},WORKING_DIR=${WORKING_DIR} split_by_chr.slurm -D $WORKING_DIR
  elif [[ $RUNTYPE == "broad" ]] && [[ $study != "study" ]] && [ $include_broad -eq 1 ]; then
    if [[ $study == "UKBB" ]] || [[ $study == "VETS" ]]; then
      infile=${study}_broad_eur_allcase.txt.gz
    else
      infile=${study}_broad_eur.txt.gz
    fi
    sbatch --time=00:15:00 --error errandout/"$infile"_split_by_chr.e --output errandout/"$infile"_split_by_chr.o  \
    --export=ALL,chrcol=${chrcol},infile=${infile},study=${study},WORKING_DIR=${WORKING_DIR} split_by_chr.slurm -D $WORKING_DIR
  fi
done < config/sumstat_studies.tsv

### 3) Metal script template generation

#User: Make a meta-analysis script. Run the next snippet to create the male/female/broad template
#Once the male/female/broad template has been created from the base template, go through and
#Make sure to open the respective template (f3_gwas_broad/f3_gwas_males/f3_gwas_females)
#after this runs to comment out any studies that should be excluded. Then you can run
#the script to split the metal script by chromosome and queue the slurm jobs

#IMPORTANT: Weights will be automatically filled in the metal script from the values contained in config/study_weights.tsv
#as well as frequency column names if the column name is different for broad/male/female (for the ones like FRQ_U_###)
if [ $RUNTYPE == "broad" ]; then
  # Add no extension to genotyped results and add _broad extension to summary stat files
  sed s/{EXT1}//g metal_templates/f3_gwas_case_control.mi | \
    sed s/{EXT2}/broad/g > metal_templates/f3_gwas_broad_allcase.mi
  # Fill in weights
  while IFS= read -r line; do
    study_id=$(echo "$line" | awk 'BEGIN{FS="\t"}{print $2}')
    if [[ $study_id != "study_id" ]]; then
      weight=$(echo "$line" | awk 'BEGIN{FS="\t"}{print $3}')
      #since the template METAL file only has {study_id_FREQ} for the ones that
      #need it, it's okay that this will be empty for some of the studies
      freq_col=$(echo "$line" | awk 'BEGIN{FS="\t"}{print $6}')
      echo {${study_id}_WEIGHT} $weight $freq_col
      sed s/{${study_id}_WEIGHT}/${weight}/g metal_templates/f3_gwas_broad_allcase.mi \
        | sed s/{${study_id}_FREQ}/${freq_col}/g \
        > metal_templates/f3_gwas_broad_allcase_tmp.mi
      mv metal_templates/f3_gwas_broad_allcase_tmp.mi metal_templates/f3_gwas_broad_allcase.mi
    fi
  done < config/study_weights_allcase.tsv
elif [ $RUNTYPE == "males" ]; then
  # Add _males extension to genotyped results and to summary stat files
  sed s/{EXT1}/_males/g metal_templates/f3_gwas_case_control.mi | \
    sed s/{EXT2}/males/g > metal_templates/f3_gwas_males_allcase.mi
  # Fill in weights
  while IFS= read -r line; do
    study_id=$(echo "$line" | awk 'BEGIN{FS="\t"}{print $2}')
    if [[ $study_id != "study_id" ]]; then
      weight=$(echo "$line" | awk 'BEGIN{FS="\t"}{print $4}')
      #since the template METAL file only has {study_id_FREQ} for the ones that
      #need it, it's okay that this will be empty for some of the studies
      freq_col=$(echo "$line" | awk 'BEGIN{FS="\t"}{print $7}')
      echo {${study_id}_WEIGHT} $weight $freq_col
      sed s/{${study_id}_WEIGHT}/${weight}/g metal_templates/f3_gwas_males_allcase.mi \
        | sed s/{${study_id}_FREQ}/${freq_col}/g \
        > metal_templates/f3_gwas_males_allcase_tmp.mi
      mv metal_templates/f3_gwas_males_allcase_tmp.mi metal_templates/f3_gwas_males_allcase.mi
    fi
  done < config/study_weights_allcase.tsv
elif [ $RUNTYPE == "females" ]; then
  # Add _females extension to genotyped results and to summary stat files
  sed s/{EXT1}/_females/g metal_templates/f3_gwas_case_control.mi | \
    sed s/{EXT2}/females/g > metal_templates/f3_gwas_females_allcase.mi
  # Fill in weights
  while IFS= read -r line; do
    study_id=$(echo "$line" | awk 'BEGIN{FS="\t"}{print $2}')
    if [[ $study_id != "study_id" ]]; then
      weight=$(echo "$line" | awk 'BEGIN{FS="\t"}{print $5}')
      #since the template METAL file only has {study_id_FREQ} for the ones that
      #need it, it's okay that this will be empty for some of the studies
      freq_col=$(echo "$line" | awk 'BEGIN{FS="\t"}{print $8}')
      echo {${study_id}_WEIGHT} $weight $freq_col
      sed s/{${study_id}_WEIGHT}/${weight}/g metal_templates/f3_gwas_females_allcase.mi \
        | sed s/{${study_id}_FREQ}/${freq_col}/g \
        > metal_templates/f3_gwas_females_allcase_tmp.mi
      mv metal_templates/f3_gwas_females_allcase_tmp.mi metal_templates/f3_gwas_females_allcase.mi
    fi
  done < config/study_weights_allcase.tsv
fi


### 4) Split METAL file by chromosome

#Before running next step, check the broad/males/females template and comment out
#any studies that should be excluded

#Adjust meta-analysis script for each chromosome (replaces "{CHR_NUM}" with a chromosome number)
for chr in {1..22} #X
do
  if [ $RUNTYPE == "broad" ]; then
    sed s/{CHR_NUM}/${chr}/g metal_templates/f3_gwas_broad_allcase.mi > metal_scripts/f3_gwas_broad_allcase_chr${chr}.mi
  elif [ $RUNTYPE == "males" ]; then
    sed s/{CHR_NUM}/${chr}/g metal_templates/f3_gwas_males_allcase.mi > metal_scripts/f3_gwas_broad_allcase_chr${chr}.mi
  elif [ $RUNTYPE == "females" ]; then
    sed s/{CHR_NUM}/${chr}/g metal_templates/f3_gwas_females_allcase.mi > metal_scripts/f3_gwas_broad_allcase_chr${chr}.mi
  fi
done

#List all chromosome meta-analysis files into metafilelist.txt file
ls  metal_scripts/f3_gwas_${RUNTYPE}_allcase_chr*.mi > metafilelist.txt

#User: Give a name for the error log files
dataset=PTSD_F3_v4_allcase_nov15_2021

sbatch -t 00:45:00  --error errandout/"$dataset".e --output errandout/"$dataset".o   --export=ALL,metafile=metafilelist.txt -D /home/pgca1pts/freeze3_gwas run_meta_v2_loo_v2.slurm


#Concatenate meta-analysis results
cat metal_results/eur_ptsd_pcs_v4_chr*_nov15_2021_1.tbl | awk '{if (NR == 1 || ($3 != "MarkerName" )) print}' > metal_results/eur_ptsd_pcs_v4_${RUNTYPE}_allcase_nov15_2021.tbl
#Format meta results for fuma #Note: COLUMNS WILL CHANGE if ANALYZE HET IS ON!
#User: find out what the max N is make sure only decently covered loci are returned (by default, 90% of total N)
#TBD: Don't do this for the x chromosome! N will definitely be less!
totalN=293992
grep -v ??? metal_results/eur_ptsd_pcs_v4_${RUNTYPE}_allcase_nov15_2021.tbl | LC_ALL=C sort -g -k 12 | gzip > results_filtered/eur_ptsd_pcs_v4_${RUNTYPE}_allcase_nov15_2021.fuma.gz

grep -v ??? metal_results/eur_ptsd_pcs_v4_${RUNTYPE}_allcase_nov15_2021.tbl | awk -v totalN=$totalN '{if(NR==1){$14 = "OBS_CT"}else{$14 = totalN} print}' > results_filtered/eur_ptsd_pcs_v4_${RUNTYPE}_allcase_nov15_2021.filtered.tbl
# create premunge.gz
LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/w_hm3.snplist.sorted \
  | LC_ALL=C sort -k1b,1 ) <(cat results_filtered/eur_ptsd_pcs_v4_${RUNTYPE}_allcase_nov15_2021.filtered.tbl \
  | awk '{if (NR == 1) $3="SNP"; print}' | LC_ALL=C sort -k3b,3 ) \
  | gzip > results_filtered/eur_ptsd_pcs_v4_${RUNTYPE}_allcase_nov15_2021.tbl.premunge.gz

# create munge.gz
python2 /home/maihofer/trauma_gwas/ldsc-master/munge_sumstats.py \
  --sumstats results_filtered/eur_ptsd_pcs_v4_${RUNTYPE}_allcase_nov15_2021.tbl.premunge.gz \
  --N-col OBS_CT --out results_filtered/eur_ptsd_pcs_v4_${RUNTYPE}_allcase_nov15_2021.tbl.munge.gz #add --N-col OBS_CT for the sample size

###LDSC analysis
# conda activate ldsc
#
# #Filter data down to just LDSC SNPs then munge
# LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(cat metal_results/eur_ptsd_pcs_v4_nov15_2021_${sex}.gz1.tbl | awk '{if (NR == 1) $3="SNP"; print}'   | LC_ALL=C sort -k3b,3 ) | gzip > results_filtered/eur_ptsd_pcs_v4_nov15_2021_${sex}.gz1.tbl.premunge.gz
# python2 /home/maihofer/trauma_gwas/ldsc-master/munge_sumstats.py --sumstats results_filtered/eur_ptsd_pcs_v4_nov15_2021_${sex}.gz1.tbl.premunge.gz  --N-col Weight --out results_filtered/eur_ptsd_pcs_v4_nov15_2021_${sex}.gz1.tbl.munge.gz #add --N-col OBS_CT for the sample size
#
# #Get format for LDhub
# zcat results_filtered/eur_ptsd_pcs_v4_nov15_2021_${sex}.fuma.gz | awk '{if (NR==1) print "CHR","BP","SNP","A1","A2","FRQ","N","Z","P"; else print $0}' > results_filtered/eur_ptsd_pcs_v4_nov15_2021_${sex}.gz1.tbl.munge.txt
# zip results_filtered/eur_ptsd_pcs_v4_nov15_2021_${sex}.gz1.tbl.munge.txt.zip results_filtered/eur_ptsd_pcs_v4_nov15_2021_${sex}.gz1.tbl.munge.txt
#
# #Run LDSC
# python2 /home/maihofer/trauma_gwas/ldsc-master/ldsc.py \
# --h2 results_filtered/eur_ptsd_pcs_v4_nov15_2021_${sex}.gz1.tbl.munge.gz.sumstats.gz \
# --ref-ld-chr /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
# --w-ld-chr  /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
# --out results_filtered/eur_ptsd_pcs_v4_nov15_2021_${sex}.gz1.tbl.munge.gz.tbl.ldsc
