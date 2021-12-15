#!/bin/bash
#-------------------------------------------------------------------------------
# SLURM command to run
#-------------------------------------------------------------------------------
# sbatch --time=24:00:00 --error errandout/copy_sumstats.e --output errandout/copy_sumstats.o 00_copy_sumstats.sh

#-------------------------------------------------------------------------------
# Description
#-------------------------------------------------------------------------------
# This script copies all of the summary statistic files from the pgc DAC account
# into a local sumstats directory (originally created to be run on pgca1pts on
# the LISA supercomputer system)

# change this to be the directory you are running the gwas in (sumstats will be
# copied to the "sumstats" folder inside this working directory)
WORKING_DIR=/home/pgca1pts/freeze3_gwas
cd $WORKING_DIR

# don't change this
ss_dir=${WORKING_DIR}/sumstats
mkdir -p ${ss_dir}
DAC_DIR=/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts



#-------------------------------------------------------------------------------
# VETS & UKBB
#-------------------------------------------------------------------------------
# VETS & UKBB - run these two studies with BOLT LMM
# VETS - ptsd_qt_vetsa_may12_2021_related_filtered.imputed.stats.gz
# UKBB - pts_ukbb_may13_2021_unrelated.bgen.stats.gz

#-------------------------------------------------------------------------------
# VETS
#-------------------------------------------------------------------------------
# Study only has males - use same file for both
echo VETS
dir=/home/pgca1pts/freeze3_bolt_data
zcat ${dir}/ptsd_qt_vetsa_may12_2021_related_filtered.imputed.stats.gz \
  | awk 'BEGIN{OFS="\t"}{if (NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/VETS_broad_eur.txt.gz
# male gwas
zcat ${dir}/ptsd_qt_vetsa_may12_2021_related_filtered.imputed.stats.gz \
  | awk 'BEGIN{OFS="\t"}{if (NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/VETS_males_eur.txt.gz
# broad case control gwas
zcat ${dir}/ptsd_qt_vetsa_may12_2021_related_filtered.imputed.stats.gz | awk 'BEGIN{OFS="\t"}{prev=77/(1012+77); if(NR==1){$11="OR"; print}if (NR>1 && $7 > 0.20 && $7 < 0.80) {$11=exp($11 /( prev *(1-prev))); $12 = $12  /( prev *(1-prev)); print} }' \
  | gzip > ${ss_dir}/VETS_broad_eur_allcase.txt.gz

#-------------------------------------------------------------------------------
# UKBB
#-------------------------------------------------------------------------------
# Study only has males - use same file for both
echo UKBB
dir=/home/maihofer/freeze3_gwas/sumstats
# base gwas
cp ${dir}/pts_ukbb_may13_2021_unrelated.bgen.stats.gz ${TMPDIR}/UKBB_broad_eur.txt.gz
zcat ${TMPDIR}/UKBB_broad_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/UKBB_broad_eur.txt.gz
# male gwas
dir=/home/pgca1pts/freeze3_bolt_data
cp ${dir}/pts_ukbb_may13_2021_unrelated_males.bgen.stats.gz \
  ${TMPDIR}/UKBB_males_eur.txt.gz
zcat ${TMPDIR}/UKBB_males_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/UKBB_males_eur.txt.gz
# female gwas
cp ${dir}/pts_ukbb_may13_2021_unrelated_females.bgen.stats.gz \
  ${TMPDIR}/UKBB_females_eur.txt.gz
zcat ${TMPDIR}/UKBB_females_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/UKBB_females_eur.txt.gz
# broad case control gwas
zcat ${dir}/pts_ukbb_may13_2021_unrelated.bgen.stats.gz | awk '{prev=10913/(124888+10913); if(NR==1) $11="OR"; if (NR>1) {$11=exp($11 /( prev *(1-prev))); $12 = $12  /( prev *(1-prev))}; print }' \
  | gzip > ${ss_dir}/UKBB_broad_eur_allcase.txt.gz


#-------------------------------------------------------------------------------
# GROUP 3
#-------------------------------------------------------------------------------
# GROUP 3 sumstats - freeze 2 summary stats


#-------------------------------------------------------------------------------
# MIRE
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
echo MIRE
dir=${DAC_DIR}/wave2/summary_stats/by_study/mire
# base gwas
cp ${dir}/MIRE_eur_analysis1_mf.gz ${TMPDIR}/MIRE_broad_eur.txt.gz
zcat ${TMPDIR}/MIRE_broad_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/MIRE_broad_eur.txt.gz
# male gwas
cp ${dir}/MIRE_eur_analysis1_males.gz ${TMPDIR}/MIRE_males_eur.txt.gz
zcat ${TMPDIR}/MIRE_males_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/MIRE_males_eur.txt.gz
# female gwas
cp ${dir}/MIRE_eur_analysis1_females.gz ${TMPDIR}/MIRE_females_eur.txt.gz
zcat ${TMPDIR}/MIRE_females_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/MIRE_females_eur.txt.gz
# clean TMPDIR
rm ${TMPDIR}/MIRE*

#-------------------------------------------------------------------------------
# INTR
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
echo INTR
dir=${DAC_DIR}/wave2/summary_stats/by_study/intr
# base gwas
cp ${dir}/INTr_eur_analysis1_mf.gz ${TMPDIR}/INTR_broad_eur.txt.gz
zcat ${TMPDIR}/INTR_broad_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/INTR_broad_eur.txt.gz
# male gwas
cp ${dir}/INTr_eur_analysis1_males.gz ${TMPDIR}/INTR_males_eur.txt.gz
zcat ${TMPDIR}/INTR_males_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/INTR_males_eur.txt.gz
# female gwas - no females file for this study
# clean TMPDIR
rm ${TMPDIR}/INTR*

#-------------------------------------------------------------------------------
# DAMI
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
echo DAMI
dir=${DAC_DIR}/wave2/summary_stats/by_study/dami
# base gwas
cp ${dir}/daner_psd_25July.gz ${TMPDIR}/DAMI_broad_eur.txt.gz
zcat ${TMPDIR}/DAMI_broad_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/DAMI_broad_eur.txt.gz
# male gwas
cp ${dir}/daner_pst_male.gz ${TMPDIR}/DAMI_males_eur.txt.gz
zcat ${TMPDIR}/DAMI_males_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/DAMI_males_eur.txt.gz
# female gwas
cp ${dir}/daner_pst_female.gz ${TMPDIR}/DAMI_females_eur.txt.gz
zcat ${TMPDIR}/DAMI_females_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/DAMI_females_eur.txt.gz
# clean TMPDIR
rm ${TMPDIR}/DAMI*

#-------------------------------------------------------------------------------
# QIMR
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
echo QIMR
dir=${DAC_DIR}/wave2/summary_stats/by_study/qimr
# base gwas
cp ${dir}/pts_qimr_mix_nm.logscale.results.gz ${TMPDIR}/QIMR_broad_eur.txt.gz
zcat ${TMPDIR}/QIMR_broad_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/QIMR_broad_eur.txt.gz
# male gwas
cp ${dir}/pts_qimr_mix_nm_males.logscale.results.gz ${TMPDIR}/QIMR_males_eur.txt.gz
zcat ${TMPDIR}/QIMR_males_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/QIMR_males_eur.txt.gz
# female gwas
cp ${dir}/pts_qimr_mix_nm_females.logscale.results.gz ${TMPDIR}/QIMR_females_eur.txt.gz
zcat ${TMPDIR}/QIMR_females_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/QIMR_females_eur.txt.gz
# clean TMPDIR
rm ${TMPDIR}/QIMR*

#-------------------------------------------------------------------------------
# NCPT
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
echo NCPT
dir=${DAC_DIR}/wave2/summary_stats/by_study/ncpt
# base gwas
cp ${dir}/N800_eur_analysis1_mf.gz ${TMPDIR}/NCPT_broad_eur.txt.gz
zcat ${TMPDIR}/NCPT_broad_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/NCPT_broad_eur.txt.gz
# male gwas
cp ${dir}/N800_eur_analysis1_males.gz ${TMPDIR}/NCPT_males_eur.txt.gz
zcat ${TMPDIR}/NCPT_males_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/NCPT_males_eur.txt.gz
# female gwas
cp ${dir}/N800_eur_analysis1_females.gz ${TMPDIR}/NCPT_females_eur.txt.gz
zcat ${TMPDIR}/NCPT_females_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/NCPT_females_eur.txt.gz
# clean TMPDIR
rm ${TMPDIR}/NCPT*

#-------------------------------------------------------------------------------
# TRAC
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
echo TRAC
dir=${DAC_DIR}/wave2/summary_stats/by_study/trac
# base gwas
cp ${dir}/TRAC_eur_analysis1_mf.gz ${TMPDIR}/TRAC_broad_eur.txt.gz
zcat ${TMPDIR}/TRAC_broad_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/TRAC_broad_eur.txt.gz
# male gwas
cp ${dir}/TRAC_eur_analysis1_males.gz ${TMPDIR}/TRAC_males_eur.txt.gz
zcat ${TMPDIR}/TRAC_males_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/TRAC_males_eur.txt.gz
# female gwas - no females for this study
# clean TMPDIR
rm ${TMPDIR}/TRAC*


#-------------------------------------------------------------------------------
# GROUP 4
#-------------------------------------------------------------------------------
# GROUP 4 sumstats - freeze 3 summary stats


#-------------------------------------------------------------------------------
# WTCS
#-------------------------------------------------------------------------------
# Annotate RS IDs
# Filter to minor allele frequency >= 0.01
echo WTCS
dir=${DAC_DIR}/wave3/summary_stats/75_wtcs
# base gwas
cp ${dir}/test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear ${TMPDIR}
cp ${dir}/test_wtc_dosage_HDS_SamClean_MAF.afreq ${TMPDIR}

# male gwas
cp ${dir}/test_wtc_dosage_HDS_SamClean_MAF_M.PHENO1.glm.linear ${TMPDIR}
cp ${dir}/test_wtc_dosage_HDS_SamClean_MAF_M.afreq ${TMPDIR}

# female gwas
cp ${dir}/test_wtc_dosage_HDS_SamClean_MAF_F.PHENO1.glm.linear ${TMPDIR}
cp ${dir}/test_wtc_dosage_HDS_SamClean_MAF_F.afreq ${TMPDIR}

# Annotate RS IDs
for run in "broad" "males" "females"; do
  echo annotating $run
  if [ $run == "broad" ]; then
    base=test_wtc_dosage_HDS_SamClean_MAF
    outfile=WTCS_broad_eur.txt.gz
  elif [ $run == "males" ]; then
    base=test_wtc_dosage_HDS_SamClean_MAF_M
    outfile=WTCS_males_eur.txt.gz
  else
    base=test_wtc_dosage_HDS_SamClean_MAF_F
    outfilie=WTCS_females_eur.txt.gz
  fi
  LC_ALL=C join -1 2 -2 3 ${TMPDIR}/${base}.afreq ${TMPDIR}/${base}.PHENO1.glm.linear | sed 's/#//g' \
    | awk '{if(NR>1){ FRQ=$5; A1=$11;} if(A1==$10) { A2=$9 }  else if(A1 != $10){ A2=$10; FRQ=1-$5}; if(NR==1) {A2="A2"; FRQ="FRQ";} ; print $0,A2,FRQ}' \
    | awk '{if(NR == 1 || ($19 >= 0.01 && $19 <= 0.99)){print}}' > ${TMPDIR}/${base}.PHENO1.glm.linear.frq
  cat vcf_header1b.txt vcf_header2.txt <(awk '{OFS="\t"}{if (NR>1) print $7,$8, ".", $3,$4 , "100", "PASS", "MVP="$1}' ${TMPDIR}/${base}.PHENO1.glm.linear.frq \
    | LC_ALL=C sort -g -k1r,1 -k2,2 ) \
    > ${TMPDIR}/${base}.PHENO1.glm.linear.frq.vcf
  bcftools view ${TMPDIR}/${base}.PHENO1.glm.linear.frq.vcf -Oz -o ${TMPDIR}/${base}.PHENO1.glm.linear.frq.vcf.gz
  tabix -f ${TMPDIR}/${base}.PHENO1.glm.linear.frq.vcf.gz
  bcftools annotate -a /home/pgca1pts/All_20180423_hg19.vcf.gz -c INFO ${TMPDIR}/${base}.PHENO1.glm.linear.frq.vcf.gz -o ${TMPDIR}/${base}.PHENO1.glm.linear.frq.vcf.gz.annotated
  echo "ID SNP"  | awk 'BEGIN {OFS="\t"} {print $1,$2}' > ${TMPDIR}/${base}_snpheader.txt
  tail -n+82 ${TMPDIR}/${base}.PHENO1.glm.linear.frq.vcf.gz.annotated  | grep RS \
    | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$2}' \
    | sed 's/RS=/rs/g' | sed 's/MVP=//g'  | cat ${TMPDIR}/${base}_snpheader.txt -  | LC_ALL=C sort -k1b,1 \
    > ${TMPDIR}/${base}.PHENO1.glm.linear.frq.vcf.gz.annotated.success
  tail -n+82 ${TMPDIR}/${base}.PHENO1.glm.linear.frq.vcf.gz.annotated | grep -v RS \
    | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$1}' \
    | sed 's/MVP=//g' | cat ${TMPDIR}/${base}_snpheader.txt - | LC_ALL=C sort -k1b,1 \
    > ${TMPDIR}/${base}.PHENO1.glm.linear.frq.vcf.gz.annotated.failed
  LC_ALL=C join -1 1 -2 1  ${TMPDIR}/${base}.PHENO1.glm.linear.frq.vcf.gz.annotated.success \
    <(LC_ALL=C sort -k1b,1 ${TMPDIR}/${base}.PHENO1.glm.linear.frq) \
    > ${TMPDIR}/${base}.PHENO1.glm.linear.frq.success
  LC_ALL=C join -1 1 -2 1  ${TMPDIR}/${base}.PHENO1.glm.linear.frq.vcf.gz.annotated.failed \
    <(LC_ALL=C sort -k1b,1 ${TMPDIR}/${base}.PHENO1.glm.linear.frq) \
    > ${TMPDIR}/${base}.PHENO1.glm.linear.frq.failed
  echo "annotation success:"
  wc -l ${TMPDIR}/${base}.PHENO1.glm.linear.frq.success
  echo "annotation failed:"
  wc -l ${TMPDIR}/${base}.PHENO1.glm.linear.frq.failed
  cat ${TMPDIR}/${base}.PHENO1.glm.linear.frq.success ${TMPDIR}/${base}.PHENO1.glm.linear.frq.failed \
    | sort -g -k 18  | awk '{if($20 >= 0.01 && $20 <= 0.99) print $2,$1,$8,$9,$13,$12,$19,$20,$14,$15,$16,$17,$18}' \
    | grep -v NA | gzip > ${TMPDIR}/${base}.PHENO1.glm.linear.maf01.gz
  echo "SNP ID CHROM POS TEST A1 A2 FRQ OBS_CT BETA SE T_STAT P" > ${TMPDIR}/${base}_header_x.txt
  zcat ${TMPDIR}/${base}.PHENO1.glm.linear.maf01.gz | cat ${TMPDIR}/${base}_header_x.txt - \
    | gzip > ${ss_dir}/${outfile}
done
# clean TMPDIR
base=test_wtc_dosage_HDS_SamClean_MAF
rm ${TMPDIR}/${base}*

#-------------------------------------------------------------------------------
# AGDS
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
echo AGDS
dir=${DAC_DIR}/wave3/summary_stats/86_93_agds_qim2
# base gwas
cp ${dir}/PTSDsum_AGDS_full_19052021.QIMRB.txt ${TMPDIR}/AGDS_broad_eur.txt
awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' ${TMPDIR}/AGDS_broad_eur.txt \
  | gzip > ${ss_dir}/AGDS_broad_eur.txt.gz
# male gwas
cp ${dir}/PTSDsum_AGDS_M_190521.QIMRB.zip ${TMPDIR}
unzip -d ${TMPDIR} ${TMPDIR}/PTSDsum_AGDS_M_190521.QIMRB.zip
mv ${TMPDIR}/PTSDsum_AGDS_M_19052021.QIMRB.txt ${TMPDIR}/AGDS_males_eur.txt
awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' ${TMPDIR}/AGDS_males_eur.txt \
  | gzip > ${ss_dir}/AGDS_males_eur.txt.gz
# female gwas
cp ${dir}/PTSDsum_AGDS_F_190521.QIMRB.zip ${TMPDIR}
unzip -d ${TMPDIR} ${TMPDIR}/PTSDsum_AGDS_F_190521.QIMRB.zip
mv ${TMPDIR}/PTSDsum_AGDS_F_19052021.QIMRB.txt ${TMPDIR}/AGDS_females_eur.txt
awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' ${TMPDIR}/AGDS_females_eur.txt \
  | gzip > ${ss_dir}/AGDS_females_eur.txt.gz
# clean TMPDIR
rm ${TMPDIR}/AGDS*


#-------------------------------------------------------------------------------
# CANA
#-------------------------------------------------------------------------------
# Already filtered to minor allele frequency >= 0.01
echo CANA
dir=${DAC_DIR}/wave3/summary_stats/92_cana
# base gwas
cp ${dir}/clsa_1-23xm2_psd_dctoff_com_out_covars_c4_20210512.PSD_DCTOFF_COM.assoc.logistic.add.gz \
  ${ss_dir}/CANA_broad_eur.txt.gz
# male gwas - no sex specific
#cp ${dir}/clsa_1-23xm2_covars_c4_20210512m.PSD_DCTOFF_COM.assoc.logistic.add.gz \
#  ${ss_dir}/CANA_males_eur.txt.gz
# female gwas - no sex specific
#cp ${dir}/clsa_1-23xm2_covars_c4_20210512f.PSD_DCTOFF_COM.assoc.logistic.add.gz \
#  ${ss_dir}/CANA_females_eur.txt.gz

#-------------------------------------------------------------------------------
# QIM2
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
echo QIM2
dir=${DAC_DIR}/wave3/summary_stats/86_93_agds_qim2
# base gwas
cp ${dir}/PTSDyn_others_full_19052021.QIMRB.txt ${TMPDIR}/QIM2_broad_eur.txt
awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' ${TMPDIR}/QIM2_broad_eur.txt \
  | gzip > ${ss_dir}/QIM2_broad_eur.txt.gz
# male gwas
cp ${dir}/PTSDyn_others_M_190521.QIMRB.zip ${TMPDIR}
unzip -d ${TMPDIR} ${TMPDIR}/PTSDyn_others_M_190521.QIMRB.zip
mv ${TMPDIR}/PTSDyn_others_M_19052021.QIMRB.txt ${TMPDIR}/QIM2_males_eur.txt
awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' ${TMPDIR}/QIM2_males_eur.txt \
  | gzip > ${ss_dir}/QIM2_males_eur.txt.gz
# female gwas
cp ${dir}/PTSDyn_others_F_190521.QIMRB.zip ${TMPDIR}
unzip -d ${TMPDIR} ${TMPDIR}/PTSDyn_others_F_190521.QIMRB.zip
mv ${TMPDIR}/PTSDyn_others_F_19052021.QIMRB.txt ${TMPDIR}/QIM2_females_eur.txt
awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' ${TMPDIR}/QIM2_females_eur.txt \
  | gzip > ${ss_dir}/QIM2_females_eur.txt.gz
# clean TMPDIR
rm ${TMPDIR}/QIM2*

#-------------------------------------------------------------------------------
# RCOG
#-------------------------------------------------------------------------------
# Already filtered to minor allele frequency >= 0.01, need SE
echo RCOG
dir=${DAC_DIR}/wave3/summary_stats/94_rcog
# base gwas
cp ${dir}/PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz ${TMPDIR}/RCOG_broad_eur.txt.gz
# male gwas - no sex specific
# female gwas - no sex specific

#-------------------------------------------------------------------------------
# WTCM
#-------------------------------------------------------------------------------
# Annotate minor allele frequency, set $1=$1 to fix weird tabs
echo WTCM
dir=/home/maihofer/freeze3_gwas/sumstats
# base gwas
cp ${dir}/PGC_WTC.W.assoc.logistic.merged.gwas2.gz ${ss_dir}/WTCM_broad_eur.txt.gz
# male gwas - no sex specific
# female gwas - no sex specific
# clean TMPDIR
rm ${TMPDIR}/WTCM*

#-------------------------------------------------------------------------------
# DAI2
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
echo DAI2
dir=${DAC_DIR}/wave3/summary_stats/74_dai2
# base gwas
cp ${dir}/daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.gz \
  ${TMPDIR}/DAI2_broad_eur.txt.gz
zcat ${TMPDIR}/DAI2_broad_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}'  \
  | gzip > ${ss_dir}/DAI2_broad_eur.txt.gz
# male gwas
cp ${dir}/daner_iPSYCH2015_PTSDbroad_males_chr1to22_HRC_MAF01.gz \
  ${TMPDIR}/DAI2_males_eur.txt.gz
zcat ${TMPDIR}/DAI2_males_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/DAI2_males_eur.txt.gz
# female gwas
cp ${dir}/daner_iPSYCH2015_PTSDbroad_females_chr1to22_HRC_MAF01.gz \
  ${TMPDIR}/DAI2_females_eur.txt.gz
zcat ${TMPDIR}/DAI2_females_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/DAI2_females_eur.txt.gz
# clean TMPDIR
rm ${TMPDIR}/DAI2*

#-------------------------------------------------------------------------------
# BIOV
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
# PLINK coding adjustment: If A1 is equal to ALT, use REF as A2. IF A1 is equal to REF, use ALT as A2. MUST also change the AF
echo BIOV
dir=${DAC_DIR}/wave3/summary_stats/79_biov_v2/
cp ${dir}/slavi_biov_sum_stats_w_freq.zip ${TMPDIR}
unzip -d ${TMPDIR} ${TMPDIR}/slavi_biov_sum_stats_w_freq.zip
mv ${TMPDIR}/ptsd_gwas_dosage_chr.all.PHENO2.glm.firth.additive.with.header.and.freqs \
  ${TMPDIR}/BIOV_broad_eur.txt
mv ${TMPDIR}/ptsd_gwas_dosage_chr.all.PHENO6.glm.firth.additive.with.header.and.freqs \
  ${TMPDIR}/BIOV_males_eur.txt
mv ${TMPDIR}/ptsd_gwas_dosage_chr.all.PHENO4.glm.firth.additive.with.header.and.freqs \
  ${TMPDIR}/BIOV_females_eur.txt
# base gwas
awk '{if(NR == 1 || ($13 >= 0.01 && $13 <= 0.99)){print}}' ${TMPDIR}/BIOV_broad_eur.txt \
  > ${TMPDIR}/BIOV_broad_eur_maf01.txt
awk '{if (NR == 1){print $0,"A2","A1_Freq"} else if ($6==$5){print $0,$4,$13} else if ($6==$4){print $0,$5,1-$13}}' \
  ${TMPDIR}/BIOV_broad_eur_maf01.txt | gzip > ${ss_dir}/BIOV_broad_eur.txt.gz
# male gwas
awk '{if(NR == 1 || ($13 >= 0.01 && $13 <= 0.99)){print}}' ${TMPDIR}/BIOV_males_eur.txt \
  > ${TMPDIR}/BIOV_males_eur_maf01.txt
awk '{if (NR == 1){print $0,"A2","A1_Freq"} else if ($6==$5){print $0,$4,$13} else if ($6==$4){print $0,$5,1-$13}}' \
  ${TMPDIR}/BIOV_males_eur_maf01.txt | gzip > ${ss_dir}/BIOV_males_eur.txt.gz
# female gwas
awk '{if(NR == 1 || ($13 >= 0.01 && $13 <= 0.99)){print}}' ${TMPDIR}/BIOV_females_eur.txt \
  > ${TMPDIR}/BIOV_females_eur_maf01.txt
awk '{if (NR == 1){print $0,"A2","A1_Freq"} else if ($6==$5){print $0,$4,$13} else if ($6==$4){print $0,$5,1-$13}}' \
  ${TMPDIR}/BIOV_females_eur_maf01.txt | gzip > ${ss_dir}/BIOV_females_eur.txt.gz
# clean TMPDIR
rm ${TMPDIR}/BIOV*

#-------------------------------------------------------------------------------
# MGBB
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
echo MGBB
dir=${DAC_DIR}/wave3/summary_stats/81_mgbb
# base gwas
cp ${dir}/pbk_eur_ptsd_gwas_broad_share.txt ${TMPDIR}/MGBB_broad_eur.txt
awk '{if(NR == 1 || ($6 >= 0.01 && $6 <= 0.99)){print}}' ${TMPDIR}/MGBB_broad_eur.txt \
  | gzip > ${ss_dir}/MGBB_broad_eur.txt.gz
# male gwas
cp ${dir}/pbk_eur_ptsd_gwas_broad_m_share.txt ${TMPDIR}/MGBB_males_eur.txt
awk '{if(NR == 1 || ($6 >= 0.01 && $6 <= 0.99)){print}}' ${TMPDIR}/MGBB_males_eur.txt \
  | gzip > ${ss_dir}/MGBB_males_eur.txt.gz
# female gwas
cp ${dir}/pbk_eur_ptsd_gwas_broad_f_share.txt ${TMPDIR}/MGBB_females_eur.txt
awk '{if(NR == 1 || ($6 >= 0.01 && $6 <= 0.99)){print}}' ${TMPDIR}/MGBB_females_eur.txt \
  | gzip > ${ss_dir}/MGBB_females_eur.txt.gz
# clean TMPDIR
rm ${TMPDIR}/MGBB*

#-------------------------------------------------------------------------------
# HUNT
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
# Annotate RS IDs
echo HUNT
dir=${DAC_DIR}/wave3/summary_stats/85_hunt
# base gwas
cp ${dir}/HUNT_ptsd_pgc_allchr_filter_info_results.txt.gz \
  ${TMPDIR}/HUNT_broad_eur.txt.gz
zcat ${TMPDIR}/HUNT_broad_eur.txt.gz | awk '{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' \
  > ${TMPDIR}/HUNT_broad_eur.txt

for run in "broad"; do
  echo annotating $run
  if [ $run == "broad" ]; then
    base=HUNT_broad_eur
  elif [ $run == "males" ]; then
    base=HUNT_males_eur
  else
    base=HUNT_females_eur
  fi
  cat vcf_header1b.txt vcf_header2.txt <(awk '{OFS="\t"}{if (NR>1) print $1,$2, ".", $4,$5 , "100", "PASS", "MVP="$3}' ${TMPDIR}/${base}.txt \
    | LC_ALL=C sort -g -k1r,1 -k2,2 ) > ${TMPDIR}/${base}.vcf
  bcftools view ${TMPDIR}/${base}.vcf -Oz -o ${TMPDIR}/${base}.vcf.gz
  tabix -f ${TMPDIR}/${base}.vcf.gz
  bcftools annotate -a /home/pgca1pts/All_20180423_hg19.vcf.gz -c INFO ${TMPDIR}/${base}.vcf.gz -o ${TMPDIR}/${base}.vcf.annotated
  # header start with # so it stays at top after sorting
  echo -e "#SNPID\tSNP" > ${TMPDIR}/${base}_snpheader.txt
  tail -n+82 ${TMPDIR}/${base}.vcf.annotated  | grep RS \
    | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$2}' \
    | sed 's/RS=/rs/g' | sed 's/MVP=//g'  | cat ${TMPDIR}/${base}_snpheader.txt -  | LC_ALL=C sort -k1b,1 \
    > ${TMPDIR}/${base}.vcf.annotated.success
  tail -n+82 ${TMPDIR}/${base}.vcf.annotated | grep -v RS \
    | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$1}' \
    | sed 's/MVP=//g' | cat ${TMPDIR}/${base}_snpheader.txt - | LC_ALL=C sort -k1b,1 \
    > ${TMPDIR}/${base}.vcf.annotated.failed
  LC_ALL=C join -1 1 -2 3  ${TMPDIR}/${base}.vcf.annotated.success \
    <(sed s/SNPID/\#SNPID/g ${TMPDIR}/${base}.txt | LC_ALL=C sort -k3b,3 ) \
    > ${TMPDIR}/${base}.txt.success
  LC_ALL=C join -1 1 -2 3  ${TMPDIR}/${base}.vcf.annotated.failed \
    <(sed s/SNPID/\#SNPID/g ${TMPDIR}/${base}.txt | LC_ALL=C sort -k3b,3 ) \
    > ${TMPDIR}/${base}.txt.failed
  echo "annotation success:"
  wc -l ${TMPDIR}/${base}.txt.success
  echo "annotation failed:"
  wc -l ${TMPDIR}/${base}.txt.failed
  cat ${TMPDIR}/${base}.txt.success ${TMPDIR}/${base}.txt.failed \
    | LC_ALL=C sort -g -k 3  | awk '{if(NR == 2  || ($8 >= 0.01 && $8 <= 0.99)) print $2,$1,$3,$4,$5,$6,$8,$9,$10,$11,$12,$13}' \
    | grep -v NA | sed s/\#SNPID/SNPID/g | gzip > ${ss_dir}/${base}.txt.gz
done
# clean TMPDIR
rm ${TMPDIR}/HUNT*

#-------------------------------------------------------------------------------
# SWED
#-------------------------------------------------------------------------------
# Already filtered to minor allele frequency >= 0.01
echo SWED
dir=${DAC_DIR}/wave3/summary_stats/89_swed/SWE_STAGE_PTSD_sumstats
# base gwas
cp ${dir}/SWE_STAGE_PTSD_saige_info_sumstats.gz ${ss_dir}/SWED_broad_eur.txt.gz
# male gwas - no sex specific
# female gwas - no sex specific

#-------------------------------------------------------------------------------
# FING
#-------------------------------------------------------------------------------
# Annotate RS IDs, liftover to hg19
# Filter to minor allele frequency >= 0.01
echo FING
dir=${DAC_DIR}/wave3/summary_stats/90_fing
# base gwas
cp ${dir}/file_download_22122020_chia_yen_PTSD_broad_bothsex.gz \
  ${TMPDIR}/FING_broad_eur.txt.gz
zcat ${TMPDIR}/FING_broad_eur.txt.gz | awk '{if(NR == 1 || ($8 >= 0.01 && $8 <= 0.99)){print}}' \
  | gzip > ${TMPDIR}/FING_broad_eur_maf01.txt.gz
mv ${TMPDIR}/FING_broad_eur_maf01.txt.gz ${TMPDIR}/FING_broad_eur.txt.gz
# male gwas
cp ${dir}/file_download_22122020_chia_yen_PTSD_broad_male.gz \
  ${TMPDIR}/FING_males_eur.txt.gz
zcat ${TMPDIR}/FING_males_eur.txt.gz | awk '{if(NR == 1 || ($8 >= 0.01 && $8 <= 0.99)){print}}' \
  | gzip > ${TMPDIR}/FING_males_eur_maf01.txt.gz
mv ${TMPDIR}/FING_males_eur_maf01.txt.gz ${TMPDIR}/FING_males_eur.txt.gz
# female gwas
cp ${dir}/file_download_22122020_chia_yen_PTSD_broad_female.gz \
  ${TMPDIR}/FING_females_eur.txt.gz
zcat ${TMPDIR}/FING_females_eur.txt.gz | awk '{if(NR == 1 || ($8 >= 0.01 && $8 <= 0.99)){print}}' \
  | gzip > ${TMPDIR}/FING_females_eur_maf01.txt.gz
mv ${TMPDIR}/FING_females_eur_maf01.txt.gz ${TMPDIR}/FING_females_eur.txt.gz

for run in "broad" "males" "females"; do
  echo annotating $run
  if [ $run == "broad" ]; then
    base=FING_broad_eur
  elif [ $run == "males" ]; then
    base=FING_males_eur
  else
    base=FING_females_eur
  fi
  zcat ${TMPDIR}/${base}.txt.gz | awk '{SNPID=$3; SNPID2=$3; gsub ("_"," ",SNPID) ;gsub ("_",":",SNPID2) ;  if(NR==1) print $0,"SNPID", "CHR","POS","A1","A2"; else print $0, SNPID2,SNPID}' > ${TMPDIR}/${base}.fixed.txt
  cat vcf_header1b.txt vcf_header2.txt <(awk '{OFS="\t"}{if (NR>1) print $26,$27, ".", $28,$29 , "100", "PASS", "MVP="$25}' ${TMPDIR}/${base}.fixed.txt \
    | sed 's/chr//g' | sort -n -k1r,1 -k2,2  ) > ${TMPDIR}/${base}.vcf
  bcftools view ${TMPDIR}/${base}.vcf -Oz -o ${TMPDIR}/${base}.vcf.gz
  tabix -f ${TMPDIR}/${base}.vcf.gz
  #Some will not be annotated, need to identify these for my merge steps..
  bcftools annotate -a /home/pgca1pts/All_20180418_hg38.vcf.gz -c INFO ${TMPDIR}/${base}.vcf.gz -o ${TMPDIR}/${base}.vcf.annotated
  echo "SNPID SNP"  | awk 'BEGIN {OFS="\t"} {print $1,$2}' > ${TMPDIR}/${base}_snpheader.txt
  tail -n+82 ${TMPDIR}/${base}.vcf.annotated  | grep RS | awk '{print $8}' \
    | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$2}' | sed 's/RS=/rs/g' | sed 's/MVP=//g' \
    | cat ${TMPDIR}/${base}_snpheader.txt -  | LC_ALL=C sort -k1b,1 > ${TMPDIR}/${base}.vcf.annotated.success
  tail -n+82 ${TMPDIR}/${base}.vcf.annotated | grep -v RS | awk '{print $8}' \
    | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$1}' | sed 's/MVP=//g' | cat ${TMPDIR}/${base}_snpheader.txt - \
    | LC_ALL=C sort -k1b,1   > ${TMPDIR}/${base}.vcf.annotated.failed
  wc -l ${TMPDIR}/${base}.fixed.txt
  wc -l ${TMPDIR}/${base}.vcf.annotated
  wc -l ${TMPDIR}/${base}.vcf.annotated.success
  wc -l ${TMPDIR}/${base}.vcf.annotated.failed
  LC_ALL=C join -1 1 -2 25  <(awk '{print "chr"$1,$2}' ${TMPDIR}/${base}.vcf.annotated.success \
    | cat ${TMPDIR}/${base}_snpheader.txt -)  <(LC_ALL=C sort -k25b,25 ${TMPDIR}/${base}.fixed.txt) \
    > ${TMPDIR}/${base}.success

  #We are at the point where all markers are annotated, time to lift over.
  #Convert to UCSC bed format
  #I am only interseted in the rs-id markers
  #Note that I do -1 for the starting position of the marker, this is what was done in  # BEGIN{OFS="\t"} {if(NR==1) print "chrom","chromStart","chromEnd";  #chrID is some junk that i accidentally merged
  tail -n+83 ${TMPDIR}/${base}.vcf.annotated.success  | sed 's/:/ /g' | awk '{print "chr"$1,$2-1,$2,$5}' | grep -v chrSNPID > ${TMPDIR}/${base}.annotated.bed

  #Liftover positions
  /home/pgca1pts/michael_liftover/liftOver ${TMPDIR}/${base}.annotated.bed \
    /home/pgca1pts/michael_liftover/hg38ToHg19.over.chain.gz \
    ${TMPDIR}/${base}.annotated.bed.lift \
    ${TMPDIR}/${base}.annotated.bed.unmapped
  #liftover rsid - not going to lift back - if an rs-id was deleted, its not real - rsids seem to be updated within builds also , so no telling which is used
  #merge in the new positions
  echo "CHRNEW POSNEWMIN1 BPNEW SNP" > ${TMPDIR}/${base}_liftheader.txt
  LC_ALL=C join -1 4 -2 2 <(cat ${TMPDIR}/${base}_liftheader.txt ${TMPDIR}/${base}.annotated.bed.lift \
    | LC_ALL=C sort -k 4b,4) <(LC_ALL=C sort -k2b,2 ${TMPDIR}/${base}.success ) \
    > ${TMPDIR}/${base}.success_lifted
  wc -l ${TMPDIR}/${base}.success_lifted

  #Get rid of 'chr' prefix..
  # file_download_22122020_chia_yen_PTSD_broad_male_july122021.gz.fuma.gz
  cat ${TMPDIR}/${base}.success_lifted | sort -g -k 23 \
    | awk '{if(NR == 1 || ($13 >= 0.01 && $13 <= 0.99)) print $1,$2,$4,$10,$11,$13,$18,$19,$20,$21,$22,$23}' \
    | sed 's/chr//g' | gzip > ${ss_dir}/${base}.txt.gz
done
# clean TMPDIR
rm ${TMPDIR}/FING*

#-------------------------------------------------------------------------------
# UKB2
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
echo UKB2
dir=${DAC_DIR}/wave3/summary_stats/91_ukb2
# base gwas
cp ${dir}/PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_WG.txt.regenie.gz \
  ${TMPDIR}/UKB2_broad_eur.txt.gz
zcat ${TMPDIR}/UKB2_broad_eur.txt.gz | awk '{if(NR == 1 || ($6 >= 0.01 && $6 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/UKB2_broad_eur.txt.gz
# male gwas
cp ${dir}/PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_Male_WG.txt.regenie.gz \
  ${TMPDIR}/UKB2_males_eur.txt.gz
zcat ${TMPDIR}/UKB2_males_eur.txt.gz | awk '{if(NR == 1 || ($6 >= 0.01 && $6 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/UKB2_males_eur.txt.gz
# female gwas
cp ${dir}/PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_Female_WG.txt.regenie.gz \
  ${TMPDIR}/UKB2_females_eur.txt.gz
zcat ${TMPDIR}/UKB2_females_eur.txt.gz | awk '{if(NR == 1 || ($6 >= 0.01 && $6 <= 0.99)){print}}' \
  | gzip > ${ss_dir}/UKB2_females_eur.txt.gz
# clean TMPDIR
rm ${TMPDIR}/UKB2*

#-------------------------------------------------------------------------------
# BIOM
#-------------------------------------------------------------------------------
echo BIOM
dir=${DAC_DIR}/wave3/summary_stats/95_biom
# base gwas
gzip -c ${dir}/PTSD_broad_EA_clean.head.assoc.logistic > ${TMPDIR}/BIOM_broad_eur.txt.gz
# male gwas
zcat ${dir}/PTSD_broad_EA_MALE.assoc.logistic.gz | awk '{if(NR == 1 || $5 == "ADD"){print}}' \
  | gzip > ${TMPDIR}/BIOM_males_eur.txt.gz
# female gwas
zcat ${dir}/PTSD_broad_EA_FEMALE.assoc.logistic.gz | awk '{if(NR == 1 || $5 == "ADD"){print}}' \
  | gzip > ${TMPDIR}/BIOM_females_eur.txt.gz
# Annotate minor allele frequencies
# SNP A1 A2 FRQ_U_###
zcat ${ss_dir}/DAI2_broad_eur.txt.gz | awk '{print $2,$4,$5,$7}' | LC_ALL=C sort -k1b,1 \
  > ${TMPDIR}/DAI2_eur_afs.txt
for run in "broad" "males" "females"; do
  if [ $run == "broad" ]; then
    base=BIOM_broad_eur
  elif [ $run == "males" ]; then
    base=BIOM_males_eur
  elif [ $run == "females" ]; then
    base=BIOM_females_eur
  fi
  LC_ALL=C join -1 1 -2 2 ${TMPDIR}/DAI2_eur_afs.txt <(zcat ${TMPDIR}/${base}.txt.gz | LC_ALL=C sort -k2b,2 ) \
    > ${TMPDIR}/${base}_merged.txt
  awk '{FRQ=$4; A1=$7; if(A1==$2) { A2=$3}  else if(A1 != $2){ A2=$2; FRQ=1-$4}; print $0,A2,FRQ}' \
    ${TMPDIR}/${base}_merged.txt | grep -v -E "A T$|T A$|C G$|G C$" > ${TMPDIR}/${base}_fixed.txt
  awk '{if(NR == 1 || ($8 == "ADD")){print $1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}}' \
    ${TMPDIR}/${base}_fixed.txt | gzip > ${ss_dir}/${base}.txt.gz
done
# clean TMPDIR
rm ${TMPDIR}/BIOM*

#-------------------------------------------------------------------------------
# ESBB
#-------------------------------------------------------------------------------
dir=${DAC_DIR}/wave3/summary_stats/98_esbb
# base gwas
cp ${dir}/PTSD_broad_EstBB_GWAS_results.txt.gz ${TMPDIR}/ESBB_broad_eur.txt.gz
zcat ${TMPDIR}/ESBB_broad_eur.txt.gz  | awk '{if (NR==1){$8="infNA"}; print}' \
  | awk '{if(NR==1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' | gzip > ${ss_dir}/ESBB_broad_eur.txt.gz
# male gwas
cp ${dir}/PTSD_broad_males_EstBB_GWAS_results.txt.gz ${TMPDIR}/ESBB_males_eur.txt.gz
zcat ${TMPDIR}/ESBB_males_eur.txt.gz  | awk '{if (NR==1){$8="infNA"}; print}' \
  | awk '{if(NR==1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' | gzip > ${ss_dir}/ESBB_males_eur.txt.gz
# female gwas
cp ${dir}/PTSD_broad_females_EstBB_GWAS_results.txt.gz ${TMPDIR}/ESBB_females_eur.txt.gz
zcat ${TMPDIR}/ESBB_females_eur.txt.gz  | awk '{if (NR==1){$8="infNA"}; print}' \
  | awk '{if(NR==1 || ($7 >= 0.01 && $7 <= 0.99)){print}}' | gzip > ${ss_dir}/ESBB_females_eur.txt.gz
# clean TMPDIR
rm ${TMPDIR}/ESBB*

#-------------------------------------------------------------------------------
# MAYO
#-------------------------------------------------------------------------------
dir=${DAC_DIR}/wave3/summary_stats/99_mayo
# base gwas
cp ${dir}/ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic.gz ${TMPDIR}/MAYO_broad_eur.txt.gz
zcat ${TMPDIR}/MAYO_broad_eur.txt.gz | awk '{if(NR==1){A2="A2"}; if(NR > 1 ){ if ($6==$5){A2=$4}; if($6==$4) {A2=$5}}; print $2,$3,$26,$13,$6,A2,$16,$17,$18,$19,$20,$21}' \
  | gzip > ${ss_dir}/MAYO_broad_eur.txt.gz
# male gwas
cp ${dir}/ptsd_def2.male.pcAdj.MAF0.01.glm.logistic.gz ${TMPDIR}/MAYO_males_eur.txt.gz
zcat ${TMPDIR}/MAYO_males_eur.txt.gz | awk '{if(NR==1){A2="A2"}; if(NR > 1 ){ if ($6==$5){A2=$4}; if($6==$4) {A2=$5}}; print $2,$3,$26,$13,$6,A2,$16,$17,$18,$19,$20,$21}' \
  | gzip > ${ss_dir}/MAYO_males_eur.txt.gz
# female gwas
cp ${dir}/ptsd_def2.female.pcAdj.MAF0.01.glm.logistic.gz ${TMPDIR}/MAYO_females_eur.txt.gz
zcat ${TMPDIR}/MAYO_females_eur.txt.gz | awk '{if(NR==1){A2="A2"}; if(NR > 1 ){ if ($6==$5){A2=$4}; if($6==$4) {A2=$5}}; print $2,$3,$26,$13,$6,A2,$16,$17,$18,$19,$20,$21}' \
  | gzip > ${ss_dir}/MAYO_females_eur.txt.gz
# clean TMPDIR
rm ${TMPDIR}/MAYO*

#-------------------------------------------------------------------------------
# GROUP 5
#-------------------------------------------------------------------------------
# GROUP 5 sumstats - MVP summary stats


#-------------------------------------------------------------------------------
# MVP
#-------------------------------------------------------------------------------
dir=${DAC_DIR}/wave3/summary_stats/82_mvp
# base gwas
cp ${dir}/dbGAP_totalPCL_eur.gz ${ss_dir}/MVP_broad_eur.txt.gz
# male gwas - no sex specific
# female gwas - no sex specific
