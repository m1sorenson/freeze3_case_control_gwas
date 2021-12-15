# Freeze3 case/control GWAS meta-analysis
This repository contains code to run the case/control meta-analysis on PGC PTSD freeze 3. This is different from the overall meta-analysis in a few ways:
1.	Each study is ran as a case/control study, meaning that continuous phenotype information is converted into a 1/0 diagnosis of PTSD
2.	Study weights are adjusted to reflect case/control effective N, so continuous studies do not have as much weight

## Directory structure & file descriptions

```
freeze3_gwas
├──config - study weights, study codes, and genotyped gwas configuration
├──errandout – error and output of the slurm scripts
├──metal_results
│   ├──eur_ptsd_pcs_v4_chr{CHR}_nov15_2021_1.tbl – metal results for chromosome CHR
│   └──eur_ptsd_pcs_v4_{RUNTYPE}_allcase_nov15_2021.tbl – concatenated results of all chromosomes for runtype broad/males/females
├──metal_scripts – metal scripts built from the relevant metal_templates file, split by chromosome
├──metal_templates
│   ├──f3_gwas_broad_allcase.mi – semi-final template for case/control analysis (studies with low N's may need to be commented out)
│   └──f3_gwas_case_control.mi – base template for case/control analysis
├──pheno
│   ├──p2_{study1}_{ancestry}_{study2}_pcs.cov – covariate file for study1_study2 for individuals of the specified ancestry
│   └──p2_{study1}_{study2}_{ancestry}.pheno – phenotype file for study1_study2 for individuals of the specified ancestry
├──results_cat – contains gwas results of each study we have genotype data for (result of step 1)
├──results_filtered – final meta-analysis results filtered to only SNPs contained in a weighted >80% of the studies (weighted by how large the studies are)
└──sumstats – summary statistics for studies we don’t have genotype data
    └──bychr – summary statistics split by chromosome
```

## Prerequisites
### Software/Programs
-	plink2 installed and in PATH variable
-	metal installed and in PATH variable (or hard coded path in run_meta_v2_loo_v2.slurm script)

### Files

1. must have read access to all the freeze 2/3 summary statistics in the DAC directory
2. must have read access to all the genotyped freeze3 studies in the DAC directory
3. 01_f3_gwas_v1_allcase.sh – this file contains all the commands to run the GWAS for genotyped studies and the meta-analysis across all studies (though it cannot be run all at once as a script as each code block contains SLURM jobs that need to finish before the next block can be run)
4. config:
   - dosage_locations_f3.csv – this file contains study codes and ancestry codes, and determines which studies to be analyzed and included in the meta-analysis (only VETS is excluded because the GWAS is done using BOLT LMM rather than plink)
   - study_weights_allcase.tsv – this file contains the effective N for all the case/control studies we have summary statistics for (and for studies with a differing allele frequency column for male/female/all summary stats, the file also includes the respective names of the allele frequency column)
   - sumstat_studies.tsv – this file keeps track of which studies with summary statistics we are including in each analysis type (broad, male, female), this is important because some studies don’t have enough males/females to be included in one of the sex-stratified analyses, or we don’t have any sex-stratified summary statistics for the study.
5. f3_gwas_case_control.mi – this file is a template METAL file containing placeholders for male/female-specific file extensions, frequency column names, and effective N (weight). These placeholders are replaced with their values found in study_weights_allcase.tsv to create a semi-final METAL file f3_gwas_broad_allcase.mi (it is semi-final because some studies may need to be commented out if they should be excluded from the respective meta-analysis due to missing/low N)
6. pheno folder – this folder was copied from /home/maihofer/freeze3_gwas, and it contains the phenotype and covariate info for each genotyped individual

## Usage
### STEP 1: Run the GWAS step for genotyped studies

1. Review the studies in config/dosage_locations_f3.csv and make sure each study that should be included has the value in the exclude column set to 0 (note that each of the included studies must have a related fam file inside the datadir defined in run_trauma_gwas_v2_freeze3_allcase.slurm as well as .pheno and .cov files in the pheno directory)
2. Run the header block of the 01_f3_gwas_v1_allcase.sh script that sets the RUNTYPE and sex variables
3. In the 01_f3_gwas_v1_allcase.sh script, look for the code block that begins with ### 1)Study Level Analysis steps and run it.
4. Wait for the slurm jobs for each study to finish running (about 3 hours)

### STEP 2: Split the summary stat files by chromosome

1. Run the next block in the 01_f3_gwas_v1_allcase.sh script, starting with ### 2) Summary stat splitting
2. Wait for the slurm jobs for each study to finish running (about 10-20 min)

### STEP 3: Generate and fix the template metal file

1. Run the next block in the script, starting with ### 3) Metal script template generation
2. If you are doing a male/female-specific meta-analysis, adjust the generated f3_gwas_[males/females]_allcase.mi file in the metal_templates folder by doing the following:
   - for each study that has no male/female-specific summary statistics, comment or delete the study out of the meta-analysis
   - for each study that has too few (less than 45) males/females, comment or delete the study out of the meta-analysis (do this for both summary statistics & genotyped datasets)
   - double check the study weights in config/study_weights_allcase.tsv and make sure they were filled in the .mi file correctly

### STEP 4: Run the meta-analysis (METAL)

1. Run the next block in the script, starting with ### 4) Split METAL file by chromosome
2. Wait for the slurm job to finish, and check the error logs in errandout/f3_gwas_{females/males/broad}_chr10.mi_errorlogs

### STEP 5: Combine meta-analysis results and run FUMA

1. Run the next block in the script, starting with ### 5) Combine METAL results and generate final output. This filters the results to SNPs only contained in a weighted >80% of the studies (studies are weighted by their effective N).
2. Once this block is done, you should have the following files in the results_filtered folder:
   - .fuma.gz file – this should be uploaded to FUMA snp2gene
   - .tbl.premunge.gz file – this is the uncleaned version of sumstats that needs to be munged to be used for LDhub and LDSC
   - .tbl.munge.gz.sumstats.gz file – this is used for LDhub and to run LDSC
   - .tbl.munge.gz.tbl.ldsc file – this is the results from LDSC
