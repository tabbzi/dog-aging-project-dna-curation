#!/bin/bash
#$ -q broad
#$ -l h_vmem=8g
#$ -l h_rt=72:00:00
#$ -pe smp 10
#$ -R y
#$ -binding linear:10
umask 002

source /broad/software/scripts/useuse
  reuse UGER
  reuse GCC-5.2
  reuse .htslib-1.8
  reuse Bcftools
  reuse Tabix
  reuse VCFtools
  reuse R-3.5
  reuse Anaconda
  plink=/seq/vgb/software/plink2/current/plink
  admixture=/seq/vgb/software/ADMIXTURE/current

dir=${dir:-'/seq/vgb/dap/data/plink/bfiles/sets/'}
input=${input:-'DogAgingProject_extras_001_gp-0.7_biallelic_snps-only_BreedReferenceSNPs_BreedReferenceMerge'}
K=${K:-114}
ref=${ref:-'/seq/vgb/dd/ancestry/reference/AncestryReferencePanel/AncestryReferencePanel_GlobalAncestry'}
cd ${dir}


/seq/vgb/software/plink2/dev --dog --bfile ${dir}/${input} --extract ${ref}'.bim' --make-bed --out ${dir}/${input}'_GlobalAncestrySNPs'
/seq/vgb/software/plink2/current/plink --dog --bfile ${dir}/${input}'_GlobalAncestrySNPs' --bmerge ${ref} --make-bed --out ${dir}/${input}'_GlobalAncestrySNPs_RefMerge'

# make .pop file
awk 'NR==FNR{a[$1]=$1; next}{$1=a[$1]; print $1}' ${ref}'.fam' ${dir}/${input}'_GlobalAncestrySNPs_RefMerge.fam' | sed -e 's/^$/-/' > ${dir}/${input}'_GlobalAncestrySNPs_RefMerge.pop'
# awk '{print $1}' ${dir}/${input}'_GlobalAncestrySNPs_RefMerge.fam' | sed s/"0"/"-"/g > ${dir}/${input}'_GlobalAncestrySNPs_RefMerge.pop'

# perform global ancestry inference
  ${admixture} --supervised ${dir}/${input}'_GlobalAncestrySNPs_RefMerge.bed' ${K} -j10

# paste sample IDs to ADMIXTURE .Q output file
  awk '{print $1,$2}' ${dir}/${input}'_GlobalAncestrySNPs_RefMerge.fam' | paste - ${input}'_GlobalAncestrySNPs_RefMerge.'${K}'.Q' > ${dir}/${input}'_GlobalAncestrySNPs_RefMerge.'${K}'.Q.labeled'

# parse outputs
  python3 /seq/vgb/dap/bin/parseGlobalAncestryADM.py -S <(awk '{print $2}' ${dir}/${input}'.fam') -Q ${dir}/${input}'_GlobalAncestrySNPs_RefMerge.'${K}'.Q.labeled'
