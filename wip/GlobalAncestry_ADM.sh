#!/bin/bash
#$ -q broad
#$ -l h_vmem=8g
#$ -l h_rt=72:00:00
#$ -pe smp 10
#$ -R y
#$ -binding linear:10
umask 002

# source /broad/software/scripts/useuse
#   reuse UGER
#   reuse GCC-5.2
#   reuse .htslib-1.8
#   reuse Bcftools
#   reuse Tabix
#   reuse VCFtools
#   reuse R-3.5
#   reuse Anaconda
#   plink=/seq/vgb/software/plink2/current/plink
#   admixture=/seq/vgb/software/ADMIXTURE/current

set -euo pipefail
dir=${dir:-'testing/'}
input=${input:-'testing/DogAgingProject_2022-08-15_gp-0.7'}
K=${K:-114}
ref=${ref:-'testing/AncestryReferencePanel_GlobalAncestry'}
imgdir="/home/tcomi/projects/dog-aging-project-dna-curation/images"


# # /seq/vgb/software/plink2/dev --dog --bfile ${dir}/${input} --extract ${ref}'.bim' --make-bed --out ${dir}/${input}'_GlobalAncestrySNPs'
# singularity exec ${imgdir}/plink2_2.00a3.3--hb2a7ceb_0.sif \
# plink2 \
#   --dog \
#   --bfile ${input}_biallelic-snps \
#   --extract ${ref}.bim \
#   --make-bed \
#   --out ${input}_GlobalAncestrySNPs
#
# # /seq/vgb/software/plink2/current/plink --dog --bfile ${dir}/${input}'_GlobalAncestrySNPs' --bmerge ${ref} --make-bed --out ${dir}/${input}'_GlobalAncestrySNPs_RefMerge'
# singularity exec ${imgdir}/plink_1.90b6.21--hec16e2b_2.sif \
# plink \
#   --dog \
#   --bfile ${input}_GlobalAncestrySNPs \
#   --bmerge ${ref} \
#   --make-bed \
#   --out ${input}_GlobalAncestrySNPs_RefMerge
#
# # make .pop file
# awk \
#   'NR==FNR{a[$1]=$1; next} {print a[$1]}' \
#   ${ref}.fam \
#   ${input}_GlobalAncestrySNPs_RefMerge.fam | \
# sed -e 's/^$/-/' \
#   > ${input}_GlobalAncestrySNPs_RefMerge.pop

cd ${dir}

# perform global ancestry inference
# ${admixture} --supervised ${dir}/${input}'_GlobalAncestrySNPs_RefMerge.bed' ${K} -j10
singularity exec ${imgdir}/admixture_1.3.0--0.sif \
admixture \
  --supervised ${input##*/}'_GlobalAncestrySNPs_RefMerge.bed' ${K} \
  -j8

cd - > /dev/null

exit 0

# paste sample IDs to ADMIXTURE .Q output file
awk \
  '{print $1,$2}'\
  ${dir}/${input}'_GlobalAncestrySNPs_RefMerge.fam' |\
paste - ${input}'_GlobalAncestrySNPs_RefMerge.'${K}'.Q' \
  > ${dir}/${input}'_GlobalAncestrySNPs_RefMerge.'${K}'.Q.labeled'

# parse outputs
python3 /seq/vgb/dap/bin/parseGlobalAncestryADM.py \
  -S <(awk '{print $2}' ${dir}/${input}'.fam') \
  -Q ${dir}/${input}'_GlobalAncestrySNPs_RefMerge.'${K}'.Q.labeled'
