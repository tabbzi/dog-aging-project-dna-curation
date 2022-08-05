#!/bin/bash
#$ -q broad
#$ -l h_vmem=4g
#$ -l h_rt=1:00:00
#$ -M kmorrill@broadinstitute.org
#$ -o /seq/vgb/dap/log/
#$ -e /seq/vgb/dap/log/
PATH=$PATH:/seq/vgb/software/anaconda3/current
source activate /seq/vgb/dd/prediction/phenotypes/env/
cd /seq/vgb/dap/

# it would actually be better to run this on a dog's unfiltered VCF because then it can account for low confidence calls and won't miss variants

python3 /seq/vgb/dd/prediction/phenotypes/scripts/GenoToPhenoNew.py -F ${1}.vcf.gz -O ${1}
