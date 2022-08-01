#!/bin/bash
#$ -q broad
#$ -l h_vmem=1g
#$ -l h_rt=1:00:00
#$ -o /seq/vgb/dap/log/
#$ -e /seq/vgb/dap/log/
PATH=$PATH:/seq/vgb/software/anaconda3/current
source activate /seq/vgb/dd/prediction/phenotypes/env/

/seq/vgb/dap/bin/pushToWebhook.sh ${json}
