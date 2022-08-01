#!/bin/bash

# set up model
bin=${bin:-'/seq/vgb/dap/predict'}
mod=${mod:-'BodySize'}
snp=${bin}/${mod}.snp.bed
rscript=${bin}/predict.r
rdata=${bin}/Model.${mod}.RData
topsnp=${bin}/${mod}.topsnp.csv

echo ${bin}
echo ${mod}
echo ${snp}
echo ${rscript}
echo ${rdata}
echo ${topsnp}

# input PLINK2 file
input=${input:-''}

# make TRAW file
/seq/vgb/software/plink2/dev --dog --pfile ${input} --extract bed1 ${snp} --export A-transpose --out ${input}'_predictionModel-'${mod}

# run prediction
${rscript} ${input}'_predictionModel-'${mod}'.traw' ${rdata} ${topsnp} ${mod}
