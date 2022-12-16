#!/bin/bash

# set up model
bin=${bin:-'predict'}
mod=${mod:-'BodySize'}
imgdir="/home/tcomi/projects/dog-aging-project-dna-curation/images"

echo ${mod}
input=${input:-''}

# make TRAW file
singularity exec ${imgdir}/plink2_2.00a3.3--hb2a7ceb_0.sif \
plink2 \
    --dog \
    --pfile ${input} \
    --extract bed1 ${bin}/${mod}.snp.bed \
    --export A-transpose \
    --out ${input}'_predictionModel-'${mod}

# run prediction
singularity exec ${imgdir}/predict.sif \
${bin}/predict.r \
    ${input}'_predictionModel-'${mod}'.traw' \
    ${bin}/Model.${mod}.RData \
    ${bin}/${mod}.topsnp.csv \
    ${mod}
