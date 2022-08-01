#!/bin/bash
plink='../bin/plink2/plink2'

# input: PLINK1 file set of biallelic SNPs
input=${input:-''}
autoLen=${autoLen:-2203765}

# SETTINGS DIFFER FROM DATA RELEASE (see curated data release README.md)
${plink} --bfile ${input} \
         --dog \
         --chr 1-38 \
         --homozyg \
         --homozyg-density 10 \
         --homozyg-gap 500 \
         --homozyg-het 3 \
         --homozyg-kb 100 \
         --homozyg-window-threshold 0.10 \
         --out ${input}'.ROH'

echo "id,nSeg,kbTot,kbAvg,coi" > ${input}'.ROH.COI.csv'
tail -n+2 ${input}'.ROH.hom.indiv' | awk -v l=${autoLen} 'OFS=","{print $2,$4,$5,$6,$5/l}' >> ${input}'.ROH.COI.csv'
