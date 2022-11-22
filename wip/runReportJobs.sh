#!/bin/bash

set -euo pipefail
# Set variables
# gcpath=${gcpath:-'gs://fc-6f3f8275-c9b4-4dcf-b2de-70f8d74f0874/2d80b263-7d01-485f-851c-2e684b13c450/DAP_SequencingData_MergeVCF/0f8bc383-9c6e-4f05-8fc7-b5934636e7f3/call-MergeData/DogAgingProject_DogAgingProject_2022_Merge_13_gp-0.7.p*'}
dir=${dir:-'testing/'} # this matters for ADMIXTURE which outputs to the working directory, which is annoying
input=${input:-"testing/DogAgingProject_2022-08-15_gp-0.7"}
# fasta=${fasta:-'gs://fc-6f3f8275-c9b4-4dcf-b2de-70f8d74f0874/ref/Canis_lupus_familiaris_assembly3.fasta'}
fasta=${fasta:-'testing/canis_assembly3.fasta'}
imgdir="/home/tcomi/projects/dog-aging-project-dna-curation/images"

# input: PLINK2 pfile set file prefix

# Obtain merged data (PLINK2 file set)
# gsutil -m cp ${gcpath} .

# plink2/dev PLINK2 Alpha 3.6 final (14 Aug) https://www.cog-genomics.org/plink/2.0/
# Generate PLINK1 file set of biallelic SNPs
# plink2/dev --dog --pfile ${input} --snps-only just-acgt --ref-from-fa --fa ${fasta} --max-alleles 2 --set-all-var-ids '@:#:$r:$a' --make-bed --out ${input}_biallelic-snps
singularity exec ${imgdir}/plink2_2.00a3.3--hb2a7ceb_0.sif \
plink2 \
  --dog \
  --pfile ${input} \
  --snps-only just-acgt \
  --ref-from-fa \
  --fa ${fasta} \
  --max-alleles 2 \
  --set-all-var-ids '@:#:$r:$a' \
  --make-bed \
  --out ${input}_biallelic-snps

# Generate VCF of simple trait prediction variants
# plink/dev --dog --pfile ${input} --extract bed0 VariantsOfInterest.bed --export vcf bgz --out ${input}_trait-predictions
singularity exec ${imgdir}/plink2_2.00a3.3--hb2a7ceb_0.sif \
plink2 \
  --dog \
  --pfile ${input}  \
  --extract bed0 /home/tcomi/projects/dog-aging-project-dna-curation/ref/VariantsOfInterest.bed \
  --export vcf bgz \
  --out ${input}_trait-predictions

# Submit PLINK1 bfile set for global ancestry inference using ADMIXTURE

./GlobalAncestry_ADM.sh
# qsub -v dir=${dir},input=${input}'_biallelic-snps' /seq/vgb/dap/bin/GlobalAncestry_ADM.sh

# Submit trait predictions

# Run simple traits:
# qsub /seq/vgb/dap/bin/runTraitPredictions.sh ${input}'_trait-predictions'
/tigress/tcomi/.conda/dap/bin/python genotype_to_phenotype.py \
  --vcf ${input}_trait-predictions.vcf.gz \
  --output testing/trait/DogAgingProject_2022-08-15_gp-0.7_trait-predictions \
  --variants /home/tcomi/projects/dog-aging-project-dna-curation/ref/DogAgingProject_VariantsOfInterest.csv \
  --reference CanFam3.1 \
  --imputation impute-v2


# Run complex traits:
export input=${input}
export mod=BodySize
./predict/run.sh
export mod=WhiteSpotting
# ./predict/run.sh

# Submit PLINK1 bfile set for inbreeding
# qsub -v input=${input} CoefficientOfInbreedingROH.sh

# SETTINGS DIFFER FROM DATA RELEASE (see curated data release README.md)
singularity exec ${imgdir}/plink_1.90b6.21--hec16e2b_2.sif \
plink --bfile ${input}_biallelic-snps \
  --dog \
  --chr 1-38 \
  --homozyg \
  --homozyg-density 10 \
  --homozyg-gap 500 \
  --homozyg-het 3 \
  --homozyg-kb 100 \
  --homozyg-window-threshold 0.10 \
  --out ${input}'.ROH'

autoLen=${autoLen:-2203765}
awk -v l=${autoLen} \
  'BEGIN{OFS=","; print "id,nSeg,kbTot,kbAvg,coi"}
   NR>1{print $2,$4,$5,$6,$5/l} ' \
  ${input}'.ROH.hom.indiv' \
  > ${input}'.ROH.COI.csv' 

# Generate genomic report JSONs for each sample
python3 ReportsToJSON.py \
  --study-ids <(awk '{print $2}' ${input}'_biallelic-snps.fam') \
  --ancestry 'testing/parseADM/StudyID-{id}.GlobalAncestry.json' \
  --inbreeding ${input}'.ROH.COI.csv' \
  --phenotypes ${input}'_trait-predictions_phenotypeTable.csv' \
  --genotypes ${input}'_trait-predictions_jsonTable.csv' \
  --body-size-genotypes ${input}'_predictionModel-BodySize.traw.BodySize.genotypes.csv' \
  --body-size-phenotypes ${input}'_predictionModel-BodySize.traw.BodySize.phenotypes.csv' \
  --white-spotting-genotypes ${input}'_predictionModel-WhiteSpotting.traw.WhiteSpotting.genotypes.csv' \
  --white-spotting-phenotypes ${input}'_predictionModel-WhiteSpotting.traw.WhiteSpotting.phenotypes.csv' \
  --output "testing/json/StudyID-{id}_GenomicReport.json"
