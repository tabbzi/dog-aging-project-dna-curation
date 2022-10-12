
# Set variables
gcpath=${gcpath:-'gs://fc-6f3f8275-c9b4-4dcf-b2de-70f8d74f0874/2d80b263-7d01-485f-851c-2e684b13c450/DAP_SequencingData_MergeVCF/0f8bc383-9c6e-4f05-8fc7-b5934636e7f3/call-MergeData/DogAgingProject_DogAgingProject_2022_Merge_13_gp-0.7.p*'}
dir=${dir:-''} # this matters for ADMIXTURE which outputs to the working directory, which is annoying
input=${input:-"DogAgingProject_DogAgingProject_2022_Merge_13_gp-0.7"}
fasta=${fasta:-'gs://fc-6f3f8275-c9b4-4dcf-b2de-70f8d74f0874/ref/Canis_lupus_familiaris_assembly3.fasta'}

# input: PLINK2 pfile set file prefix

# Obtain merged data (PLINK2 file set)
gsutil -m cp ${gcpath} .

# Generate PLINK1 file set of biallelic SNPs
/seq/vgb/software/plink2/dev --dog --pfile ${input} --snps-only just-acgt --ref-from-fa --fa ${fasta} --max-alleles 2 --set-all-var-ids '@:#:$r:$a' --make-bed --out ${input}_biallelic-snps

# Generate VCF of simple trait prediction variants
/seq/vgb/software/plink2/dev --dog --pfile ${input} --extract bed0 VariantsOfInterest.bed --export vcf bgz --out ${input}_trait-predictions

# Submit PLINK1 bfile set for global ancestry inference using ADMIXTURE
qsub -v dir=${dir},input=${input}'_biallelic-snps' /seq/vgb/dap/bin/GlobalAncestry_ADM.sh

# Parse ADMIXURE outputs
python3 /seq/vgb/dap/bin/parseGlobalAncestryADM.py -S <(awk '{print $2}' ${input}'_biallelic-snps.fam') -Q ${input}'_biallelic-snps_GlobalAncestrySNPs_RefMerge.114.Q.labeled'

# Submit trait predictions

# Run simple traits:
qsub /seq/vgb/dap/bin/runTraitPredictions.sh ${input}'_trait-predictions'

# Run complex traits:
export input=${input}
export mod=BodySize
/seq/vgb/dap/predict/run.sh
export mod=WhiteSpotting
/seq/vgb/dap/predict/run.sh

# Submit PLINK1 bfile set for inbreeding
qsub -v input=${input} CoefficientOfInbreedingROH.sh

# Generate genomic report JSONs for each sample
while read S
do
  pytho3 ReportsToJSON.py -S ${S} -A 'StudyID.'${S}'.json' -I ${input}'.ROH.COI.csv' -P ${input}'_trait-predictions_phenotypeTable.csv' -G ${input}'_trait-predictions_jsonTable.csv' -SG ${input}'_predictionModel-BodySize.traw.BodySize.genotypes.csv' -SP ${input}'_predictionModel-BodySize.traw.BodySize.phenotypes.csv' -WG ${input}'_predictionModel-WhiteSpotting.traw.WhiteSpotting.genotypes.csv' -WP ${input}'_predictionModel-WhiteSpotting.traw.WhiteSpotting.phenotypes.csv'
done < <(awk '{print $2}' ${input}'_biallelic-snps.fam')