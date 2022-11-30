# Namespace: DAP_SequencingData
# Workflow: DAP_SequencingData_PlatformInference

# Curator: Kathleen Morrill (kmorrill@broadinstitute.org)
# Description:
#   This workflow...
#   1. Curate samples and sequencing runs (task: CurateProjectData)
#      1.1. Checks sequencing platform for status on sequencing runs (entity: platform_id)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

workflow DAP_SequencingData_PlatformInference {
  # inputs
  String base_filename
  String sample_name
  File reference_fasta = "canis_assembly3.fasta"
  String webhook_url
  File variants_of_interest = "VariantsOfInterest.bed"
  String plink_reference = "AncestryReferencePanel_GlobalAncestry"
  String population_key = "REDCap_Population_Key.tsv"
  File reference_variants = "DogAgingProject_VariantsOfInterest.csv"
  String predict_basename = "wip/predict"
  String webhook = ""  # will not signal if unset

  # containers
  String plink1_docker_image = "quay.io/biocontainers/plink:1.90b6.21--hec16e2b_2"
  String plink2_docker_image = "quay.io/biocontainers/plink2:2.00a3.3--hb2a7ceb_0"
  String admixture_docker_image = "quay.io/biocontainers/admixture:1.3.0--0"
  String python_docker_image = "tabbzi/terra-dogagingproject-genomic-reports"
  String simple_predict_docker_image = "tabbzi/terra-dogagingproject-simple-predict"
  String predict_docker_image = "tabbzi/terra-dogagingproject-complex-predict"

  # scripts
  File python_parseADM = "parseGlobalAncestryADM.py"
  File python_genotype_to_phenotype = "genotype_to_phenotype.py"
  File python_reports_to_json = "ReportsToJSON.py"

  call GenerateBiallilicSnps {
    input:
      dockerImage = plink2_docker_image,
      base_filename = base_filename,
      sample_name = sample_name,
      fasta = reference_fasta,
  }

  call GenerateTraitPredictions {
    input:
      dockerImage = plink2_docker_image,
      base_filename = base_filename,
      sample_name = sample_name,
      variants_of_interest = variants_of_interest,
  }

  call GenerateGlobalAncestrySNPs {
    input:
      dockerImage = plink2_docker_image,
      sample_name = sample_name,
      reference = plink_reference,
      input_bed = GenerateBiallilicSnps.bed,
      input_bim = GenerateBiallilicSnps.bim,
      input_fam = GenerateBiallilicSnps.fam,
  }

  call AncestryMergeReference {
    input:
      dockerImage = plink1_docker_image,
      sample_name = sample_name,
      reference = plink_reference,
      input_bed = GenerateGlobalAncestrySNPs.bed,
      input_bim = GenerateGlobalAncestrySNPs.bim,
      input_fam = GenerateGlobalAncestrySNPs.fam,
  }

  call AncestryInference  {
    input:
      dockerImage = admixture_docker_image,
      sample_name = sample_name,
      merge_bed = AncestryMergeReference.bed,
      merge_fam = AncestryMergeReference.fam,
      merge_bim = AncestryMergeReference.bim,
      merge_pop = AncestryMergeReference.pop,
  }

  call ParseAdmixture {
    input:
      dockerImage = python_docker_image,

      merge_fam = AncestryMergeReference.fam,
      snps_fam = GenerateBiallilicSnps.fam,
      admixed_q = AncestryInference.q,
      script = python_parseADM,
      population_key = population_key,
  }

  call SimpleTraitPrediction {
    input:
      dockerImage = simple_predict_docker_image,
      sample_name = sample_name,

      script = python_genotype_to_phenotype,
      trait_vcf = GenerateTraitPredictions.vcf,
      variants = reference_variants,
  }

  call ProduceTRAW as body_traw {
    input:
      dockerImage = plink2_docker_image,
      sample_name = sample_name,
      base_filename = base_filename,

      predict_basename = predict_basename,
      module = "BodySize",
  }

  call PredictComplex as body_predict {
    input:
      dockerImage = predict_docker_image,

      predict_basename = predict_basename,
      module = "BodySize",
      traw = body_traw.traw,
  }

  call ProduceTRAW as white_traw {
    input:
      dockerImage = plink2_docker_image,
      sample_name = sample_name,
      base_filename = base_filename,

      predict_basename = predict_basename,
      module = "WhiteSpotting",
  }

  call PredictComplex as white_predict {
    input:
      dockerImage = predict_docker_image,

      predict_basename = predict_basename,
      module = "WhiteSpotting",
      traw = white_traw.traw,
  }

  call GenerateInbreeding {
    input:
      dockerImage = plink1_docker_image,
      sample_name = sample_name,
      input_bed = GenerateBiallilicSnps.bed,
      input_bim = GenerateBiallilicSnps.bim,
      input_fam = GenerateBiallilicSnps.fam,
  }

  call SendReports {
    input:
      dockerImage = python_docker_image,
      script = python_reports_to_json,
      webhook = webhook,

      ancestry = ParseAdmixture.json,
      snps_fam = GenerateBiallilicSnps.fam,
      inbreeding = GenerateInbreeding.coi,
      phenotypes = SimpleTraitPrediction.phenotypes,
      genotypes = SimpleTraitPrediction.json_table,
      size_genotypes = body_predict.genotypes,
      size_phenotypes = body_predict.phenotypes,
      white_genotypes = white_predict.genotypes,
      white_phenotypes = white_predict.phenotypes,
  }


}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# TASKS

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

task GenerateBiallilicSnps {
  String dockerImage
  String base_filename
  String sample_name
  File fasta
  File pgen = "${base_filename}.pgen"
  File psam = "${base_filename}.psam"
  File pvar = "${base_filename}.pvar"

  command <<<
    plink2 \
      --dog \
      --pgen ${pgen} \
      --pvar ${pvar} \
      --psam ${psam} \
      --snps-only just-acgt \
      --ref-from-fa \
      --fa ${fasta} \
      --max-alleles 2 \
      --set-all-var-ids '@:#:$r:$a' \
      --make-bed \
      --out ${sample_name}_biallelic-snps
  >>>

  parameter_meta {
    base_filename: "Full path base name to .p files"
    fasta: "Reference fasta file"
  }

  runtime {
    docker: dockerImage
    disks: "local-disk 5 HDD"
    memory: "2G"
    CPU: 1
    preemptible: 1
  }

  output {
    File log = "${sample_name}_biallelic-snps.log"
    File fam = "${sample_name}_biallelic-snps.fam"
    File bim = "${sample_name}_biallelic-snps.bim"
    File bed = "${sample_name}_biallelic-snps.bed"
  }
}

task GenerateTraitPredictions {
  String dockerImage
  String base_filename
  String sample_name
  File variants_of_interest
  File pgen = "${base_filename}.pgen"
  File psam = "${base_filename}.psam"
  File pvar = "${base_filename}.pvar"

  command <<<
    plink2 \
      --dog \
      --pgen ${pgen} \
      --pvar ${pvar} \
      --psam ${psam} \
      --extract bed0 ${variants_of_interest} \
      --export vcf bgz \
      --out ${sample_name}_trait-predictions
  >>>

  parameter_meta {
    base_filename: "Full path base name to .p files"
    sample_name: "Prefix for outputs"
    variants_of_interest: "Reference variant sites"
  }

  runtime {
    docker: dockerImage
    disks: "local-disk 5 HDD"
    memory: "2G"
    CPU: 1
    preemptible: 1
  }

  output {
    File log = "${sample_name}_trait-predictions.log"
    File vcf = "${sample_name}_trait-predictions.vcf.gz"
  }
}

task GenerateGlobalAncestrySNPs {
  String dockerImage
  String sample_name
  String reference
  File reference_bim = reference + ".bim"
  File input_bed
  File input_bim
  File input_fam

  command <<<
    plink2 \
      --dog \
      --bed ${input_bed} \
      --bim ${input_bim} \
      --fam ${input_fam} \
      --extract ${reference_bim} \
      --make-bed \
      --out ${sample_name}_GlobalAncestrySNPs
  >>>

  parameter_meta {
    sample_name: "Prefix for outputs"
    reference_bim: "Reference bim file"
  }

  runtime {
    docker: dockerImage
    disks: "local-disk 5 HDD"
    memory: "2G"
    CPU: 1
    preemptible: 1
  }

  output {
    File log = "${sample_name}_GlobalAncestrySNPs.log"
    File fam = "${sample_name}_GlobalAncestrySNPs.fam"
    File bim = "${sample_name}_GlobalAncestrySNPs.bim"
    File bed = "${sample_name}_GlobalAncestrySNPs.bed"
  }
}

task AncestryMergeReference {
  String dockerImage
  String sample_name
  String reference

  File input_bed
  File input_bim
  File input_fam

  File ref_bed = reference + ".bed"
  File ref_bim = reference + ".bim"
  File ref_fam = reference + ".fam"

  command <<<
    plink \
      --dog \
      --bed ${input_bed} \
      --bim ${input_bim} \
      --fam ${input_fam} \
      --bmerge ${ref_bed} ${ref_bim} ${ref_fam} \
      --make-bed \
      --out ${sample_name}_GlobalAncestrySNPs_RefMerge

    # make .pop file
    awk \
      'NR==FNR{a[$1]=$1; next} {print a[$1]}' \
      ${ref_fam} \
      ${sample_name}_GlobalAncestrySNPs_RefMerge.fam | \
    sed -e 's/^$/-/' \
      > ${sample_name}_GlobalAncestrySNPs_RefMerge.pop
  >>>

  runtime {
    docker: dockerImage
    disks: "local-disk 5 HDD"
    memory: "2G"
    CPU: 1
    preemptible: 1
  }

  output {
    File log = "${sample_name}_GlobalAncestrySNPs_RefMerge.log"
    File fam = "${sample_name}_GlobalAncestrySNPs_RefMerge.fam"
    File bim = "${sample_name}_GlobalAncestrySNPs_RefMerge.bim"
    File bed = "${sample_name}_GlobalAncestrySNPs_RefMerge.bed"
    File pop = "${sample_name}_GlobalAncestrySNPs_RefMerge.pop"
    File nosex = "${sample_name}_GlobalAncestrySNPs_RefMerge.nosex"
  }
}

task AncestryInference {
  String dockerImage
  String sample_name
  String K = "114"

  File merge_bed
  File merge_fam
  File merge_bim
  File merge_pop

  command <<<
    admixture \
      --supervised ${merge_bed} ${K} \
      -j8
  >>>

  runtime {
    docker: dockerImage
    disks: "local-disk 5 HDD"
    memory: "16G"
    CPU: 8
    preemptible: 1
  }

  output {
    File q = "${sample_name}_GlobalAncestrySNPs_RefMerge.${K}.Q"
    File p = "${sample_name}_GlobalAncestrySNPs_RefMerge.${K}.P"
  }

}

task ParseAdmixture {
  String dockerImage

  File merge_fam
  File snps_fam
  File admixed_q
  File script
  File population_key

  command <<<
    # paste sample IDs to ADMIXTURE .Q output file
    cut -d' ' -f 1,2 \
      ${merge_fam} |\
    paste - ${admixed_q} \
      > 'Q.labeled'

    # parse outputs
    python3 ${script} \
      -P ${population_key} \
      -S <(cut -f 2 ${snps_fam}) \
      -Q 'Q.labeled' \
      -O ./

  >>>

  runtime {
    docker: dockerImage
    disks: "local-disk 5 HDD"
    memory: "2G"
    CPU: 1
    preemptible: 1
  }

  output {
    File q = "Q.labeled"
    File json = "GlobalAncestry.json"
    File csv = "GlobalAncestry.csv"
    File fullResults = "fullResults.GlobalAncestry.csv" 
  }

}

task SimpleTraitPrediction {
  String dockerImage
  String sample_name

  File script
  File trait_vcf
  File variants

  String reference = "CanFam3.1"
  String imputation = "impute-v2"

  command <<<
    python3 ${script} \
      --vcf ${trait_vcf} \
      --output ${sample_name}_trait-predictions \
      --variants ${variants} \
      --reference ${reference} \
      --imputation ${imputation}
  >>>

  runtime {
    docker: dockerImage
    disks: "local-disk 5 HDD"
    memory: "2G"
    CPU: 1
    preemptible: 1
  }

  output {
    File genotypes = "${sample_name}_trait-predictions_genotypeTable.csv"
    File json_table = "${sample_name}_trait-predictions_jsonTable.csv"
    File phenotypes = "${sample_name}_trait-predictions_phenotypeTable.csv"
    File tb_genotype = "${sample_name}_trait-predictions_trailblazerGenotypeTable.csv"
    File tb_phenotype = "${sample_name}_trait-predictions_trailblazerPhenotypeTable.tsv"
  }
}

task ProduceTRAW {
  String dockerImage
  String base_filename
  String sample_name
  String module
  String predict_basename

  File pgen = "${base_filename}.pgen"
  File psam = "${base_filename}.psam"
  File pvar = "${base_filename}.pvar"

  File snp_bed = "${predict_basename}/${module}.snp.bed"

  command <<<
    plink2 \
        --dog \
        --pgen ${pgen} \
        --pvar ${pvar} \
        --psam ${psam} \
        --extract bed1 ${snp_bed} \
        --export A-transpose \
        --out ${sample_name}'_predictionModel-'${module}
  >>>

  runtime {
    docker: dockerImage
    disks: "local-disk 5 HDD"
    memory: "2G"
    CPU: 1
    preemptible: 1
  }

  output {
    File log = "${sample_name}_predictionModel-${module}.log"
    File traw = "${sample_name}_predictionModel-${module}.traw"
  }
}

task PredictComplex {
  String dockerImage
  String module
  String predict_basename

  File traw
  File script = "${predict_basename}/predict.r"
  File r_data = "${predict_basename}/Model.${module}.RData"
  File top_snp = "${predict_basename}/${module}.topsnp.csv"

  command <<<
    ${script} \
        ${traw} \
        ${r_data} \
        ${top_snp} \
        ${module}
  >>>

  runtime {
    docker: dockerImage
    disks: "local-disk 5 HDD"
    memory: "2G"
    CPU: 1
    preemptible: 1
  }

  output {
    File genotypes = "${module}.genotypes.csv"
    File phenotypes = "${module}.phenotypes.csv"
  }
}

task GenerateInbreeding {
  String dockerImage
  String sample_name
  String autoLen = "2203765"

  File input_bed
  File input_bim
  File input_fam

  command <<<
    plink \
      --dog \
      --bed ${input_bed} \
      --bim ${input_bim} \
      --fam ${input_fam} \
      --chr 1-38 \
      --homozyg \
      --homozyg-density 10 \
      --homozyg-gap 500 \
      --homozyg-het 3 \
      --homozyg-kb 100 \
      --homozyg-window-threshold 0.10 \
      --out ${sample_name}'.ROH'

    awk -v l=${autoLen} \
      'BEGIN{OFS=","; print "id,nSeg,kbTot,kbAvg,coi"}
       NR>1{print $2,$4,$5,$6,$5/l} ' \
      ${sample_name}'.ROH.hom.indiv' \
      > ${sample_name}'.ROH.COI.csv' 

  >>>

  runtime {
    docker: dockerImage
    disks: "local-disk 5 HDD"
    memory: "2G"
    CPU: 1
    preemptible: 1
  }

  output {
    File log = "${sample_name}.ROH.log"
    File nosex = "${sample_name}.ROH.nosex"
    File hom_indiv = "${sample_name}.ROH.hom.indiv"
    File hom = "${sample_name}.ROH.hom"
    File hom_summary = "${sample_name}.ROH.hom.summary"
    File coi = "${sample_name}.ROH.COI.csv"
  }
}

task SendReports {
  String dockerImage
  File script

  String webhook

  File ancestry
  File snps_fam
  File inbreeding
  File phenotypes
  File genotypes
  File size_genotypes
  File size_phenotypes
  File white_genotypes
  File white_phenotypes

  command <<<
    python3 ${script} \
      --study-ids <(awk '{print $2}' ${snps_fam}) \
      --ancestry ${ancestry} \
      --inbreeding ${inbreeding} \
      --phenotypes ${phenotypes} \
      --genotypes ${genotypes} \
      --body-size-genotypes ${size_genotypes} \
      --body-size-phenotypes ${size_phenotypes} \
      --white-spotting-genotypes ${white_genotypes} \
      --white-spotting-phenotypes ${white_phenotypes} \
      --output "StudyID-{id}_GenomicReport.json" \
      --webhook "${webhook}"

  >>>

  runtime {
    docker: dockerImage
    disks: "local-disk 5 HDD"
    memory: "2G"
    CPU: 1
    preemptible: 1
  }

  output {
    Array[File] json = glob("StudyID-*_GenomicReport.json")
  }
}
