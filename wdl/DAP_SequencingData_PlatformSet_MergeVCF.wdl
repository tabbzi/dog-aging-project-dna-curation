version 1.0
# Namespace: DAP_SequencingData
# Workflow: DAP_SequencingData_PlatformSet_MergeVCF

# Curator: Kathleen Morrill (kmorrill@broadinstitute.org)
#   This workflow...
#   1. Merges VCFs from a set of platform entities
#      1.1. ...
#      1.2. ...

# Entity: platform_set

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

workflow DAP_SequencingData_PlatformSet_MergeVCF {
  input {
    Array[File] input_vcf
    Array[File] input_vcf_ind
    String set_id
    String? workflow_docker_image_override
    String workflow_docker_image = select_first([workflow_docker_image_override, "tabbzi/terra-dogagingproject-genomic-reports"])
    Float? genotype_probability_override
    Float genotype_probability = select_first([genotype_probability_override,0.7])
    Int? disk_pad_override
    Int disk_pad = select_first([disk_pad_override,20])
    Int? merge_disk_size_override
    Int merge_disk_size = select_first([merge_disk_size_override,(ceil(size(input_vcf,"GB"))*2)+disk_pad])
    Int? merge_mem_size_override
    Int merge_mem_size = select_first([merge_mem_size_override,8])
    Int? threads_override
    Int threads = select_first([threads_override,4])
  }

  call MergeData {
    input:
    mem_size = merge_mem_size,
    disk_size = merge_disk_size,
    num_threads = threads,
    docker_image = workflow_docker_image,
    filterGP = genotype_probability,
    setID = set_id,
    inputSetVCF=input_vcf,
    inputSetVCFind=input_vcf_ind
  }
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# TASKS

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

task MergeData {
  input {
    Int mem_size
    Int disk_size
    Int num_threads
    String docker_image
    Array[File] inputSetVCF
    Array[File] inputSetVCFind
    String filterGP
    String setID
  }

  command {
    echo ${sep=' ' inputSetVCF}
    ../bin/bcftools-1.9/bcftools merge --force-samples --threads ${num_threads} -Oz -o ${setID}_gp-${filterGP}.vcf.gz ${sep=' ' inputSetVCF}
    ../bin/bcftools-1.9/bcftools index -t '${setID}_gp-${filterGP}.vcf.gz'

    # make a PLINK2 set of multiallelic variants
    ../bin/plink2/plink2 --dog --vcf ${setID}_gp-${filterGP}.vcf.gz --make-pfile --out ${setID}_gp-${filterGP}

    # make a PLINK1 set of biallelic variants
    ../bin/plink2/plink2 --dog --pfile ${setID}_gp-${filterGP} --max-alleles 2 --snps-only just-acgt --set-all-var-ids '@:#:$r:$a' --make-bed --out ${setID}_gp-${filterGP}_biallelic_all-vars
  }

  runtime {
    docker: docker_image
    memory: "${mem_size} GB"
    cpu: "${num_threads}"
    preemptible: 1
    disks: "local-disk ${disk_size} HDD"
  }

  output {
    File mergedVCF = '${setID}_gp-${filterGP}.vcf.gz'
    File mergedVCFind = '${setID}_gp-${filterGP}.vcf.gz.tbi'
    File plinkPVAR = '${setID}_gp-${filterGP}.pvar'
    File plinkPSAM = '${setID}_gp-${filterGP}.psam'
    File plinkPGEN = '${setID}_gp-${filterGP}.pgen'
    File plinkBIM = '${setID}_gp-${filterGP}_biallelic_all-vars.bim'
    File plinkFAM = '${setID}_gp-${filterGP}_biallelic_all-vars.fam'
    File plinkBED = '${setID}_gp-${filterGP}_biallelic_all-vars.bed'
  }
}
