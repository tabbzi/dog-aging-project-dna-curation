# Namespace: DAP_SequencingData
# Workflow: DAP_SequencingData_PlatformFilterVCF

# Curator: Kathleen Morrill (kmorrill@broadinstitute.org)
#   This workflow...
#   1. Filter individual VCFs at given genotyping probability (GP) (task: FilterData)
#      1.1. Normalizes indels and splits multi-allelic variants into biallelic records
#      1.2. Filters for max(GP) greater than given threshold (default: 0.7)
#      1.3. Renames internal ID from Swab ID (sample) to Study ID (participant)

# Entity: platform

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

workflow DAP_SequencingData_PlatformFilterVCF {
  File input_vcf
  File input_vcf_ind
  String input_sample
  String input_participant
  String? workflow_docker_image_override
  String workflow_docker_image = select_first([workflow_docker_image_override, "tabbzi/terra-dogagingproject-genomic-reports"])
  Float genotype_probability = 0.7

  call FilterData {
    input:
    docker_image = workflow_docker_image,
    querySample = input_sample,
    queryParticipant = input_participant,
    queryVCF = input_vcf,
    queryVCFind = input_vcf_ind,
    filterGP = genotype_probability
  }
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# TASKS

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

task FilterData {
  String docker_image
  File queryVCF
  File queryVCFind
  String querySample
  String queryParticipant
  String filterGP
  String filterExp = "'MAX(GP[*])<" + filterGP + "'"

  command <<<
    # filter by GP field, split multiallelic sites, and annotate variant IDs:
    echo "check in bin, is there a bin?"
    ls ../bin/
    ../bin/bcftools-1.9/bcftools view -e ${filterExp} -Ou ${queryVCF} | \
    ../bin/bcftools-1.9/bcftools norm -m-any -Ou | \
    ../bin/bcftools-1.9/bcftools annotate -x 'ID' -I '%CHROM:%POS:%REF:%ALT' -Oz -o ${querySample}.vcf.gz
    ../bin/bcftools-1.9/bcftools index -t ${querySample}.vcf.gz

    # change sample ID to study ID:
    echo ${queryParticipant} > ${queryParticipant}.txt
    ../bin/bcftools-1.9/bcftools reheader -s ${queryParticipant}.txt -o "${queryParticipant}_gp-${filterGP}.vcf.gz" ${querySample}.vcf.gz
    ../bin/bcftools-1.9/bcftools index -t "${queryParticipant}_gp-${filterGP}.vcf.gz"
  >>>

  runtime {
    docker: docker_image
    memory: "2 GB"
    disks: "local-disk " + ceil(2*size(queryVCF,"GB"))+20 + " HDD"
  }

  output {
    File vcf_filt = '${queryParticipant}_gp-${filterGP}.vcf.gz'
    File vcf_filt_index = '${queryParticipant}_gp-${filterGP}.vcf.gz.tbi'
  }
}
