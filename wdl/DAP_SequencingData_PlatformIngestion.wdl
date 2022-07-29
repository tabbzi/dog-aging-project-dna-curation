# Namespace: DAP_SequencingData
# Workflow: DAP_SequencingData_PlatformIngestion

# WARNING: This workflow will only query data for a platform ID ('{gencove ID}') already in the platform table.
# WARNING: This workflow does not check if the data cannot be downloaded and the output may erase existing data model entries if no outputs found. Ensure correct platform ID is selected, data is available, and model is backed up BEFORE running.

# Curator: Kathleen Morrill (kmorrill@broadinstitute.org)
#   This workflow...
#   1. Downloads deliverables from the Gencove platform (task: GencoveAPI_Download)
#      1.1. Downloads data deliverables
#      1.2. Renames index files that should have the same prefix for genomics software to work on them (e.g. "Doggo_impute-vcf.vcf.gz" and "Doggo_impute-tbi.vcf.gz.tbi" will not pair in bcftools but "Doggo_impute-vcf.vcf.gz" and "Doggo_impute-vcf.vcf.gz.tbi" will, yet platform decided to name deliverables in the former manner).
#      1.3. If the sequencing run failed qc and only FASTQ files are available, then outputs will be NULL for other files

# Entity: platform

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

workflow DAP_SequencingData_PlatformIngestion {
  Int workflow_disk_size
  Int workflow_mem_size
  String? workflow_docker_image_override
  String workflow_docker_image = select_first([workflow_docker_image_override, "tabbzi/terra-dogagingproject-genomic-reports"])
  String? platform_id
  String? gencove_API_key_override
  String gencove_API_key = select_first([gencove_API_key_override, "${{ secrets.PLATFORM_PROJECT_API_KEY }}"])
  String? gencove_project_ID_override
  String gencove_project_ID = select_first([gencove_project_ID_override, "${{ secrets.PLATFORM_PROJECT_ID }}"])

  call GencoveAPI_Download {
    input:
    disk_size = workflow_disk_size,
    mem_size = workflow_mem_size,
    dockerImage = workflow_docker_image,
    keyAPI = gencove_API_key,
    projectID = gencove_project_ID,
    query_platform_id = platform_id
  }
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# TASKS

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

task GencoveAPI_Download {
  Int disk_size
  Int mem_size
  String dockerImage
  String keyAPI
  String projectID
  String query_platform_id
  String download_template = "'{gencove_id}_{file_type}.{file_extension}'"

  command {
    gencove download . --sample-ids ${query_platform_id} --api-key ${keyAPI} --file-types alignment-bam,alignment-bai,impute-vcf,impute-tbi,fastq-r1,fastq-r2 --download-template ${download_template}

    mv ${query_platform_id}_alignment-bai.bam.bai ${query_platform_id}_alignment-bam.bam.bai
    mv ${query_platform_id}_impute-tbi.vcf.gz.tbi ${query_platform_id}_impute-vcf.vcf.gz.tbi

  }

  runtime {
    docker: dockerImage
    memory: mem_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    continueOnReturnCode: [0, 1]
  }

  output {
    File? vcf = select_first(["${query_platform_id}_impute-vcf.vcf.gz","NULL"])
    File? vcf_index = select_first(["${query_platform_id}_impute-vcf.vcf.gz.tbi","NULL"])
    File? bam = select_first(["${query_platform_id}_alignment-bam.bam","NULL"])
    File? bam_index = select_first(["${query_platform_id}_alignment-bam.bam.bai","NULL"])
    File? fastqr1 = select_first(["${query_platform_id}_fastq-r1.fastq.gz","NULL"])
    File? fastqr2 = select_first(["${query_platform_id}_fastq-r2.fastq.gz","NULL"])
  }
}
