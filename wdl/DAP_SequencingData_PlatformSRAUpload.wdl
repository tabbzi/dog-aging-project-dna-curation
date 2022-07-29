version 1.0
# Namespace: DAP_SequencingData
# Workflow: DAP_SequencingData_PlatformSRAUpload

# Curator: Kathleen Morrill (kmorrill@broadinstitute.org)
# Description:
#   This workflow...
#   1. Localize FASTQ files for a sequencing run
#   2. Rename files to {platformID}_{sampleID}_fastq-{r1/r2}.fastq.gz
#   3. Aspera transfer to NCBI SRA preload directory

# Entity: platform

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

workflow DAP_SequencingData_PlatformSRAUpload {
  input {
    File? aspera_keyfile_override
    File aspera_keyfile = select_first([aspera_keyfile_override, ""${{ secrets.SRA_ASPERA_KEY }}"])
    String? aspera_directory_override
    String aspera_directory = select_first([aspera_directory_override, "${{ secrets.SRA_ASPERA_DIR }}"])
    File input_fastqr1
    File input_fastqr2
    String platform_id
    String sample_id
    Int? disk_padding_override
    Int disk_padding = select_first([disk_padding_override,4])
    Int? disk_size_override
    Int disk_size = select_first([disk_size_override,ceil(size(input_fastqr1,"GB"))+ceil(size(input_fastqr2,"GB"))+disk_padding])
  }

  call AsperaUpload {
    input:
    fastqr1 = input_fastqr1,
    fastqr2 = input_fastqr2,
    platform = platform_id,
    sample = sample_id,
    asperaKeyFile = aspera_keyfile,
    asperaDirectory = aspera_directory,
    diskSize = disk_size
  }
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# TASKS

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

task AsperaUpload {
  input {
    File fastqr1
    File fastqr2
    String platform
    String sample
    File asperaKeyFile
    String asperaDirectory
    Int diskSize
  }

  command {
    mkdir upload
    mv ${fastqr1} upload/${platform}_${sample}_fastq-r1.fastq.gz
    mv ${fastqr2} upload/${platform}_${sample}_fastq-r2.fastq.gz
    ascp -i ${asperaKeyFile} -QT -l100m -k1 -d upload/ ${asperaDirectory}
  }

  runtime {
    docker: "tabbzi/terra-dogagingproject-data-curation:latest"
    disks: "local-disk ${diskSize} HDD"
    memory: "4G"
  }

  parameter_meta {
    asperaKeyFile: "Key file provided by NCBI SRA for use in Aspera CLI for data preloading"
    asperaDirectory: "Directory path provided by NCBI SRA for use in Aspera CLI for data preloading"
  }
}
