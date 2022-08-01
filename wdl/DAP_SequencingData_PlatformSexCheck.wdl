# Namespace: DAP_SequencingData
# Workflow: DAP_SequencingData_PlatformSexCheck

# Curator: Kathleen Morrill (kmorrill@broadinstitute.org)
#   This workflow...
#   1. Checks coverage of genome by sequencing reads from BAM
#      1.1. If no BAM index, then index BAM (just to fix problem at start of project in which BAM index files not downloaded)
#      1.2. Calculate X chromosome coverage and divide by autosomal coverage
#      1.3. Calculate average autosomal coverage

# Entity: platform

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

workflow DAP_SequencingData_PlatformSexCheck {
  File? bam
  File? bam_index
  Boolean? indexed_override
  Boolean indexed = select_first([indexed_override,true])
  File? bash_checksex # bin/sex_xratio_bam.sh
  File? bash_autcov # bin/avg_aut_cov_bam.sh
  String? workflow_docker_image_override
  String workflow_docker_image = select_first([workflow_docker_image_override, "tabbzi/terra-dogagingproject-genomic-reports"])
  String? callback_docker_image_override
  String callback_docker_image = select_first([callback_docker_image_override, "tabbzi/terra-dogagingproject-genomic-reports-callback"])

  if (indexed) {
    call SexCheck {
      input:
        queryBAM = bam,
        queryBAMind = bam_index,
        dockerImage = workflow_docker_image,
        checkSex = bash_checksex,
        autCov = bash_autcov
    }
  }

  if (!indexed) {
    call IndexSexCheck {
    input:
      queryBAM = bam,
      dockerImage = workflow_docker_image,
        checkSex = bash_checksex,
        autCov = bash_autcov
    }
  }
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# TASKS

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

task SexCheck {
  File queryBAM
  File queryBAMind
  String dockerImage
  File checkSex
  File autCov
  String outbase = basename(queryBAM, ".bam")

  command {
    bash ${checkSex} ${queryBAM} > ${outbase}.sex.txt
    bash ${autCov} ${queryBAM} > ${outbase}.avg_aut_cov.txt
  }

  runtime {
    docker: dockerImage
    disks: "local-disk " + ceil(2*size(queryBAM,"GB"))+5 + " HDD"
    memory: "8 GB"
    CPU: 1
    preemptible: 1
    maxRetries: 1
  }

  output {
    Float sex_xratio_bam = read_string('${outbase}.sex.txt')
    Float avg_aut_cov_bam = read_string('${outbase}.avg_aut_cov.txt')
    File bam_index = '${queryBAMind}'
  }
}

task IndexSexCheck {
  File queryBAM
  String dockerImage
  File checkSex
  File autCov
  String outbase = basename(queryBAM, ".bam")

  command {
    ../bin/samtools-1.11/samtools index -b ${queryBAM}
    ls -lth
    bash ${checkSex} ${queryBAM} > ${outbase}.sex.txt
    bash ${autCov} ${queryBAM} > ${outbase}.avg_aut_cov.txt
  }

  runtime {
    docker: dockerImage
    disks: "local-disk " + ceil(2*size(queryBAM,"GB"))+5 + " HDD"
    memory: "8 GB"
    CPU: 1
    preemptible: 1
    maxRetries: 1
  }

  output {
    Float sex_xratio_bam = read_string('${outbase}.sex.txt')
    Float avg_aut_cov_bam = read_string('${outbase}.avg_aut_cov.txt')
    File bam_index = '${queryBAM}.bai'
  }
}
