# Namespace: DAP_SequencingData
# Workflow: DAP_SequencingData_PlatformStatus

# Curator: Kathleen Morrill (kmorrill@broadinstitute.org)
# Description:
#   This workflow...
#   1. Curate samples and sequencing runs (task: CurateProjectData)
#      1.1. Checks sequencing platform for status on sequencing runs (entity: platform_id)
#      1.2. Checks project platform for sample assignments (entity: sample_id)
#      1.3. Curates updates to data model tables (participant, sample, platform)
#      1.4. Update `participant` data model with new participant_id entities
#      1.5. Update `sample` data model with new sample_id entities
#      1.6. Update `platform` data model to reflect sequencing run statuses
#      1.7. Generate status update JSONs
#   2. Notifies project platform on status updates (task: PostStatus)
#   By the end of this workflow, the data model will be ready to ingest data.

# Entity: none

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

workflow DAP_SequencingData_PlatformStatus {
  File? sample_kits_table # key: SampleKitTable
  String this_workspace_name
  String this_project_name
  String? workflow_docker_image_override
  String workflow_docker_image = select_first([workflow_docker_image_override, "tabbzi/terra-dogagingproject-data-curation"])
  String? callback_docker_image_override
  String callback_docker_image = select_first([callback_docker_image_override, "tabbzi/terra-dogagingproject-genomic-reports-callback"])
  String? webhook_url_override
  String webhook_url = select_first([webhook_url_override, "${{ secrets.DAP_WEBHOOK_URL }}"])
  String? gencove_API_key_override
  String gencove_API_key = select_first([gencove_API_key_override, "${{ secrets.PLATFORM_PROJECT_API_KEY }}"])
  String? gencove_project_ID_override
  String gencove_project_ID = select_first([gencove_project_ID_override, "${{ secrets.PLATFORM_PROJECT_ID }}"])
  File python_GetDataModel # bin/get_data_model.py
  File rscript_SetDataModel # bin/set_data_model.R
  File python_UpdateDataModel # bin/update_data_model.py

  call CurateProjectData {
    input:
      dockerImage = workflow_docker_image,
      keyAPI = gencove_API_key,
      projectID = gencove_project_ID,
      dockerImage = workflow_docker_image,
      workspaceName = this_workspace_name,
      projectName = this_project_name,
      kitsTable = sample_kits_table,
      pythonGetDataModel = python_GetDataModel,
      rscriptSetDataModel = rscript_SetDataModel,
      pythonUpdateDataModel = python_UpdateDataModel
  }

  scatter (J in CurateProjectData.jsons) {
    call PostStatus {
      input:
        dockerImage = callback_docker_image,
        thisJSON = J,
        webhookURL = webhook_url
    }
  }
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# TASKS

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

task CurateProjectData {
  String dockerImage
  String keyAPI
  String projectID
  String workspaceName
  String projectName
  File kitsTable
  File pythonGetDataModel
  File rscriptSetDataModel
  File pythonUpdateDataModel

  command <<<
    # get sequencing run status table from platform:
    gencove projects list-samples --api-key ${keyAPI} ${projectID} >> statusTable.tsv

    # get data model tables from workspace using FireCloud API:
    python3 ${pythonGetDataModel} -W ${workspaceName} -P ${projectName}

    # set data model tables with updated data:
    Rscript ${rscriptSetDataModel} -g statusTable.tsv -p platformTable.tsv -k ${kitsTable} -s sampleTable.tsv

    # update data model using FireCloud API and generate JSONs for status updates:
    python3 ${pythonUpdateDataModel} -W ${workspaceName} -P ${projectName} -A participant.tsv -B sample.tsv -C platform.tsv

  >>>

  parameter_meta {
    keyAPI: "API Key obtained for organization from Gencove platform at web.gencove.com"
    projectID: "ID for project on Gencove platform at web.gencove.com"
    workspaceName: "The sequencing data workspace"
    projectName: "The billing project hosting the sequencing data workspace"
    kitsTable: "TSV file with columns: entity:sample_id, dog_id, cohort, date_swab_arrival_laboratory, sample_type"
    pythonGetDataModel: "Python script for getting current data models using Firecloud API"
    rscriptSetDataModel: "R script for make data model tables"
    pythonUpdateDataModel: "Python script for updating data models using Firecloud API"
  }

  runtime {
    docker: dockerImage
    disks: "local-disk 5 HDD"
    memory: "2G"
    CPU: 1
    preemptible: 1
    maxRetries: 1
  }

  output {
    File participant_table = "participant.tsv"
    File sample_table = "sample.tsv"
    File platform_table = "platform.tsv"
    Array[File]+ jsons = glob("*.JSON")
  }
}

task PostStatus {
  String dockerImage
  String webhookURL
  File thisJSON

  command { /opt/callback.sh "${webhookURL}" "${thisJSON}"}

  output {
        String response_stdout = read_string(stdout())
        String response_stderr = read_string(stderr())
  }

  runtime {
    docker: dockerImage
    disks: "local-disk 5 HDD"
    memory: "2G"
    CPU: 1
    preemptible: 1
    maxRetries: 1
  }

  parameter_meta {
    webhookURL: "The webhook URL to post JSON"
  }
}
