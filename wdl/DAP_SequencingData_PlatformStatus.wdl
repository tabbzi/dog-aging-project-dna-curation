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
  File sample_kits_table # key: SampleKitTable
  String this_workspace_name
  String this_project_name
  String workflow_docker_image = "tabbzi/terra-dogagingproject-data-curation"
  String webhook_url
  String gencove_API_key
  String gencove_project_ID
  File python_UpdateDataModel # bin/update_data_model.py

  call CurateProjectData {
    input:
      dockerImage = workflow_docker_image,
      keyAPI = gencove_API_key,
      projectID = gencove_project_ID,
      workspaceName = this_workspace_name,
      projectName = this_project_name,
      webhookURL = webhook_url,
      kitsTable = sample_kits_table,
      pythonUpdateDataModel = python_UpdateDataModel
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
  String webhookURL
  File kitsTable
  File pythonUpdateDataModel

  command <<<
    # get sequencing run status table from platform:
    gencove projects list-samples \
      --api-key ${keyAPI} \
      ${projectID} >> statusTable.tsv

    # update data model using FireCloud API and post status to webhook:
    python3 ${pythonUpdateDataModel} \
      --workspace ${workspaceName} \
      --project ${projectName} \
      --status statusTable.tsv \
      --kittable ${kitsTable} \
      --webhook ${webhookURL}

  >>>

  parameter_meta {
    keyAPI: "API Key obtained for organization from Gencove platform at web.gencove.com"
    projectID: "ID for project on Gencove platform at web.gencove.com"
    workspaceName: "The sequencing data workspace"
    projectName: "The billing project hosting the sequencing data workspace"
    kitsTable: "TSV file with columns: entity:sample_id, dog_id, cohort, date_swab_arrival_laboratory, sample_type"
    pythonUpdateDataModel: "Python script for updating data models using Firecloud API"
    webhookURL: "The webhook URL to post JSON"
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
  }
}
