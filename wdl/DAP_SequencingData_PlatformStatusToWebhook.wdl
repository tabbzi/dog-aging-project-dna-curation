# Namespace: DAP_SequencingData
# Workflow: DAP_SequencingData_PlatformStatusToWebhook

# If status JSON failed to post to webhook, then run this workflow to generate and send a new JSON

# Entity: platform

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

workflow DAP_SequencingData_PlatformStatusToWebhook {
  String sample_id
  String participant_id
  String gencove_id
  String date_time
  String sample_status
  String? workflow_docker_image_override
  String workflow_docker_image = select_first([workflow_docker_image_override, "tabbzi/terra-dogagingproject-genomic-reports"])
  String? callback_docker_image_override
  String callback_docker_image = select_first([callback_docker_image_override, "tabbzi/terra-dogagingproject-genomic-reports-callback"])
  String? webhook_url_override
  String webhook_url = select_first([webhook_url_override, "${{ secrets.DAP_WEBHOOK_URL }}"])
  File python_SampleStatusJSON # bin/make_status_json.py

  call MakeJSON {
    input:
      dockerImage = workflow_docker_image,
      sampleID = sample_id,
      participantID = participant_id,
      gencoveID = gencove_id,
      datetime = date_time,
      status = sample_status,
      pythonSampleStatusJSON = python_SampleStatusJSON
  }

  call PostJSON {
    input:
      dockerImage = callback_docker_image,
      thisJSON = MakeJSON.sample_status_json,
      webhookURL = webhook_url
  }
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# TASKS

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

task MakeJSON {
  String dockerImage
  String sampleID
  String participantID
  String gencoveID
  String datetime
  String status
  File pythonSampleStatusJSON

  command <<<
    python3 ${pythonSampleStatusJSON} -B "${sampleID}" -P "${participantID}" -G "${gencoveID}" -D "${datetime}" -S "${status}"
  >>>

  output {
    File sample_status_json = '${sampleID}.JSON'
  }

  runtime {
    docker: dockerImage
    disks: "local-disk 5 HDD"
    memory: "2G"
  }
}

task PostJSON {
  String dockerImage
  String webhookURL
  File thisJSON

  command { /opt/callback.sh "${webhookURL}" "${thisJSON}"}

  output {
        #String? response_stdout = read_string(stdout())
        #String? response_stderr = read_string(stderr())
  }

  runtime {
    docker: dockerImage
    disks: "local-disk 5 HDD"
    memory: "2G"
  }

  parameter_meta {
    webhookURL: "The webhook URL to post JSON"
  }
}
