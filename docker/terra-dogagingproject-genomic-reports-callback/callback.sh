#!/bin/bash

WEBHOOK_URL="${1}"
INPUT_JSON="${2}"

postToWebhook() {
    WEBHOOK_URL="${1}"
    INPUT_JSON="${2}"

    # --silent - don't show download progress
    # --write-out - Print HTTP status code
    # --output - sending HTML output to '/dev/null'
    # --show-error - Show error if failure

    curl --header "Content-Type: application/json" \
        --silent --show-error \
        --write-out "%{http_code}" --output /dev/null \
        --request POST --data @${INPUT_JSON} \
        ${WEBHOOK_URL}
}

echo "Preparing to call DAP webhook"
http_code=$(postToWebhook ${WEBHOOK_URL} ${INPUT_JSON})
status=$?

if [ $status -ne 0 ]; then
    echo "Command failed to POST to webhook! Respose code $http_code" 1>&2
    exit 1
fi
echo "Call to DAP webhook yielded HTTP code $http_code"
exit 0