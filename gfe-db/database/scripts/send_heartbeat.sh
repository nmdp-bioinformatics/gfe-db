#!/bin/bash -x

# This script sends heartbeats back to the StepFunctions API during the task execution.

set -e

while true
do
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Sending task heartbeat"
    aws stepfunctions send-task-heartbeat \
        --task-token "$TASK_TOKEN" \
        --region $AWS_REGION

    # Send TaskSuccess token to StepFunctions
    if [[ $? != "0" ]]; then
        exit 1
    fi

    # TODO exit if StepFunctions returns activity timeout

    sleep $HEARTBEAT_INTERVAL
done
