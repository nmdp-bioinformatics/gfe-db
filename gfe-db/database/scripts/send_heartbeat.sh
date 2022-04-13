#!/bin/bash

# This script sends heartbeats back to the StepFunctions API during the task execution.

while [ 1 ]
do
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Sending task heartbeat"
    aws stepfunctions send-task-heartbeat \
        --task-token "$TASK_TOKEN" \
        --region $REGION

    # Send TaskSuccess token to StepFunctions
    if [[ $? != "0" ]]; then
        exit 1
    fi

    sleep $HEARTBEAT_INTERVAL
done
