#!/bin/bash

set -e

# Send task failure if script errors
send_result () {
    if [[ $status = "SUCCESS" ]]; then
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Sending task success"
        aws stepfunctions send-task-success \
            --task-token "$TASK_TOKEN" \
            --task-output "{\"status\":\"$status\"}" \
            --region $REGION
    else
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Sending task failure"
        aws stepfunctions send-task-failure \
            --task-token "$TASK_TOKEN" \
            --cause "$cause" \
            --error "$error" \
            --region $REGION
    fi
}

trap 'cause="Error on line $LINENO" && error=$? && send_result' ERR
trap 'kill 0' EXIT

export REGION=$(curl --silent http://169.254.169.254/latest/dynamic/instance-identity/document | jq -r .region)

export RELEASE=$1

# Check for release argument
if [[ -z $RELEASE ]]; then
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Release version not found"
    exit 1
else
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Starting load process for $RELEASE"
fi

# Poll StepFunctions API for new activities
echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Polling for new activities..."
export ACTIVITY=$(aws stepfunctions get-activity-task \
    --activity-arn arn:aws:states:us-east-1:531868584498:activity:load-neo4j \
    --worker-name gfe-db \
    --region us-east-1)

echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Activity found"
# echo $ACTIVITY | jq -r

export TASK_TOKEN=$(echo $ACTIVITY | jq -r '.taskToken')
export TASK_INPUT=$(echo $ACTIVITY | jq -r '.input')

# TODO: parameterize heartbeat and set interval / 2
export HEARTBEAT_INTERVAL=10
bash send_heartbeat.sh &
send_heartbeat_pid=$!

# Invoke load script 
bash load_db.sh $RELEASE
TASK_EXIT_STATUS=$?
echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Task exit status: $TASK_EXIT_STATUS"

# Send TaskSuccess token to StepFunctions
if [[ $TASK_EXIT_STATUS != "0" ]]; then
    status="FAILED"
    error="$TASK_EXIT_STATUS"
    cause="Error on line $LINENO"
    send_result 
else
    status="SUCCESS"
    send_result
    kill $send_heartbeat_pid
fi

exit 0
