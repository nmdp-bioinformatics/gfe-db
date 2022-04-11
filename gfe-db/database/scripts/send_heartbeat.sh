#!/bin/bash

# This script watches a process and sends heartbeats back to the StepFunctions API during the task execution.
# When the process has exited the script will return either a success or failure signal, depending on the exit code.
# The task state will be reflected in the StepFunctions API and visible in the console.

# TODO: refactor so that this script only sends a heartbeat. It does not need to monitor the processed. 
# The load script will call this script as a chld process and kill it when it completes. the Load script
# will also handle returning the TaskFailure/Success tokens

# SELF_NAME=$(basename "$0")
# INTERVAL=$1

# # TODO: get name of state machine, tasktoken
# echo "Polling for new activities..."
# ACTIVITY=$(aws stepfunctions get-activity-task \
#     --activity-arn arn:aws:states:us-east-1:531868584498:activity:load-neo4j \
#     --worker-name gfe-db \
#     --region us-east-1)

# TASK_TOKEN=$(echo $ACTIVITY | jq -r '.taskToken')
# TASK_INPUT=$(echo $ACTIVITY | jq -r '.input')

echo "Activity found (send_heartbeat):"
echo "Task token: $TASK_TOKEN"
echo "Task input: $TASK_INPUT"

# heartbeat_interval=10
# counter=0

while [ 1 ]
do
    echo "Sending task heartbeat"
    aws stepfunctions send-task-heartbeat \
        --task-token "$TASK_TOKEN" \
        --region $REGION

    sleep $heartbeat_interval

    # PID=$(ps ax | grep "$COMMAND" | grep -v grep | grep -v "$SELF_NAME")
    # echo "$PID"
    # [ -n "$PID" ] && sleep $heartbeat_interval || break

    # counter=$((counter + 1))
    # [[ $counter -eq 2 ]] && break || true
done

# # TODO: check the exit code of the load script and update MESSAGE for success or failure
# message="complete"

# # Send TaskSuccess token to StepFunctions
# aws stepfunctions send-task-success \
#     --task-token "$TASK_TOKEN" \
#     --task-output "{\"message\":\"$message\"}" \
#     --region us-east-1

# echo "Process for command \"$COMMAND\" are no longer running" && exit 0