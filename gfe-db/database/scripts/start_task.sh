#!/bin/bash -x

set -e

# TODO remove the application logic to make this script agnostic

# Send task failure if script errors
send_result () {
    if [[ $status = "SUCCESS" ]]; then
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Sending task success"
        aws stepfunctions send-task-success \
            --task-token "$TASK_TOKEN" \
            --task-output "{\"status\":\"$status\"}" \
            --region $AWS_REGION
    else
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Sending task failure"
        aws stepfunctions send-task-failure \
            --task-token "$TASK_TOKEN" \
            --cause "$cause" \
            --error "$error" \
            --region $AWS_REGION
    fi
}

trap 'cause="Script failed due to error on line $LINENO. Please see logs in System Manager Run Command history for more details" && error=$? && send_result && kill 0' ERR

export AWS_REGION=$(curl --silent http://169.254.169.254/latest/dynamic/instance-identity/document | jq -r '.region')

export PARAMS=$1
if [[ -z $PARAMS ]]; then
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - No parameters found"
    exit 1
else
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Found parameters:"
    echo "$PARAMS"
fi

export ACTIVITY_ARN=$(echo $PARAMS | jq -r '.params.activity_arn')

# TODO can now source APP_NAME and STAGE from env.sh
export APP_NAME=$(echo $PARAMS | jq -r '.params.app_name')

echo "ACTIVITY_ARN=$ACTIVITY_ARN"
echo "APP_NAME=$APP_NAME"

# Poll StepFunctions API for new activities
echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Polling for new activities..."
export ACTIVITY=$(aws stepfunctions get-activity-task \
    --activity-arn $ACTIVITY_ARN \
    --worker-name $APP_NAME \
    --region $AWS_REGION)

echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Activity found"

export TASK_TOKEN=$(echo $ACTIVITY | jq -r '.taskToken')
export TASK_INPUT=$(echo $ACTIVITY | jq -r '.input')

echo "TASK_TOKEN=$TASK_TOKEN"
echo "TASK_INPUT=$TASK_INPUT"

export RELEASE=$(echo $TASK_INPUT | jq -r '.RELEASES')
export ALIGN=$(echo $TASK_INPUT | jq -r '.ALIGN')
export KIR=$(echo $TASK_INPUT | jq -r '.KIR')

echo "RELEASE=$RELEASE"
echo "ALIGN=$ALIGN"
echo "KIR=$KIR"

# Check for release argument
if [[ -z $RELEASE ]]; then
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Release version not found"
    kill -1 $$
else
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Starting load process for $RELEASE"
fi

# TODO: parameterize heartbeat and set interval / 2
export HEARTBEAT_INTERVAL=30
bash send_heartbeat.sh &
send_heartbeat_pid=$!

# Run task - invoke load script 
bash load_db.sh $RELEASE
TASK_EXIT_STATUS=$?
echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Task exit status: $TASK_EXIT_STATUS"

# Send TaskSuccess token to StepFunctions
if [[ $TASK_EXIT_STATUS != "0" ]]; then
    status="FAILED"
    error="$TASK_EXIT_STATUS"
    cause="Task failed due to error on line $LINENO. Please see logs in System Manager Run Command history for more details."
    send_result
    kill 0
else
    status="SUCCESS"
    send_result
    kill $send_heartbeat_pid
fi

exit 0
