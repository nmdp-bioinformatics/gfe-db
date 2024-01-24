#!/bin/bash -x

# Exit immediately if a command exits with a non-zero status
set -e

# get APP_NAME and STAGE from arguments
APP_NAME=$1
STAGE=$2
ERR_MSG=null

if [[ -z $APP_NAME ]]; then
    ERR_MSG="APP_NAME environment variable not set"
    echo $ERR_MSG >&2
    exit 1
fi

if [[ -z $STAGE ]]; then
    ERR_MSG="STAGE environment variable not set"
    echo $ERR_MSG >&2
    exit 1
fi

# check that AWS_REGION is set from the environment
if [[ -z $AWS_REGION ]]; then
    export AWS_REGION=$(curl --silent http://169.254.169.254/latest/dynamic/instance-identity/document | jq -r '.region')
fi

QUEUE_URL=$(aws ssm get-parameter \
    --name "/${APP_NAME}/${STAGE}/${AWS_REGION}/GfeDbLoadQueueUrl" \
    --query "Parameter.Value" \
    --output text \
    --region "$AWS_REGION")

if [[ -z $QUEUE_URL ]]; then
    ERR_MSG="Queue URL not found"
    echo $ERR_MSG >&2
    exit 1
fi

# Send task failure if script errors
send_result () {
    if [[ $status = "SUCCESS" ]]; then
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Sending task success"
        aws stepfunctions send-task-success \
            --task-token "$task_token" \
            --task-output "{\"status\":\"$status\"}" \
            --region "$AWS_REGION"
    else
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Sending task failure"
        aws stepfunctions send-task-failure \
            --task-token "$task_token" \
            --cause "$cause" \
            --error "$error" \
            --region "$AWS_REGION"
    fi
}

# This will send task failure if any error is encountered, but sometimes errors can occur that do not affect that actual loading process
trap 'cause="Error on line $LINENO: $ERR_MSG" && error=$? && send_result && kill 0' ERR

# start_time=$(date +%s)
# timeout=120
while true; do

    # Poll SQS for new messages
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Polling for new messages..."
    export message=$(aws sqs receive-message \
        --queue-url "$QUEUE_URL" \
        --region "$AWS_REGION" \
        --wait-time-seconds 20 \
        --max-number-of-messages 1 \
        | jq -r '.Messages[0]')

    if [[ -z $message ]]; then
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - No messages found"
        break
    else
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Message found"

        export task_token=$(echo "$message | jq -r '.Body.TaskToken'")
        export task_input=$(echo "$message | jq -r '.Body.Input'")

        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - task_token=$task_token"
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - task_input=$task_input"

        # Check for release argument
        release=$(echo $task_input | jq -r '.version')
        if [[ -z $task_input || "$task_input" == "" ]]; then
            ERR_MSG="Release version not found, or is empty."
            echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - $ERR_MSG" >&2
            status="FAILED"
            error="1"
            cause="$ERR_MSG"
            send_result
            kill 0
        else
            echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Loading release version $release for commit $(echo $task_input | jq -r '.commit_sha')"
        fi

        # export HEARTBEAT_INTERVAL=30
        # bash send_heartbeat.sh &
        # send_heartbeat_pid=$!

        # Run task - invoke load script 
        # TODO get s3 path from step functions payload
        bash load_db.sh $RELEASE
        TASK_EXIT_STATUS=$?
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Task exit status: $TASK_EXIT_STATUS"

        # Send TaskSuccess token to StepFunctions
        if [[ $TASK_EXIT_STATUS != "0" ]]; then
            status="FAILED"
            error="$TASK_EXIT_STATUS"
            cause="Error on line $LINENO"
            send_result
            kill 0

        else
            status="SUCCESS"
            send_result
            kill $send_heartbeat_pid
        fi
    fi
done

exit 0
