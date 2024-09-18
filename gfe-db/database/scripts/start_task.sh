#!/bin/bash -x

# Exit immediately if a command exits with a non-zero status
set -e

# source /home/ec2-user/env.sh

# # get APP_NAME and STAGE from arguments
# APP_NAME=$1
# STAGE=$2
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

ACTIVITY_ARN=$(aws ssm get-parameter \
    --name "/${APP_NAME}/${STAGE}/${AWS_REGION}/LoadReleaseActivityArn" \
    --query "Parameter.Value" \
    --output text \
    --region "$AWS_REGION")

if [[ -z $ACTIVITY_ARN ]]; then
    ERR_MSG="ACTIVITY_ARN environment variable not set"
    echo $ERR_MSG >&2
    exit 1
fi

# Send task failure if script errors
send_result () {
    if [[ $status = "SUCCESS" ]]; then
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Sending task success"
        aws stepfunctions send-task-success \
            --task-token "$TASK_TOKEN" \
            --task-output "{\"status\":\"$status\"}" \
            --region "$AWS_REGION"
    else
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Sending task failure"
        aws stepfunctions send-task-failure \
            --task-token "$TASK_TOKEN" \
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

    # Poll StepFunctions API for new activities
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Polling for new activities..."
    export ACTIVITY=$(aws stepfunctions get-activity-task \
        --activity-arn "$ACTIVITY_ARN" \
        --worker-name "$APP_NAME" \
        --region "$AWS_REGION" \
        --cli-connect-timeout 65)

    if [[ -z $ACTIVITY ]]; then
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - No activities found"
        break

    else
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Activity found"
        echo "$ACTIVITY"

        # Check for task token
        export TASK_TOKEN=$(echo "$ACTIVITY" | jq -r '.taskToken')
        if [[ -z $TASK_TOKEN ]]; then
            ERR_MSG="TASK_TOKEN not found"
            echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - $ERR_MSG" >&2
            exit 1
        else
            echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Found task token"
        fi

        # Check for release argument
        release=$(echo "$ACTIVITY" | jq -r '.input' | jq -r '.input.version')
        if [[ -z $release || "$release" == "null" || ! $release =~ ^[0-9]{1,4}$ ]]; then
            ERR_MSG="Release version not found, is 'null', or is not a 1-4 digit integer"
            echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - $ERR_MSG" >&2
            exit 1
        else
            echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Starting load process for $release"
        fi

        export HEARTBEAT_INTERVAL=30
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Begin task heartbeat"
        bash send_heartbeat.sh &
        send_heartbeat_pid=$!

        # Run task - invoke load script 
        # TODO get s3 path from step functions payload
        bash load_db.sh $release
        TASK_EXIT_STATUS=$?
        
        # Send TaskSuccess token to StepFunctions
        if [[ $TASK_EXIT_STATUS != "0" ]]; then
            status="FAILED"
            error="Non-zero exit code: $TASK_EXIT_STATUS"
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