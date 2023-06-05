#!/bin/bash -x

# check that APP_NAME and AWS_REGION are set from the environment
if [[ -z $APP_NAME ]]; then
    echo "APP_NAME environment variable not set"
    exit 1
fi

if [[ -z $AWS_REGION ]]; then
    echo "AWS_REGION environment variable not set"
    exit 1
fi

ACTIVITY_ARN=$(aws ssm get-parameter \
    --name "/${APP_NAME}/${STAGE}/${AWS_REGION}/LoadReleaseActivityArn" \
    --query "Parameter.Value" \
    --output text \
    --region "$AWS_REGION")

if [[ -z $ACTIVITY_ARN ]]; then
    echo "ACTIVITY_ARN environment variable not set"
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

trap 'cause="Error on line $LINENO" && error=$? && send_result && kill 0' ERR
while true; do

    # Poll StepFunctions API for new activities
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Polling for new activities..."
    export ACTIVITY=$(aws stepfunctions get-activity-task \
        --activity-arn "$ACTIVITY_ARN" \
        --worker-name "$APP_NAME" \
        --region "$AWS_REGION")

    if [[ -z $ACTIVITY ]]; then
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - No activities found"
        break

    else
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Activity found"
        echo "$ACTIVITY" | jq -r

        export TASK_TOKEN=$(echo "$ACTIVITY" | jq -r '.taskToken')
        export RELEASE=$(echo "$ACTIVITY" | jq -r '.input' | jq '.version')

        echo "TASK_TOKEN=$TASK_TOKEN"
        echo "RELEASE=$RELEASE"

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
