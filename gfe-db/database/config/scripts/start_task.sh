#!/bin/bash -x

set -e

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

trap 'cause="Error on line $LINENO" && error=$? && send_result && kill 0' ERR

export AWS_REGION=$(curl --silent http://169.254.169.254/latest/dynamic/instance-identity/document | jq -r '.region')

# while loop using aws sqs receive-message to fetch a message from the queue, process it, and delete it
QUEUE_URL="https://sqs.us-east-1.amazonaws.com/531868584498/dev-gfe-db-pipeline-GfeDbLoadReleaseQueue-Lauwwh5Y8Uir"

while true; do

    message=$(aws sqs receive-message \
        --queue-url $QUEUE_URL \
        --max-number-of-messages 1 \
        --region us-east-1)

    if [[ -z $message ]]; then
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - No messages found"
        break

    else
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Message found"
        echo "$message" | jq -r

        export RECEIPT_HANDLE=$(echo "$message" | jq -r '.Messages[0].ReceiptHandle')
        export PARAMS=$(echo "$message" | jq -r '.Messages[0].Body')
        export TASK_TOKEN=$(echo "$PARAMS" | jq -r '.task_token')
        export RELEASE=$(echo "$PARAMS" | jq -r '.input.version')

        # Debug: create for loop to print RECEIPT_HANDLE, PARAMS, TASK_TOKEN, RELEASE
        for var in RECEIPT_HANDLE PARAMS TASK_TOKEN RELEASE; do
            echo "$var=${!var}"
        done

        echo "TASK_TOKEN=$TASK_TOKEN"
        echo "RELEASE=$RELEASE"

        # Check for release argument
        if [[ -z $RELEASE ]]; then
            echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Release version not found"
            kill -1 $$
        else
            echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Starting load process for $RELEASE"
        fi

        # TODO BOOKMARK
        # Run task - invoke load script 
        # TODO get s3 path from message
        bash load_db.sh $RELEASE
        TASK_EXIT_STATUS=$?
        echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Task exit status: $TASK_EXIT_STATUS"

        # Send TaskSuccess token to StepFunctions
        if [[ $TASK_EXIT_STATUS != "0" ]]; then
            status="FAILED"
            error="$TASK_EXIT_STATUS"
            cause="Error on line $LINENO"
            send_result

            # Return message to queue after failed processing
            aws sqs change-message-visibility \
                --queue-url $QUEUE_URL \
                --receipt-handle $RECEIPT_HANDLE \
                --visibility-timeout 0 \
                --region us-east-1
            echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Message returned to queue"

        else
            status="SUCCESS"
            send_result

            # Delete message from queue after successful processing
            aws sqs delete-message \
                --queue-url $QUEUE_URL \
                --receipt-handle $RECEIPT_HANDLE \
                --region us-east-1
            echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Message deleted"

        fi
    fi

### end while loop
done
exit 0
