#!/bin/bash

# This script checks for SQS messages containing environment variables. 
# If no messages are retrieved within the timeout period, the server shuts itself down.

BIN_DIR="./bin"
TIMEOUT=30
start=$SECONDS
end=$(($SECONDS+$TIMEOUT))

echo "Cloning most recent application version..."
git clone https://github.com/abk7777/gfe-db.git
cd gfe-db
git checkout fix/optimize-build 

echo "Creating Python virtual environment..."
python3.8 -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -r requirements.txt

echo "SQS Queue: $SQS_BUILD_QUEUE_URL"
echo "Start: $start"
echo "End: $end"

while [ $SECONDS -lt $end ]; do

    # Poll SQS every 5 seconds and save message
    echo "Checking for message..."
    aws sqs receive-message \
        --queue-url $SQS_BUILD_QUEUE_URL | \
            jq -r > sqs-message.json
    sleep 5

    # If the response has content, run the build
    if [ -s sqs-message.json ]; then
        echo "Message found."

        export RELEASES=$(cat sqs-message.json | jq '.Messages[].Body | fromjson | .release')
        export ALIGN=$(cat sqs-message.json | jq '.Messages[].Body | fromjson | .align')
        export KIR=$(cat sqs-message.json | jq '.Messages[].Body | fromjson | .kir')
        export MEM_PROFILE=$(cat sqs-message.json | jq '.Messages[].Body | fromjson | .mem_profile')
        export LIMIT=$(cat sqs-message.json | jq '.Messages[].Body | fromjson | .limit')

        echo "$(pwd)"
        echo "$BIN_DIR""/build.sh"
        sh "$BIN_DIR"/build.sh 100 # $LIMIT
    fi

done

echo "Operation timed out. Shutting down server."
exit 0
