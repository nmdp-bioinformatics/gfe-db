#!/bin/bash -x

# Exit immediately if a command exits with a non-zero status
set -e

ERR_MSG=null

source /home/ec2-user/env.sh

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

# Check for load event argument from command line
if [[ -z $1 ]]; then
    ERR_MSG="No load event provided"
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - $ERR_MSG" >&2
    exit 1
fi

load_event=$1
# Log the load event
echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Load event: $load_event"

release=$(echo "$load_event" | jq -r '.sqs.Body.input.version')
if [[ -z $release || "$release" == "null" || ! $release =~ ^[0-9]{1,4}$ ]]; then
    ERR_MSG="Release version \"$release\" not found, is \"null\", or is not a 1-4 digit integer"
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - $ERR_MSG" >&2
    exit 1
else
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Starting load process for $release"
fi

# Run task - invoke load script 
bash load_db.sh $release
TASK_EXIT_STATUS=$?

# Exit with the status of the load_db.sh script
exit $TASK_EXIT_STATUS