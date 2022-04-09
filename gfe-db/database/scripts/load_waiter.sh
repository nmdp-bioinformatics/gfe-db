#!/bin/bash

# TODO: get name of state machine, tasktoken, heartbeatseconds etc

SELF_NAME=$(basename "$0")
COMMAND=$1

while [ 1 ]
do
    # TODO: Send HeartBeat token every 10 seconds (60 in prod)
    sleep 1
    PID=$(ps ax | grep "$COMMAND" | grep -v grep | grep -v "$SELF_NAME")
    echo $PID
    [ -n "$PID" ] && sleep 5 || break
done

# TODO send TaskCompleted token to StepFunctions
echo "Process for command \"$COMMAND\" are no longer running" && exit 0