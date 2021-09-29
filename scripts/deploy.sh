#!/usr/bin/env bash

echo "Set AWS credentials and region:"
aws configure

# Set in root Makefile
STAGE=${1:-dev}
APP_NAME=${2:-gfe-db}
NEO4J_USERNAME=${3:-neo4j}
NEO4J_PASSWORD=${4:-gfedb}

CFN_DIR=cfn
CFN_OUTPUT_DIR=$CFN_DIR/output
REGION=$(aws ec2 describe-availability-zones \
    --output text \
    --query 'AvailabilityZones[0].[RegionName]')
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
EC2_KEY_PAIR=$STAGE-$APP_NAME-$REGION-ec2-key
DATA_BUCKET=$STAGE-$APP_NAME-$ACCOUNT_ID-$REGION
CFN_LOG_FILENAME=$STAGE-$APP_NAME-$(date +%s)

# Check if EC2 key pair exists for region: gfe-db-<region>, if not create one
CURRENT_EC2_KEY_PAIR=$(aws ec2 describe-key-pairs --key-name $EC2_KEY_PAIR | jq '.KeyPairs[0].KeyName')

if [ -z "$CURRENT_EC2_KEY_PAIR" ]; then
    echo "Creating EC2 key pair \"$EC2_KEY_PAIR\" ..."
    aws ec2 create-key-pair --key-name $EC2_KEY_PAIR | jq -r '.KeyMaterial' > $EC2_KEY_PAIR.pem
    sed -i '' "s/<ec2-key-pair-name>/$EC2_KEY_PAIR/" $CFN_DIR/setup.yml
else
    echo "Key pair found: $CURRENT_EC2_KEY_PAIR"
fi

echo "Deploying stacks..."

mkdir -p $CFN_OUTPUT_DIR/

# Deploy setup stack
setup_stack_name=$STAGE-$APP_NAME-setup
aws cloudformation deploy \
    --template-file $CFN_DIR/setup.yml \
    --stack-name $setup_stack_name	\
    --parameter-overrides \
        Stage=$STAGE \
        AppName=$APP_NAME \
        EC2KeyPairName=$EC2_KEY_PAIR \
        DataBucketName=$DATA_BUCKET \
        Neo4jUsername=$NEO4J_USERNAME \
        Neo4jPassword=$NEO4J_PASSWORD

aws cloudformation list-stack-resources --stack-name $setup_stack_name > $CFN_OUTPUT_DIR/$setup_stack_name.json

# Sync templates to S3
aws s3 cp --recursive $CFN_DIR/ s3://$DATA_BUCKET/templates/

# Deploy nested stacks
master_stack_name=$STAGE-$APP_NAME
aws cloudformation deploy \
    --template-file $CFN_DIR/master-stack.yml \
    --stack-name $master_stack_name \
    --capabilities "CAPABILITY_NAMED_IAM" \
    --parameter-overrides \
        Stage=$STAGE \
        AppName=$APP_NAME > master_stack_name.json \
        DataBucketName=$DATA_BUCKET

aws cloudformation list-stack-resources --stack-name $master_stack_name > $CFN_OUTPUT_DIR/$master_stack_name.json

# Describe all resources
for stack in $(aws cloudformation list-stacks \
    --output text \
    --query "StackSummaries[?contains(StackName, $master_stack_name) && (StackStatus==`CREATE_COMPLETE`||StackStatus==`UPDATE_COMPLETE`)].[StackName]") ; do 
    aws cloudformation describe-stack-resources \
        --stack-name $stack \
        --output json > $STAGE-$APP_NAME-$stack.json
     ; 
    done


echo "Finished"
exit 0