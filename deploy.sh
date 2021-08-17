#!/usr/bin/env bash

echo "Set AWS credentials and region:"
aws configure

APP_NAME=gfe-db
CFN_DIR=cfn
REGION=$(aws ec2 describe-availability-zones \
    --output text \
    --query 'AvailabilityZones[0].[RegionName]'
)
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)

# Check if EC2 key pair exists for region: gfe-db-<region>, if not create one
EC2_KEY_PAIR=$(aws ec2 describe-key-pairs --key-name $APP_NAME-$REGION | jq '.KeyPairs[0].KeyName')

if [ -z "$EC2_KEY_PAIR" ]; then
    echo "Creating EC2 key pair \"$APP_NAME-$REGION\" ..."
    aws ec2 create-key-pair --key-name $APP_NAME-$REGION | jq -r '.KeyMaterial' > $APP_NAME-$REGION.pem
    sed -i '' "s/<ec2-key-pair-name>/$APP_NAME-$REGION/" $CFN_DIR/setup.yml
else
    echo "Key pair found: $EC2_KEY_PAIR"
fi

echo "Deploying stacks..."

# Deploy setup stack
aws cloudformation deploy \
   --template-file $CFN_DIR/setup.yml \
   --stack-name $APP_NAME-setup

# Sync templates to S3
aws s3 cp --recursive $CFN_DIR/ s3://$APP_NAME-$ACCOUNT_ID-$REGION/templates/

# Deploy nested stacks
aws cloudformation deploy \
    --template-file $CFN_DIR/master-stack.yml \
    --stack-name $APP_NAME \
    --capabilities "CAPABILITY_NAMED_IAM"

echo "Finished"
exit 0