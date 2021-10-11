#!/usr/bin/env bash

echo "Set AWS credentials and region:"
aws configure

# Variables set in root Makefile
STAGE=${1:-dev}
APP_NAME=${2:-gfe-db}
REGION=${3:-$(aws ec2 describe-availability-zones \
    --output text \
    --query 'AvailabilityZones[0].[RegionName]')}

# Deployment variables
export ROOT=$(dirname $(dirname "$0"))
export BUILD_SERVICE_DIR=$ROOT/build
export LOAD_SERVICE_DIR=$ROOT/load

NEO4J_USERNAME=${4:-neo4j}
NEO4J_PASSWORD=${5:-gfedb}
CFN_DIR=cfn
# CFN_OUTPUT_DIR=$CFN_DIR/output
CFN_LOG_FILENAME=$STAGE-$APP_NAME-$(date +%s)
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
EC2_KEY_PAIR=$STAGE-$APP_NAME-$REGION-ec2-key
DATA_BUCKET=$STAGE-$APP_NAME-$ACCOUNT_ID-$REGION

# Check if EC2 key pair exists for region: gfe-db-<region>, if not create one
CURRENT_EC2_KEY_PAIR=$(aws ec2 describe-key-pairs --key-name $EC2_KEY_PAIR | jq '.KeyPairs[0].KeyName')

if [ -z "$CURRENT_EC2_KEY_PAIR" ]; then
    echo "Creating EC2 key pair \"$EC2_KEY_PAIR\" ..."
    aws ec2 create-key-pair --key-name $EC2_KEY_PAIR | jq -r '.KeyMaterial' > $EC2_KEY_PAIR.pem
    sed -i '' "s/<ec2-key-pair-name>/$EC2_KEY_PAIR/" $CFN_DIR/setup.yml
else
    echo "Key pair found: $CURRENT_EC2_KEY_PAIR"
fi

# Create CloudFormation stacks
echo "Deploying stacks..."

# mkdir -p $CFN_OUTPUT_DIR/

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

# Sync templates to S3
aws s3 cp --recursive $CFN_DIR/ s3://$DATA_BUCKET/templates/

# Deploy nested stacks (database, ECR repos, Batch environment and StepFunctions)
master_stack_name=$STAGE-$APP_NAME
aws cloudformation deploy \
    --template-file $CFN_DIR/master-stack.yml \
    --stack-name $master_stack_name \
    --capabilities "CAPABILITY_NAMED_IAM" \
    --parameter-overrides \
        Stage=$STAGE \
        AppName=$APP_NAME \
        DataBucketName=$DATA_BUCKET

# # Login to docker/ECR
# aws ecr get-login-password --region $REGION | docker login --username AWS --password-stdin $ACCOUNT_ID.dkr.ecr.$REGION.amazonaws.com

# # Deploy build service image to ECR
# echo "Deploying container image: build service"
# docker build -t $STAGE-$APP_NAME-build-service $BUILD_SERVICE_DIR/
# docker tag $STAGE-$APP_NAME-build-service:latest $ACCOUNT_ID.dkr.ecr.$REGION.amazonaws.com/$STAGE-$APP_NAME-build-service:latest
# docker push $ACCOUNT_ID.dkr.ecr.$REGION.amazonaws.com/$STAGE-$APP_NAME-build-service:latest

# # Deploy load service image to ECR
# echo "Deploying container image: load service"
# docker build -t $STAGE-$APP_NAME-load-service $LOAD_SERVICE_DIR
# docker tag $STAGE-$APP_NAME-load-service:latest $ACCOUNT_ID.dkr.ecr.$REGION.amazonaws.com/$STAGE-$APP_NAME-load-service:latest
# docker push $ACCOUNT_ID.dkr.ecr.$REGION.amazonaws.com/$STAGE-$APP_NAME-load-service:latest

echo "Neo4j can be accessed at this URL: $(aws ssm get-parameter --name "/$APP_NAME/$STAGE/$REGION/Neo4jDatabaseEndpoint" | jq -r '.Parameter.Value'):7474"

exit 0