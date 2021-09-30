#!/usr/bin/env bash

# Set in root Makefile
STAGE=${1:-dev}
APP_NAME=${2:-gfe-db}
REGION=${3:-$(aws ec2 describe-availability-zones \
    --output text \
    --query 'AvailabilityZones[0].[RegionName]')}

setup_stack_name=$STAGE-$APP_NAME-setup
master_stack_name=$STAGE-$APP_NAME

echo "This will delete all data and resources in stacks:"
echo $setup_stack_name
echo $master_stack_name

read -p "Continue? " -n 1 -r
echo && echo

if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    echo "Exiting"
    [[ "$0" = "$BASH_SOURCE" ]] && exit 0 || return 0
fi

echo
echo -e "******* Deleting S3 objects... *******"
echo
# Get and delete objects on cfn bucket
DATA_BUCKET=$(aws ssm get-parameter \
	 		--name "/${APP_NAME}/${STAGE}/${REGION}/DataBucketName" | jq -r '.Parameter.Value')
aws s3 rm --recursive s3://$DATA_BUCKET
echo "Deleted bucket contents: $DATA_BUCKET"

echo
echo -e "******* Deleting ECR images... *******"
echo
# Get and delete objects on cfn bucket
# for service in build load database; do

# # Get the repo name from SSM
# ECR_REPOS=$(aws ssm get-parameter \
# 	 		--name "/${APP_NAME}/${STAGE}/${REGION}/DataBucketName" | jq -r '.Parameter.Value')

# TODO: add for loop for build, load and database
# Delete build service
build_repo=$(aws ssm get-parameter \
    --name "/${APP_NAME}/${STAGE}/${REGION}/BuildServiceRepositoryName" | jq -r '.Parameter.Value')

# Delete all images
for image_digest in $(aws ecr list-images --repository-name $build_repo | jq -r '.imageIds[].imageDigest'); do
    echo "Deleting ECR image: $image_digest"
    aws ecr batch-delete-image --repository-name $build_repo --image-ids imageDigest=$image_digest
done

# Delete load service
load_repo=$(aws ssm get-parameter \
    --name "/${APP_NAME}/${STAGE}/${REGION}/LoadServiceRepositoryName" | jq -r '.Parameter.Value')

# Delete all images
for image_digest in $(aws ecr list-images --repository-name $load_repo | jq -r '.imageIds[].imageDigest'); do
    echo "Deleting ECR image: $image_digest"
    aws ecr batch-delete-image --repository-name $load_repo --image-ids imageDigest=$image_digest
done

# Delete stacks
echo
echo -e "******* Deleting CloudFormation stacks... *******"
echo

# Delete CloudFormation nested stack
aws cloudformation delete-stack \
    --stack-name $master_stack_name 
echo "Deleted stack: $master_stack_name"

aws cloudformation delete-stack \
    --stack-name $setup_stack_name
echo "Deleted stack: $setup_stack_name"

echo "Finished"
exit 0