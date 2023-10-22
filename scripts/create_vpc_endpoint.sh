#!/bin/bash -x

set -e

# Deploy a VPC endpoint for either SSM, Secrets Manager, or both
services=$1

# throw error if services is empty, or if anything other than `ssm` or `secretsmanager` is passed
if [[ -z $services ]]; then
    echo "No services specified"
    exit 1
elif [[ $services != "ssm" ]] && [[ $services != "secretsmanager" ]] && [[ $services != "ssm secretsmanager" ]] && [[ $services != "secretsmanager ssm" ]]; then
    echo "Invalid service specified. Valid services are 'ssm', 'secretsmanager', or 'ssm secretsmanager'"
    exit 1
fi

# get the Neo4j sg ID from SSM Parameter store
neo4j_sg_id=$(aws ssm get-parameter \
    --name /${APP_NAME}/${STAGE}/${AWS_REGION}/Neo4jDatabaseSecurityGroupId \
    --query Parameter.Value \
    --output text)

# if neo4j_sg_id is empty, exit 1
if [[ -z $neo4j_sg_id ]]; then
    echo "Neo4j sg ID not found in SSM Parameter store"
    exit 1
else
    echo "Found Neo4j security group id: $neo4j_sg_id"
fi

# iterate over services separated by space and run `aws ec2 create-vpc-endpoint --dry-run --vpc-id ${VPC_ID} --service-name com.amazonaws.${AWS_REGION}.$service --vpc-endpoint-type Interface` for each
for service in $services; do
    echo "Creating VPC endpoint for $service..."
    res=$(aws ec2 create-vpc-endpoint \
        --vpc-id ${VPC_ID} \
        --service-name com.amazonaws.${AWS_REGION}.$service \
        --security-group-ids $neo4j_sg_id \
        --vpc-endpoint-type Interface)
    echo $res | jq -r
done
