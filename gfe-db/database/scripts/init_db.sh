#!/bin/bash

CYPHER_PATH=$1

# gfe_bucket=$(aws ssm get-parameters \
#     --region $(curl -s http://169.254.169.254/latest/meta-data/placement/availability-zone | sed 's/\(.*\)[a-z]/\1/') \
#     --names "/gfe-db/dev/us-east-1/DataBucketName" \
#     | jq -r '.Parameters | map(select(.Version == 1))[0].Value')

# if [[ -z $gfe_bucket ]]; then
#     echo "S3 bucket not found."
#     exit 1
# else
#     echo "Found S3 bucket: $gfe_bucket"
# fi

# Run Cypher


# Test

# Exit